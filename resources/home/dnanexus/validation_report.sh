#!/bin/bash
# set -e
set -x
#####################################################################
# VERSION
#####################################################################
VERSION="20150622-1"

#####################################################################
# SOFTWARE AND DATA LOCATIONS
#####################################################################

#
# Software & Resources
#
RSCRIPT="/home/dnanexus/bin/Rscript"
JAVA=`which java`
GIT=`which git`
BEDTOOLS=`which bedtools`
TABIX=`which tabix`
BGZIP=`which bgzip`
RESOURCES_HEAD="/home/dnanexus/resources"
SCRATCH_DEFAULT="/tmp"

RTG_CORE="${RESOURCES_HEAD}/rtg-core/rtg-core.jar"
RTG_THREADS=`nproc`

# Calculate 80% of memory size, for java
mem_in_mb=`head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}'`
RTG_VCFEVAL="${JAVA} -Xmx${mem_in_mb}m -jar ${RTG_CORE} vcfeval -T ${RTG_THREADS}"


#
# Data
#
GOLD_CALLS_VCFGZ="${RESOURCES_HEAD}/gold_standard/calls-2.19.vcf.gz"
GOLD_CALLS_VCFGZTBI="${RESOURCES_HEAD}/gold_standard/calls-2.19.vcf.gz.tbi"
GOLD_HARDMASK_VALID_REGIONS_BEDGZ="${RESOURCES_HEAD}/gold_standard/valid_regions-2.19.bed.gz"
REFERENCE_SDF="${RESOURCES_HEAD}/reference/ref.sdf/"
FUNCTIONAL_REGIONS_BEDGZ_PREFIX="${RESOURCES_HEAD}/functional_regions/"
MASK_REGIONS_BEDGZ_PREFIX="${RESOURCES_HEAD}/mask_regions/"
REFERENCE_BSGENOME="BSgenome.HSapiens.1000g.37d5"		# This is a custom package, available at /share/ClusterShare/biodata/contrib/marpin/reference/hs37d5/build/BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz


EXEC_DIR=$(pwd)


SCRATCH=$(mktemp -d --tmpdir=${SCRATCH_DEFAULT} valrept.XXXXXXXXXX)
INPUT_SCRATCH="${SCRATCH}/input"
RTG_OVERLAP_SCRATCH="${SCRATCH}/overlap"
KNITR_SCRATCH="${SCRATCH}/knitr"


#####################################################################
# COMMAND LINE PARSING
#####################################################################
print_usage() {
cat << EOF
Usage: ${0##*/} [-o OUTFILE] [-r BEDFILE] [-x] <INFILE>

Create a WGS validation report.

    INFILE       Input NA12878 genotype calls, in vcf.gz format.
    -o OUTFILE   Write the report to OUTFILE (default: report.pdf)
    -r BEDFILE   Restrict analysis to the regions in BEDFILE only.
                 Default: the full genome is considered.
    -x           Generate an extended report, with threshold and
                 score diagnostics appended to the standard report.
                 Default: generate the standard report only.
    -h           Display this help and exit.

Version ${VERSION}

Mark Pinese
EOF
}

# Head location for resources bundle
OPTIND=1
input_vcfgz_path=""
region_bed_supplied=0
region_bed_path="NA"
output_pdf_path="${EXEC_DIR}/report.pdf"
extended=0

while getopts "r:o:hx" opt; do
	case "$opt" in
		h)
			print_usage
			exit 0
			;;
		o)
			output_pdf_path=$(readlink -f $OPTARG)
			;;
		r)
			region_bed_supplied=1
			region_bed_path=$(readlink -f $OPTARG)
			;;
		x)
			extended=1
			;;
		'?')
			print_usage >&2
			exit 1
			;;
	esac
done

shift $((OPTIND-1))

if [ $# -ne 1 ]; then
	print_usage >&2
	exit 1
fi

input_vcfgz_path=$1


#####################################################################
# PARAMETER CHECKING
#####################################################################

if [ ! -d "${RESOURCES_HEAD}" ]; then
  echo >&2 "Error: Resources path ${RESOURCES_HEAD} not found."
  exit 2
fi

if [ ! -e ${input_vcfgz_path} ]; then
	echo >&2 "Error: Input file ${input_vcfgz_path} not found."
	exit 3
fi

if [ -e ${output_pdf_path} ]; then
	echo >&2 "Error: Output file ${output_pdf_path} already exists."
	exit 4
fi

if [ ${region_bed_supplied} -eq 1 ] && [ ! -e ${region_bed_path} ]; then
  echo >&2 "Error: Region file ${region_bed_path} not found."
  exit 5
fi


#####################################################################
# SOFTWARE CHECKING
#####################################################################
if [ ! -e ${RSCRIPT} ]; then
  echo >&2 "Error: Rscript executable not found."
  exit 6
fi

if [ ! -e ${JAVA} ]; then
  echo >&2 "Error: java executable not found."
  exit 7
fi

if [ ! -e ${GIT} ]; then
  echo >&2 "Error: git executable not found."
  exit 8
fi

#set -e should catch this.
[[ -f report.Rnw ]] && [[ -f report_functions.R ]] && [[ -f report_debug.Rnw ]] && [[ -f report_calculations.R ]]


#####################################################################
# R PACKAGE CHECKING
#####################################################################
R --vanilla -e "suppressWarnings(require(\"${REFERENCE_BSGENOME}\", quiet=TRUE)) || stop(\"${REFERENCE_BSGENOME} package not installed\")"


#####################################################################
# VERSIONING
#####################################################################
VERSION_EXEC_HOST=$(uname -a)


#####################################################################
# CREATE TEMPORARY DIRECTORIES
#####################################################################
mkdir -p ${SCRATCH}
mkdir -p ${KNITR_SCRATCH}
mkdir -p ${INPUT_SCRATCH}


#####################################################################
# SUBSET TO REGION BED
#####################################################################
if [ ${region_bed_supplied} -eq 1 ]; then
  echo "Subsetting input files to supplied BED..."
  # Sort the region bed
  sort -k1,1 -k2,2n ${region_bed_path} > ${INPUT_SCRATCH}/region.bed

  # Perform the intersection
  ${BEDTOOLS} intersect -wa -sorted -header -a ${input_vcfgz_path} -b ${INPUT_SCRATCH}/region.bed | ${BGZIP} > ${INPUT_SCRATCH}/test_variants.vcf.gz
  ${BEDTOOLS} intersect -wa -sorted -header -a ${GOLD_CALLS_VCFGZ} -b ${INPUT_SCRATCH}/region.bed | ${BGZIP} > ${INPUT_SCRATCH}/gold_variants.vcf.gz
  ${BEDTOOLS} intersect -wa -sorted -a ${GOLD_HARDMASK_VALID_REGIONS_BEDGZ} -b ${INPUT_SCRATCH}/region.bed | ${BGZIP} > ${INPUT_SCRATCH}/gold_regions.bed.gz

  # We need to re-index the gold variants
  ${TABIX} -p vcf ${INPUT_SCRATCH}/gold_variants.vcf.gz
else
  # No region bed supplied; copy over the files in their entirety
  cp ${input_vcfgz_path} ${INPUT_SCRATCH}/test_variants.vcf.gz
  cp ${GOLD_CALLS_VCFGZ} ${INPUT_SCRATCH}/gold_variants.vcf.gz
  cp ${GOLD_HARDMASK_VALID_REGIONS_BEDGZ} ${INPUT_SCRATCH}/gold_regions.bed.gz

  # We can use the gold standard index unchanged, so no need to 
  # index it as above.
  cp ${GOLD_CALLS_VCFGZTBI} ${INPUT_SCRATCH}/gold_variants.vcf.gz.tbi
fi

# Regardless of whether a region bed was supplied or not, we still
# need to index the test variant vcf.
${TABIX} -p vcf ${INPUT_SCRATCH}/test_variants.vcf.gz



#####################################################################
# VCF OVERLAP EVALUATION
#####################################################################
echo "Computing VCF overlaps..."

if [ -e ${RTG_OVERLAP_SCRATCH} ]; then
	echo "Overlap scratch directory ${RTG_OVERLAP_SCRATCH} already exists.  Clearing scratch directory and continuing..."
	rm -rf ${RTG_OVERLAP_SCRATCH}
fi

${RTG_VCFEVAL} --all-records -b ${INPUT_SCRATCH}/gold_variants.vcf.gz -c ${INPUT_SCRATCH}/test_variants.vcf.gz -t ${REFERENCE_SDF} -o ${RTG_OVERLAP_SCRATCH} > /dev/null 2>&1

# TODO: Parse ${RTG_OVERLAP_SCRATCH}/vcfeval.log to identify regions to exclude
# eg Evaluation too complex (5001 unresolved paths, 18033 iterations) at reference region 2:105849275-105849281. Variants in this region will not be included in results.
# This will be required to remove the TNs in this region, but it's polish.


#####################################################################
# CALCULATIONS FOR REPORT
#####################################################################
echo "Performing calculations for report..."

# knitr doesn't play well with building knits outside of its working
# directory.  Currently we get around this with a bit of a kludge, 
# by copying the report files to ${KNITR_SCRATCH}, then executing in 
# that directory.
cp -f report.Rnw ${KNITR_SCRATCH}
cp -f report_functions.R ${KNITR_SCRATCH}
cp -f report_extended.Rnw ${KNITR_SCRATCH}
cp -f report_calculations.R ${KNITR_SCRATCH}
cd ${KNITR_SCRATCH}

# Run the script
# SECURITY WARNING: Code injection possible in VERSION_ variables.
# Ensure that git branch and execution host names can not be under
# malicious control.
${RSCRIPT} --vanilla report_calculations.R ${extended} ${input_vcfgz_path} \
  ${region_bed_supplied} ${region_bed_path} \
  ${RTG_OVERLAP_SCRATCH}/tp.vcf.gz ${RTG_OVERLAP_SCRATCH}/fp.vcf.gz ${RTG_OVERLAP_SCRATCH}/fn.vcf.gz \
  ${GOLD_CALLS_VCFGZ} ${REFERENCE_BSGENOME} ${INPUT_SCRATCH}/gold_regions.bed.gz \
  ${FUNCTIONAL_REGIONS_BEDGZ_PREFIX} ${MASK_REGIONS_BEDGZ_PREFIX} \
   "'"${VERSION}"'" "${VERSION_EXEC_HOST}"

# NB above: old ${GOLD_CALLS_VCFGZ} kept as it's not actually used in the R script, and gives the full path.
# Regions subset though, as it *is* actually used in R.
# TODO for later: Move those cmd line params to a KV file.  In doing so, give both original (for rept), and subset (for actual use) paths


#####################################################################
# REPORT GENERATION
#####################################################################
echo "Generating report..."

${RSCRIPT} --vanilla -e "library(knitr); knit('report.Rnw', output = 'report.tex')"

# Latex often 'fails' (returns a nonzero exit status), but still 
# generates a report.  Keep going when this happens, and test for
# failure explicitly later.
set +e

# Remove the report.pdf that may be present in the scratch directory,
# so we can later check whether pdflatex successfully built a report
# or not.
rm -f ${KNITR_SCRATCH}/report.pdf

# Run pdflatex
pdflatex -interaction nonstopmode report.tex
pdflatex -interaction nonstopmode report.tex

# Check  whether the report.pdf was generated
if [ ! -e ${KNITR_SCRATCH}/report.pdf ]; then
	echo >&2 "Error: pdflatex did not successfully generate report.pdf."
	echo >&2 "Check ${KNITR_SCRATCH}/report.tex and the latex log for errors."
	exit 9
fi

# Copy the completed report to the final destination
cp "${KNITR_SCRATCH}/report.pdf" "${output_pdf_path}"

echo "Report generated successfully."

echo "Done."
