#!/bin/bash
#set -e -x -o pipefail
set -x -o pipefail

#####################################################################
# VERSION
#####################################################################
export CONST_VERSION_SCRIPT="20150624-1"


#####################################################################
# SOFTWARE AND DATA LOCATIONS
#####################################################################

# dnanexus jobs run as root; assume that if the user is root, we're
# on dx.
export IS_DNANEXUS=0
if [ `whoami` == root ]; then
  IS_DNANEXUS=1
fi

# Software & Resources
if [ ${IS_DNANEXUS} -eq 1 ]; then
  PATH_RESOURCES_HEAD="/home/dnanexus/resources"
  PATH_SCRATCH_DEFAULT="/tmp"

  RSCRIPT="/home/dnanexus/bin/Rscript"
  R="/home/dnanexus/bin/R"
  JAVA=`which java`
  BEDTOOLS=`which bedtools`
  TABIX=`which tabix`
  BGZIP=`which bgzip`

  RTG_CORE="${PATH_RESOURCES_HEAD}/rtg-core/rtg-core.jar"
  RTG_THREADS=`nproc`
  mem_in_mb=`head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}'`                 # Calculate 80% of memory size, for java
  RTG_VCFEVAL="${JAVA} -Xmx${mem_in_mb}m -jar ${RTG_CORE} vcfeval -T ${RTG_THREADS}"
else
  PATH_RESOURCES_HEAD="/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources"
  PATH_SCRATCH_DEFAULT="/directflow/ClinicalGenomicsPipeline/tmp"

  RSCRIPT="/home/marpin/bin/Rscript"
  R="/home/marpin/bin/R"
  JAVA="/usr/java/latest/bin/java"
  BEDTOOLS="/home/marpin/software/bedtools2/bin/bedtools"
  TABIX="/home/marpin/software/htslib/tabix"
  BGZIP="/home/marpin/software/htslib/bgzip"

  RTG_CORE="${PATH_RESOURCES_HEAD}/rtg-core/rtg-core.jar"
  RTG_THREADS=4
  RTG_VCFEVAL="${JAVA} -Xmx8G -jar ${RTG_CORE} vcfeval -T ${RTG_THREADS}"
fi


# Set and export location variables, for easy access in R.

# Constant data
export CONST_GOLD_CALLS_VCFGZ="${PATH_RESOURCES_HEAD}/gold_standard/calls-2.19.vcf.gz"
export CONST_GOLD_CALLS_VCFGZTBI="${PATH_RESOURCES_HEAD}/gold_standard/calls-2.19.vcf.gz.tbi"
export CONST_GOLD_HARDMASK_VALID_REGIONS_BEDGZ="${PATH_RESOURCES_HEAD}/gold_standard/valid_regions-2.19.bed.gz"
export CONST_REFERENCE_SDF="${PATH_RESOURCES_HEAD}/reference/ref.sdf/"
export CONST_FUNCTIONAL_REGIONS_BEDGZ_PREFIX="${PATH_RESOURCES_HEAD}/functional_regions/"
export CONST_MASK_REGIONS_BEDGZ_PREFIX="${PATH_RESOURCES_HEAD}/mask_regions/"
export CONST_REFERENCE_BSGENOME="BSgenome.HSapiens.1000g.37d5"		# This is a custom package, available at /share/ClusterShare/biodata/contrib/marpin/reference/hs37d5/build/BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz

# Script location
export PARAM_EXEC_PATH=$(pwd)

# Scratch space
mkdir -p ${PATH_SCRATCH_DEFAULT}
export PARAM_SCRATCH=$(mktemp -d --tmpdir=${PATH_SCRATCH_DEFAULT} valrept.XXXXXXXXXX)
export PARAM_INPUT_SCRATCH="${PARAM_SCRATCH}/input"
export PARAM_RTG_OVERLAP_SCRATCH="${PARAM_SCRATCH}/overlap"
export PARAM_KNITR_SCRATCH="${PARAM_SCRATCH}/knitr"

# Program parameters
export PARAM_INPUT_VCFGZ_PATH
export PARAM_REGION_BED_SUPPLIED
export PARAM_REGION_BED_PATH
export PARAM_OUTPUT_PDF_PATH
export PARAM_EXTENDED
export PARAM_VERSION_EXEC_HOST
export PARAM_VERSION_RTG
export PARAM_VERSION_JAVA
export PARAM_VERSION_BEDTOOLS

# Temporary file locations
export PATH_TEST_VARIANTS="${PARAM_INPUT_SCRATCH}/test_variants.vcf.gz"
export PATH_TEST_VARIANTS_INDEX="${PATH_TEST_VARIANTS}.tbi"
export PATH_GOLD_VARIANTS="${PARAM_INPUT_SCRATCH}/gold_variants.vcf.gz"
export PATH_GOLD_VARIANTS_INDEX="${PATH_GOLD_VARIANTS}.tbi"
export PATH_GOLD_REGIONS="${PARAM_INPUT_SCRATCH}/gold_regions.bed.gz"

# Overlap file locations
export PATH_OVERLAP_TP="${PARAM_RTG_OVERLAP_SCRATCH}/tp.vcf.gz"
export PATH_OVERLAP_FP="${PARAM_RTG_OVERLAP_SCRATCH}/fp.vcf.gz"
export PATH_OVERLAP_FN="${PARAM_RTG_OVERLAP_SCRATCH}/fn.vcf.gz"


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
PARAM_INPUT_VCFGZ_PATH=""
PARAM_REGION_BED_SUPPLIED=0
PARAM_REGION_BED_PATH="NA"
PARAM_OUTPUT_PDF_PATH="${PARAM_EXEC_PATH}/report.pdf"
PARAM_EXTENDED=0

while getopts "r:o:hx" opt; do
	case "$opt" in
		h)
			print_usage
			exit 0
			;;
		o)
			PARAM_OUTPUT_PDF_PATH=$(readlink -f $OPTARG)
			;;
		r)
			PARAM_REGION_BED_SUPPLIED=1
			PARAM_REGION_BED_PATH=$(readlink -f $OPTARG)
			;;
		x)
			PARAM_EXTENDED=1
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

PARAM_INPUT_VCFGZ_PATH=$1


#####################################################################
# PARAMETER CHECKING
#####################################################################

if [ ! -d "${PATH_RESOURCES_HEAD}" ]; then
  echo >&2 "Error: Resources path ${PATH_RESOURCES_HEAD} not found."
  exit 2
fi

if [ ! -e ${PARAM_INPUT_VCFGZ_PATH} ]; then
	echo >&2 "Error: Input file ${PARAM_INPUT_VCFGZ_PATH} not found."
	exit 3
fi

if [ -e ${PARAM_OUTPUT_PDF_PATH} ]; then
	echo >&2 "Error: Output file ${PARAM_OUTPUT_PDF_PATH} already exists."
	exit 4
fi

if [ ${PARAM_REGION_BED_SUPPLIED} -eq 1 ] && [ ! -e ${PARAM_REGION_BED_PATH} ]; then
  echo >&2 "Error: Region file ${PARAM_REGION_BED_PATH} not found."
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

if [ ! -e ${RTG_CORE} ]; then
  echo >&2 "Error: RTG-core.jar not found."
  exit 8
fi

#set -e should catch this.
[[ -f report.Rnw ]] && [[ -f report_functions.R ]] && [[ -f report_debug.Rnw ]] && [[ -f report_calculations.R ]]


#####################################################################
# R PACKAGE CHECKING
#####################################################################
${R} --vanilla -e "if (!(\"${CONST_REFERENCE_BSGENOME}\" %in% installed.packages(.Library))) stop(\"${CONST_REFERENCE_BSGENOME} package not installed\")"


#####################################################################
# VERSIONING
#####################################################################
PARAM_VERSION_EXEC_HOST=$(uname -a)
PARAM_VERSION_RTG=$(${JAVA} -jar ${RTG_CORE} version | grep 'Core Version: ' | sed 's/.*: //g')
PARAM_VERSION_JAVA=$(${JAVA} -version 2>&1 | head -n 1 | sed -E 's/[^"]+"//;s/"$//')
PARAM_VERSION_BEDTOOLS=$(${BEDTOOLS} --version | cut -d' ' -f 2)


#####################################################################
# CREATE TEMPORARY DIRECTORIES
#####################################################################
mkdir -p ${PARAM_SCRATCH}
mkdir -p ${PARAM_KNITR_SCRATCH}
mkdir -p ${PARAM_INPUT_SCRATCH}


#####################################################################
# SUBSET TO REGION BED
#####################################################################

if [ ${PARAM_REGION_BED_SUPPLIED} -eq 1 ]; then
  echo "Subsetting input files to supplied BED..."
  # Sort the region bed
  sort -k1,1 -k2,2n ${PARAM_REGION_BED_PATH} > ${PARAM_INPUT_SCRATCH}/region.bed

  # Perform the intersection
  ${BEDTOOLS} intersect -header -wa -a ${PARAM_INPUT_VCFGZ_PATH} -b ${PARAM_INPUT_SCRATCH}/region.bed | ${BGZIP} > ${PATH_TEST_VARIANTS}
  ${BEDTOOLS} intersect -header -wa -a ${CONST_GOLD_CALLS_VCFGZ} -b ${PARAM_INPUT_SCRATCH}/region.bed | ${BGZIP} > ${PATH_GOLD_VARIANTS}
  ${BEDTOOLS} intersect -wa -a ${CONST_GOLD_HARDMASK_VALID_REGIONS_BEDGZ} -b ${PARAM_INPUT_SCRATCH}/region.bed | ${BGZIP} > ${PATH_GOLD_REGIONS}

  # We need to re-index the gold variants
  ${TABIX} -p vcf ${PATH_GOLD_VARIANTS}
else
  # No region bed supplied; copy over the files in their entirety
  cp ${PARAM_INPUT_VCFGZ_PATH} ${PATH_TEST_VARIANTS}
  cp ${CONST_GOLD_CALLS_VCFGZ} ${PATH_GOLD_VARIANTS}
  cp ${CONST_GOLD_HARDMASK_VALID_REGIONS_BEDGZ} ${PATH_GOLD_REGIONS}

  # We can use the gold standard index unchanged, so no need to 
  # index it as above.
  cp ${CONST_GOLD_CALLS_VCFGZTBI} ${PATH_GOLD_VARIANTS_INDEX}
fi

# Regardless of whether a region bed was supplied or not, we still
# need to index the test variant vcf.
${TABIX} -p vcf ${PATH_TEST_VARIANTS}



#####################################################################
# VCF OVERLAP EVALUATION
#####################################################################
echo "Computing VCF overlaps..."

if [ -e ${PARAM_RTG_OVERLAP_SCRATCH} ]; then
	echo "Overlap scratch directory ${PARAM_RTG_OVERLAP_SCRATCH} already exists.  Clearing scratch directory and continuing..."
	rm -rf ${PARAM_RTG_OVERLAP_SCRATCH}
fi

${RTG_VCFEVAL} --all-records -b ${PATH_GOLD_VARIANTS} -c ${PATH_TEST_VARIANTS} -t ${CONST_REFERENCE_SDF} -o ${PARAM_RTG_OVERLAP_SCRATCH}

#####################################################################
# DIRTY DIRTY DIRTY
# For some reason, dx is running RTG (via java) as root, and so 
# output files are owned by root.  Unfortunately, later R invocations
# are run as dnanexus, and can't access the root-owned output files
# from RTG vcfeval.  Get around this by explicitly chowning files 
# dnanexus.  We run as root at this point, so no sudo is needed.
# (My current guess is that R drops root, but it's just a guess
# right now.)
if [ ${IS_DNANEXUS} -eq 1 ]; then
  chown -R dnanexus:dnanexus ${PARAM_SCRATCH}
fi
# DIRTY DIRTY DIRTY
#####################################################################

# TODO: Parse ${PARAM_RTG_OVERLAP_SCRATCH}/vcfeval.log to identify regions to exclude
# eg Evaluation too complex (5001 unresolved paths, 18033 iterations) at reference region 2:105849275-105849281. Variants in this region will not be included in results.
# This will be required to remove the TNs in this region, but it's polish.


#####################################################################
# CALCULATIONS FOR REPORT
#####################################################################
echo "Performing calculations for report..."

# knitr doesn't play well with building knits outside of its working
# directory.  Currently we get around this with a bit of a kludge, 
# by copying the report files to ${PARAM_KNITR_SCRATCH}, then executing in 
# that directory.
cp -f report.Rnw ${PARAM_KNITR_SCRATCH}
cp -f report_functions.R ${PARAM_KNITR_SCRATCH}
cp -f report_extended.Rnw ${PARAM_KNITR_SCRATCH}
cp -f report_calculations.R ${PARAM_KNITR_SCRATCH}
cd ${PARAM_KNITR_SCRATCH}

# Run the script.  All options are passed via exported environment 
# variables.  Also save these variables to a file for later source-ing,
# to ease debugging.
export > environment
${RSCRIPT} --vanilla report_calculations.R


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
rm -f ${PARAM_KNITR_SCRATCH}/report.pdf

# Run pdflatex
pdflatex -interaction nonstopmode report.tex
pdflatex -interaction nonstopmode report.tex

# Check  whether the report.pdf was generated
if [ ! -e ${PARAM_KNITR_SCRATCH}/report.pdf ]; then
	echo >&2 "Error: pdflatex did not successfully generate report.pdf."
	echo >&2 "Check ${PARAM_KNITR_SCRATCH}/report.tex and the latex log for errors."
	exit 9
fi

# Copy the completed report to the final destination
cp "${PARAM_KNITR_SCRATCH}/report.pdf" "${PARAM_OUTPUT_PDF_PATH}"

echo "Report generated successfully."

echo "Done."
