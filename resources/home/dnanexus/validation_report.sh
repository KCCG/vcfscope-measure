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

IS_DNANEXUS=0
if [ `whoami` == dnanexus ]; then
  let IS_DNANEXUS=1
fi

#
# Software & Resources
#
if [ ${IS_DNANEXUS} ]; then
  RSCRIPT="/home/dnanexus/bin/Rscript"
  JAVA=`which java`
  GIT=`which git`
  RESOURCES_HEAD="/home/dnanexus/resources"
  SCRATCH_DEFAULT="/tmp"

  RTG_CORE="${RESOURCES_HEAD}/rtg-core/rtg-core.jar"
  RTG_THREADS=`nproc`

  # Calculate 80% of memory size, for java
  mem_in_mb=`head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}'`
  RTG_VCFEVAL="${JAVA} -Xmx${mem_in_mb}m -jar ${RTG_CORE} vcfeval -T ${RTG_THREADS}"
else
  RSCRIPT="/home/marpin/bin/Rscript"
  JAVA="/usr/java/latest/bin/java"
  GIT="/home/marpin/bin/git"
  RESOURCES_HEAD="/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources"
  SCRATCH_DEFAULT="/directflow/ClinicalGenomicsPipeline/tmp"

  RTG_CORE="${RESOURCES_HEAD}/rtg-core/rtg-core.jar"
  RTG_THREADS=4
  RTG_VCFEVAL="${JAVA} -Xmx8G -jar ${RTG_CORE} vcfeval -T ${RTG_THREADS}"
fi


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
region_file_supplied=0
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
			region_file_supplied=1
			region_file=$(readlink -f $OPTARG)
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
  echo -e >&2 "Error: Resources path ${RESOURCES_HEAD} not found."
  exit 10
fi

if [ ! -e ${input_vcfgz_path} ]; then
	echo -e "Error: Input file ${input_vcfgz_path} not found."
	exit 2
fi

if [ -e ${output_pdf_path} ]; then
	echo -e "Output file ${output_pdf_path} already exists; quitting."
	exit 3
fi


#####################################################################
# SOFTWARE CHECKING
#####################################################################
if [ ! -e ${RSCRIPT} ]; then
  echo -e "Error: Rscript executable not found."
  exit 2
fi
if [ ! -e ${JAVA} ]; then
  echo -e "Error: java executable not found."
  exit 2
fi
if [ ! -e ${GIT} ]; then
  echo -e "Error: git executable not found."
  exit 2
fi

#set -e should catch this.
[[ -f report.Rnw ]] && [[ -f report_functions.R ]] && [[ -f report_debug.Rnw ]] && [[ -f report_calculations.R ]]


#####################################################################
# R PACKAGE CHECKING
#####################################################################
R --vanilla -e "suppressWarnings(require(\"${REFERENCE_BSGENOME}\", quiet=TRUE)) || stop(\"${REFERENCE_BSGENOME} package not installed\")"


#####################################################################
# SOFTWARE VERSIONING
#####################################################################
if [ ${IS_DNANEXUS} ]; then
  VERSION_GIT_BRANCH="unknown"
  VERSION_GIT_COMMIT="unknown"
else
  VERSION_GIT_BRANCH=$(${GIT} rev-parse --abbrev-ref HEAD)
  VERSION_GIT_COMMIT=$(${GIT} rev-parse --verify HEAD)
fi
VERSION_EXEC_HOST=$(uname -a)


#####################################################################
# CREATE TEMPORARY DIRECTORIES
#####################################################################
mkdir -p ${SCRATCH}
mkdir -p ${KNITR_SCRATCH}


#####################################################################
# VCF OVERLAP EVALUATION
#####################################################################
echo -e "Computing VCF overlaps..."

if [ -e ${RTG_OVERLAP_SCRATCH} ]; then
	echo -e "Overlap scratch directory ${RTG_OVERLAP_SCRATCH} already exists.  Clearing scratch directory and continuing..."
	rm -rf ${RTG_OVERLAP_SCRATCH}
fi

${RTG_VCFEVAL} --all-records -b ${GOLD_CALLS_VCFGZ} -c ${input_vcfgz_path} -t ${REFERENCE_SDF} -o ${RTG_OVERLAP_SCRATCH} > /dev/null 2>&1

# TODO: Parse ${RTG_OVERLAP_SCRATCH}/vcfeval.log to identify regions to exclude
# eg Evaluation too complex (5001 unresolved paths, 18033 iterations) at reference region 2:105849275-105849281. Variants in this region will not be included in results.
# This will be required to remove the TNs in this region, but it's polish.


#####################################################################
# CALCULATIONS FOR REPORT
#####################################################################
echo -e "Performing calculations for report..."
#mkdir -p ${KNITR_SCRATCH}
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
${RSCRIPT} --vanilla report_calculations.R ${debug} ${debug_chrom} ${extended} ${input_vcfgz_path} \
  ${RTG_OVERLAP_SCRATCH}/tp.vcf.gz ${RTG_OVERLAP_SCRATCH}/fp.vcf.gz ${RTG_OVERLAP_SCRATCH}/fn.vcf.gz \
  ${GOLD_CALLS_VCFGZ} ${REFERENCE_BSGENOME} ${GOLD_HARDMASK_VALID_REGIONS_BEDGZ} \
  ${FUNCTIONAL_REGIONS_BEDGZ_PREFIX} ${MASK_REGIONS_BEDGZ_PREFIX} \
   "'""${VERSION}""'" "${VERSION_GIT_BRANCH}" "${VERSION_GIT_COMMIT}" "${VERSION_EXEC_HOST}"


#####################################################################
# REPORT GENERATION
#####################################################################
echo -e "Generating report..."

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
	echo -e "Error: pdflatex did not successfully generate a report.pdf."
	echo -e "Check ${KNITR_SCRATCH}/report.tex and the latex log for errors."
	exit 4
fi

# Copy the completed report to the final destination
cp "${KNITR_SCRATCH}/report.pdf" "${output_pdf_path}"

echo -e "Report generated successfully."

echo -e "Done."
