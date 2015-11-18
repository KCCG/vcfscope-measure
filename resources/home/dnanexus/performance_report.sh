#!/bin/bash
set -e -u -o pipefail
set +x
IFS=$'\n\t'

#####################################################################
# VERSION
#####################################################################
export CONST_VERSION_SCRIPT="2.0.0"


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
  PYTHON=`which python`
  BEDTOOLS=`which bedtools`
  TABIX=`which tabix`
  BGZIP=`which bgzip`
  BCFTOOLS=`which bcftools`
  GHOSTSCRIPT=`which gs`
  PARALLEL="`which parallel` --gnu --halt now,fail=1"
  SAMTOOLS=`which samtools`

  RTG_TOOLS="${PATH_RESOURCES_HEAD}/rtg-tools/RTG.jar"
  RTG_THREADS=`nproc`
  mem_in_mb=`head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}'`                 # Calculate 80% of memory size, for java
  RTG_VCFEVAL="${JAVA} -Xmx${mem_in_mb}m -jar ${RTG_TOOLS} vcfeval -T ${RTG_THREADS}"
else
  # Wolfpack settings (marpin only for now)
  PATH_RESOURCES_HEAD="/directflow/ClinicalGenomicsPipeline/projects/performance-reporter/kccg_performance_reporter_resources_bundle-2.0"
  PATH_SCRATCH_DEFAULT="/directflow/ClinicalGenomicsPipeline/tmp"

  RSCRIPT="/home/marpin/bin/Rscript"
  R="/home/marpin/bin/R"
  JAVA="/usr/java/latest/bin/java"
  PYTHON=`which python`
  BEDTOOLS="/home/marpin/software/bedtools2/bin/bedtools"
  TABIX="/home/marpin/software/htslib/tabix"
  BGZIP="/home/marpin/software/htslib/bgzip"
  BCFTOOLS="/home/marpin/software/bcftools/bcftools"
  GHOSTSCRIPT=`which gs`
  PARALLEL="/home/marpin/software/parallel/parallel --gnu --halt now,fail=1"
  SAMTOOLS=`which samtools`

  RTG_TOOLS="${PATH_RESOURCES_HEAD}/rtg-tools/RTG.jar"
  RTG_THREADS=4
  mem_in_mb=8192
  RTG_VCFEVAL="${JAVA} -Xmx${mem_in_mb}m -jar ${RTG_TOOLS} vcfeval -T ${RTG_THREADS}"
fi


# Set and export location variables, for easy access in R.

# Constant data
export CONST_GOLD_CALLS_VCFGZ="${PATH_RESOURCES_HEAD}/gold_standard/calls-2.19.vcf.gz"
export CONST_GOLD_CALLS_VCFGZTBI="${PATH_RESOURCES_HEAD}/gold_standard/calls-2.19.vcf.gz.tbi"
export CONST_GOLD_HARDMASK_VALID_REGIONS_BEDGZ="${PATH_RESOURCES_HEAD}/gold_standard/valid_regions-2.19.bed.gz"
export CONST_REFERENCE_SDF="${PATH_RESOURCES_HEAD}/reference/hs37d5.sdf/"
export CONST_GENOME_BEDGZ="${PATH_RESOURCES_HEAD}/reportable_range/genome.bed.gz"
export CONST_REFERENCE_BSGENOME="BSgenome.HSapiens.1000g.37d5"		# This is a custom package, available at /share/ClusterShare/biodata/contrib/marpin/reference/hs37d5/build/BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz

# Script location
export PARAM_SCRIPT_PATH=$(readlink -f $(dirname $0))
export PARAM_EXEC_PATH=$(pwd)

# Scratch space
mkdir -p ${PATH_SCRATCH_DEFAULT}
export PARAM_SCRATCH=$(mktemp -d --tmpdir=${PATH_SCRATCH_DEFAULT} perfrept.XXXXXXXXXX)
export PARAM_INPUT_SCRATCH="${PARAM_SCRATCH}/input"
export PARAM_RTG_OVERLAP_SCRATCH="${PARAM_SCRATCH}/overlap"
export PARAM_KNITR_SCRATCH="${PARAM_SCRATCH}/knitr"
export PARAM_DEPTH_SCRATCH="${PARAM_SCRATCH}/depth"

# Program parameters
export PARAM_INPUT_VCFGZ_PATH
export PARAM_INPUT_VCF_SAMPLES
export PARAM_REGION_BED_SUPPLIED
export PARAM_REGION_BED_PATH
export PARAM_OUTPUT_PDF_PATH
export PARAM_OUTPUT_RDS_PATH
export PARAM_EXTENDED
export PARAM_VERSION_EXEC_HOST
export PARAM_VERSION_RTG
export PARAM_VERSION_JAVA
export PARAM_VERSION_BEDTOOLS

# Temporary file locations
export PATH_TEST_VARIANTS="${PARAM_INPUT_SCRATCH}/test_variants.vcf.gz"
export PATH_TEST_VARIANTS_INDEX="${PATH_TEST_VARIANTS}.tbi"
export PATH_TEST_READS="${PARAM_INPUT_SCRATCH}/test_reads.bam"
export PATH_TEST_READS_INDEX="${PATH_TEST_READS}.bai"
export PATH_GOLD_VARIANTS="${PARAM_INPUT_SCRATCH}/gold_variants.vcf.gz"
export PATH_GOLD_VARIANTS_INDEX="${PATH_GOLD_VARIANTS}.tbi"
export PATH_GOLD_REGIONS="${PARAM_INPUT_SCRATCH}/gold_regions.bed.gz"

# Overlap file locations.
export PATH_SAMPLE_OVERLAP_TP
export PATH_SAMPLE_OVERLAP_FP
export PATH_SAMPLE_OVERLAP_FN


#####################################################################
# COMMAND LINE PARSING
#####################################################################
print_usage() {
cat << EOF
Usage: ${0##*/} [-o OUTFILE] [-d RDSOUT] [-r BEDFILE] [-t] <VCF> <BAM>

Create a single-sample WGS performance report.

    VCF          Input NA12878 genotype calls, in vcf.gz format.
    BAM          Mapped NA12878 reads used to generate the VCF calls.
                 An accompanying index (<BAM>.bai) must be present in 
                 the same directory.
    -o OUTFILE   Write the report to OUTFILE (default: report.pdf)
    -d RDSOUT    Write performance report data to RDSOUT (default: 
                 not written)
    -r BEDFILE   Restrict analysis to the regions in BEDFILE only.
                 Default: the full genome is considered.
    -t           Perfom regression tests and consistency checks 
                 prior to report generation.  Generally only for
                 development use.  Default: tests are not performed.
    -h           Display this help and exit.

Version ${CONST_VERSION_SCRIPT}

Mark Pinese
EOF
}

# Head location for resources bundle
OPTIND=1
PARAM_INPUT_VCFGZ_PATH=""
PARAM_REGION_BED_SUPPLIED=0
PARAM_REGION_BED_PATH="NA"
PARAM_OUTPUT_PDF_PATH="${PARAM_SCRIPT_PATH}/performance_report.pdf"
PARAM_OUTPUT_RDS_PATH=""
PARAM_DOTESTS=0

while getopts "r:o:d:ht" opt; do
	case "$opt" in
		h)
			print_usage
			exit 0
			;;
		o)
			PARAM_OUTPUT_PDF_PATH=$(readlink -f "${OPTARG}")
			;;
    d)
      PARAM_OUTPUT_RDS_PATH=$(readlink -f "${OPTARG}")
      ;;
		r)
			PARAM_REGION_BED_SUPPLIED=1
			PARAM_REGION_BED_PATH=$(readlink -f "${OPTARG}")
			;;
    t)
      PARAM_DOTESTS=1
      ;;
		'?')
			print_usage >&2
			exit 1
			;;
	esac
done

shift $((OPTIND-1))

if [ $# -ne 2 ]; then
	print_usage >&2
	exit 1
fi

PARAM_INPUT_VCFGZ_PATH=$(readlink -f $1)
PARAM_INPUT_BAM_PATH=$(readlink -f $2)


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

if [ ! -e ${PARAM_INPUT_BAM_PATH} ]; then
  echo >&2 "Error: Input file ${PARAM_INPUT_BAM_PATH} not found."
  exit 3
fi

if [ ! -e ${PARAM_INPUT_BAM_PATH}.bai ]; then
  echo >&2 "Error: Index not found for input BAM file (expected ${PARAM_INPUT_BAM_PATH}.bai)"
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
  echo >&2 "Error: Rscript executable not found at ${RSCRIPT}."
  exit 6
fi

if [ ! -e ${JAVA} ]; then
  echo >&2 "Error: Java executable not found at ${JAVA}."
  exit 7
fi

if [ ! -e ${RTG_TOOLS} ]; then
  echo >&2 "Error: RTG.jar not found at ${RTG_TOOLS}."
  exit 8
fi

if [ ! -f ${PARAM_SCRIPT_PATH}/report.Rnw -o ! -f ${PARAM_SCRIPT_PATH}/report_functions.R -o ! -f ${PARAM_SCRIPT_PATH}/report_calculations.R ]; then
  echo >&2 "Error: Missing at least one required R source file."
  exit 9
fi


#####################################################################
# CHECK THE INPUT FILE IS A SINGLE-SAMPLE VCF
#####################################################################
if [ $(gzip -dc ${PARAM_INPUT_VCFGZ_PATH} | grep -m1 '^#[^#]' | tr '\t' '\n' | wc -l) -ne 10 ]; then
  echo >&2 "Error: Multi-sample VCF found.  Input VCF must contain variants for only a single sample."
  exit 10
fi


#####################################################################
# R PACKAGE CHECKING
#####################################################################
${R} --vanilla --slave -e "if (!(\"${CONST_REFERENCE_BSGENOME}\" %in% installed.packages(.Library))) stop(\"${CONST_REFERENCE_BSGENOME} package not installed\")"


#####################################################################
# VERSIONING
#####################################################################
PARAM_VERSION_EXEC_HOST=$(uname -a)
PARAM_VERSION_RTG=$(${JAVA} -Xmx${mem_in_mb}m -jar ${RTG_TOOLS} version | grep -m1 'Product: ' | sed 's/.* //g')
PARAM_VERSION_RTG="${PARAM_VERSION_RTG} $(${JAVA} -Xmx${mem_in_mb}m -jar ${RTG_TOOLS} version | grep -m1 'Core Version: ' | sed 's/.*: //g')"
PARAM_VERSION_JAVA=$(${JAVA} -Xmx${mem_in_mb}m -version 2>&1 | head -n 1 | sed -E 's/[^"]+"//;s/"$//')
PARAM_VERSION_BEDTOOLS=$(${BEDTOOLS} --version | cut -d' ' -f 2)


#####################################################################
# PARAMETER LOGGING
#####################################################################
echo >&2 "performance_report.sh parameters:"
echo >&2 "  PATH_RESOURCES_HEAD=${PATH_RESOURCES_HEAD}"
echo >&2 "  PATH_SCRATCH_DEFAULT=${PATH_SCRATCH_DEFAULT}"
echo >&2 "  BEDTOOLS=${BEDTOOLS}"
echo >&2 "  BCFTOOLS=${BCFTOOLS}"
echo >&2 "  BGZIP=${BGZIP}"
echo >&2 "  GHOSTSCRIPT=${GHOSTSCRIPT}"
echo >&2 "  JAVA=${JAVA}"
echo >&2 "  PARALLEL=${PARALLEL}"
echo >&2 "  PYTHON=${PYTHON}"
echo >&2 "  R=${R}"
echo >&2 "  RSCRIPT=${RSCRIPT}"
echo >&2 "  RTG_TOOLS=${RTG_TOOLS}"
echo >&2 "  RTG_VCFEVAL=${RTG_VCFEVAL}"
echo >&2 "  SAMTOOLS=${SAMTOOLS}"
echo >&2 "  TABIX=${TABIX}"
echo >&2 "  RTG_THREADS=${RTG_THREADS}"
echo >&2 "  mem_in_mb=${mem_in_mb}"
echo >&2 "  CONST_GOLD_CALLS_VCFGZ=${CONST_GOLD_CALLS_VCFGZ}"
echo >&2 "  CONST_GOLD_CALLS_VCFGZTBI=${CONST_GOLD_CALLS_VCFGZTBI}"
echo >&2 "  CONST_GOLD_HARDMASK_VALID_REGIONS_BEDGZ=${CONST_GOLD_HARDMASK_VALID_REGIONS_BEDGZ}"
echo >&2 "  CONST_REFERENCE_SDF=${CONST_REFERENCE_SDF}"
echo >&2 "  CONST_GENOME_BEDGZ=${CONST_GENOME_BEDGZ}"
echo >&2 "  CONST_REFERENCE_BSGENOME=${CONST_REFERENCE_BSGENOME}"
echo >&2 "  PARAM_SCRIPT_PATH=${PARAM_SCRIPT_PATH}"
echo >&2 "  PARAM_SCRATCH=${PARAM_SCRATCH}"
echo >&2 "  PARAM_INPUT_SCRATCH=${PARAM_INPUT_SCRATCH}"
echo >&2 "  PARAM_RTG_OVERLAP_SCRATCH=${PARAM_RTG_OVERLAP_SCRATCH}"
echo >&2 "  PARAM_KNITR_SCRATCH=${PARAM_KNITR_SCRATCH}"
echo >&2 "  PARAM_INPUT_VCFGZ_PATH=${PARAM_INPUT_VCFGZ_PATH}"
echo >&2 "  PARAM_INPUT_BAM_PATH=${PARAM_INPUT_BAM_PATH}"
echo >&2 "  PARAM_REGION_BED_SUPPLIED=${PARAM_REGION_BED_SUPPLIED}"
echo >&2 "  PARAM_REGION_BED_PATH=${PARAM_REGION_BED_PATH}"
echo >&2 "  PARAM_OUTPUT_PDF_PATH=${PARAM_OUTPUT_PDF_PATH}"
echo >&2 "  PARAM_OUTPUT_RDS_PATH=${PARAM_OUTPUT_RDS_PATH}"
echo >&2 "  PARAM_VERSION_EXEC_HOST=${PARAM_VERSION_EXEC_HOST}"
echo >&2 "  PARAM_VERSION_RTG=${PARAM_VERSION_RTG}"
echo >&2 "  PARAM_VERSION_JAVA=${PARAM_VERSION_JAVA}"
echo >&2 "  PARAM_VERSION_BEDTOOLS=${PARAM_VERSION_BEDTOOLS}"
echo >&2 "  PATH_TEST_VARIANTS=${PATH_TEST_VARIANTS}"
echo >&2 "  PATH_TEST_VARIANTS_INDEX=${PATH_TEST_VARIANTS_INDEX}"
echo >&2 "  PATH_TEST_READS=${PATH_TEST_READS}"
echo >&2 "  PATH_TEST_READS_INDEX=${PATH_TEST_READS_INDEX}"
echo >&2 "  PATH_GOLD_VARIANTS=${PATH_GOLD_VARIANTS}"
echo >&2 "  PATH_GOLD_VARIANTS_INDEX=${PATH_GOLD_VARIANTS_INDEX}"
echo >&2 "  PATH_GOLD_REGIONS=${PATH_GOLD_REGIONS}"


#####################################################################
# CREATE TEMPORARY DIRECTORIES
#####################################################################
mkdir -p ${PARAM_SCRATCH}
mkdir -p ${PARAM_INPUT_SCRATCH}
mkdir -p ${PARAM_DEPTH_SCRATCH}


#####################################################################
# CREATE OUTPUT DIRECTORIES
#####################################################################
mkdir -p $(dirname ${PARAM_OUTPUT_PDF_PATH})
[ -e ${PARAM_OUTPUT_RDS_PATH} ] || mkdir -p $(dirname ${PARAM_OUTPUT_RDS_PATH})


#####################################################################
# SUBSET TO REGION BED
#####################################################################
if [ ${PARAM_REGION_BED_SUPPLIED} -eq 1 ]; then
  echo "Subsetting input files to supplied BED..."
  # Sort the region bed
  sort -k1,1 -k2,2n ${PARAM_REGION_BED_PATH} > ${PARAM_INPUT_SCRATCH}/region.bed

  # Perform the VCF intersection
  ${BEDTOOLS} intersect -header -wa -a ${PARAM_INPUT_VCFGZ_PATH} -b ${PARAM_INPUT_SCRATCH}/region.bed | ${BGZIP} > ${PATH_TEST_VARIANTS}
  ${BEDTOOLS} intersect -header -wa -a ${CONST_GOLD_CALLS_VCFGZ} -b ${PARAM_INPUT_SCRATCH}/region.bed | ${BGZIP} > ${PATH_GOLD_VARIANTS}
  ${BEDTOOLS} intersect -wa -a ${CONST_GOLD_HARDMASK_VALID_REGIONS_BEDGZ} -b ${PARAM_INPUT_SCRATCH}/region.bed | ${BGZIP} > ${PATH_GOLD_REGIONS}

  # We need to re-index the gold variants
  ${TABIX} -p vcf ${PATH_GOLD_VARIANTS}
else
  # No region bed supplied; link the unmodified variant files
  ln ${PARAM_INPUT_VCFGZ_PATH} ${PATH_TEST_VARIANTS}
  ln ${CONST_GOLD_CALLS_VCFGZ} ${PATH_GOLD_VARIANTS}
  ln ${CONST_GOLD_HARDMASK_VALID_REGIONS_BEDGZ} ${PATH_GOLD_REGIONS}

  # We can use the gold standard index unchanged, so no need to 
  # index it as above.
  ln ${CONST_GOLD_CALLS_VCFGZTBI} ${PATH_GOLD_VARIANTS_INDEX}
fi

# Link the reads and index unmodified.  I tried to subset here, 
# but bedtools intersect on the bam was way too slow.
ln ${PARAM_INPUT_BAM_PATH} ${PATH_TEST_READS}
ln ${PARAM_INPUT_BAM_PATH}.bai ${PATH_TEST_READS_INDEX}

# Regardless of whether a region bed was supplied or not, we still
# need to index the test variant vcf.
${TABIX} -p vcf ${PATH_TEST_VARIANTS}


#####################################################################
# VCF OVERLAP EVALUATION
#####################################################################
echo "Computing VCF overlaps..."
if [ -e ${PARAM_RTG_OVERLAP_SCRATCH} ]; then
  echo "PARAM_RTG_OVERLAP_SCRATCH directory ${PARAM_RTG_OVERLAP_SCRATCH} already exists.  Clearing scratch directory and continuing..."
  rm -rf ${PARAM_RTG_OVERLAP_SCRATCH}
fi

eval ${RTG_VCFEVAL} --all-records -b ${PATH_GOLD_VARIANTS} -c ${PATH_TEST_VARIANTS} -t ${CONST_REFERENCE_SDF} -o ${PARAM_RTG_OVERLAP_SCRATCH}

PATH_SAMPLE_OVERLAP_TP_RAW="${PARAM_RTG_OVERLAP_SCRATCH}/tp.vcf.gz"
PATH_SAMPLE_OVERLAP_FP_RAW="${PARAM_RTG_OVERLAP_SCRATCH}/fp.vcf.gz"
PATH_SAMPLE_OVERLAP_FN_RAW="${PARAM_RTG_OVERLAP_SCRATCH}/fn.vcf.gz"

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
  chown -R dnanexus:dnanexus $(dirname ${PARAM_OUTPUT_PDF_PATH})
  [ -e ${PARAM_OUTPUT_RDS_PATH} ] || chown -R dnanexus:dnanexus $(dirname ${PARAM_OUTPUT_RDS_PATH})
fi
# DIRTY DIRTY DIRTY
#####################################################################

# TODO: Parse ${PARAM_RTG_OVERLAP_SCRATCH}/vcfeval.log to identify regions to exclude
# eg Evaluation too complex (5001 unresolved paths, 18033 iterations) at reference region 2:105849275-105849281. Variants in this region will not be included in results.
# This will be required to remove the TNs in this region, but it's polish.


#####################################################################
# CALCULATE BAM DEPTH AT VARIANT LOCATIONS
#####################################################################
echo "Calculating depth at variant loci..."
# Make beds of regions for which we want to know achieved depth
gzip -dc ${PATH_SAMPLE_OVERLAP_TP_RAW} | grep -v '^#' | awk 'BEGIN { FS="\t"; OFS="\t" } { print $1, $2-1, $2+length($4)-1 }' > ${PARAM_DEPTH_SCRATCH}/tp.bed
gzip -dc ${PATH_SAMPLE_OVERLAP_FP_RAW} | grep -v '^#' | awk 'BEGIN { FS="\t"; OFS="\t" } { print $1, $2-1, $2+length($4)-1 }' > ${PARAM_DEPTH_SCRATCH}/fp.bed
gzip -dc ${PATH_SAMPLE_OVERLAP_FN_RAW} | grep -v '^#' | awk 'BEGIN { FS="\t"; OFS="\t" } { print $1, $2-1, $2+length($4)-1 }' > ${PARAM_DEPTH_SCRATCH}/fn.bed

# Shard beds and the bam by chromosome
# Shard the beds:
awk -v dir=${PARAM_DEPTH_SCRATCH} '{ print > dir"/tp_"$1".bed" }' < ${PARAM_DEPTH_SCRATCH}/tp.bed
awk -v dir=${PARAM_DEPTH_SCRATCH} '{ print > dir"/fp_"$1".bed" }' < ${PARAM_DEPTH_SCRATCH}/fp.bed
awk -v dir=${PARAM_DEPTH_SCRATCH} '{ print > dir"/fn_"$1".bed" }' < ${PARAM_DEPTH_SCRATCH}/fn.bed

# And the bam:
echo "  Sharding BAM..."
${SAMTOOLS} view -H ${PATH_TEST_READS} | grep '^@SQ' | cut -f 2 | sed 's/^SN://' > ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt ${SAMTOOLS} view -1 -o ${PARAM_DEPTH_SCRATCH}/reads_{1}.bam ${PATH_TEST_READS} {1}
eval ${PARALLEL} ${SAMTOOLS} index ::: ${PARAM_DEPTH_SCRATCH}/reads_*.bam

# Fill in empty beds for chromosomes in the BAM but not the beds
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt touch ${PARAM_DEPTH_SCRATCH}/tp_{1}.bed
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt touch ${PARAM_DEPTH_SCRATCH}/fp_{1}.bed
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt touch ${PARAM_DEPTH_SCRATCH}/fn_{1}.bed

# For each shard (chromosome), calculate depth in the BAM in the VCF regions
echo "  True positives..."
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt ${PYTHON} ${PARAM_SCRIPT_PATH}/depth.py ${PARAM_DEPTH_SCRATCH}/reads_{1}.bam ${PARAM_DEPTH_SCRATCH}/tp_{1}.bed {1} 20 20 ${PARAM_DEPTH_SCRATCH}/tp_{1}.depth
echo "  False positives..."
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt ${PYTHON} ${PARAM_SCRIPT_PATH}/depth.py ${PARAM_DEPTH_SCRATCH}/reads_{1}.bam ${PARAM_DEPTH_SCRATCH}/fp_{1}.bed {1} 20 20 ${PARAM_DEPTH_SCRATCH}/fp_{1}.depth
echo "  False negatives..."
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt ${PYTHON} ${PARAM_SCRIPT_PATH}/depth.py ${PARAM_DEPTH_SCRATCH}/reads_{1}.bam ${PARAM_DEPTH_SCRATCH}/fn_{1}.bed {1} 20 20 ${PARAM_DEPTH_SCRATCH}/fn_{1}.depth


#####################################################################
# AUGMENT OVERLAP VCFs WITH BAM DEPTH DATA
#####################################################################
PATH_SAMPLE_OVERLAP_TP="${PARAM_RTG_OVERLAP_SCRATCH}/tp-augmented.vcf.gz"
PATH_SAMPLE_OVERLAP_FP="${PARAM_RTG_OVERLAP_SCRATCH}/fp-augmented.vcf.gz"
PATH_SAMPLE_OVERLAP_FN="${PARAM_RTG_OVERLAP_SCRATCH}/fn-augmented.vcf.gz"

echo "  Augmenting VCFs..."
${PYTHON} ${PARAM_SCRIPT_PATH}/add-depth.py ${PATH_SAMPLE_OVERLAP_TP_RAW} ${PARAM_DEPTH_SCRATCH}/tp_ | ${BGZIP} -c > ${PATH_SAMPLE_OVERLAP_TP}
${PYTHON} ${PARAM_SCRIPT_PATH}/add-depth.py ${PATH_SAMPLE_OVERLAP_FP_RAW} ${PARAM_DEPTH_SCRATCH}/fp_ | ${BGZIP} -c > ${PATH_SAMPLE_OVERLAP_FP}
${PYTHON} ${PARAM_SCRIPT_PATH}/add-depth.py ${PATH_SAMPLE_OVERLAP_FN_RAW} ${PARAM_DEPTH_SCRATCH}/fn_ | ${BGZIP} -c > ${PATH_SAMPLE_OVERLAP_FN}

${TABIX} -p vcf ${PATH_SAMPLE_OVERLAP_TP}
${TABIX} -p vcf ${PATH_SAMPLE_OVERLAP_FP}
${TABIX} -p vcf ${PATH_SAMPLE_OVERLAP_FN}


#####################################################################
# CALCULATIONS FOR REPORT
#####################################################################
echo "Performing calculations for report..."

# knitr doesn't play well with building knits outside of its working
# directory.  Currently we get around this with a bit of a kludge, 
# by copying the report files to ${PARAM_KNITR_SCRATCH}, then executing in 
# that directory.
mkdir -p ${PARAM_KNITR_SCRATCH}
cp -f ${PARAM_SCRIPT_PATH}/report.Rnw ${PARAM_KNITR_SCRATCH}
cp -f ${PARAM_SCRIPT_PATH}/report_functions.R ${PARAM_KNITR_SCRATCH}
cp -f ${PARAM_SCRIPT_PATH}/report_calculations.R ${PARAM_KNITR_SCRATCH}
cp -f ${PARAM_SCRIPT_PATH}/test-calcs.R ${PARAM_KNITR_SCRATCH}
cd ${PARAM_KNITR_SCRATCH}

# Run the script.  All options are passed via exported environment 
# variables.  Also save these variables to a file for later source-ing,
# to ease debugging.
export > environment
${RSCRIPT} --vanilla report_calculations.R

if [ ${PARAM_DOTESTS} -eq 1 ]; then
  ${RSCRIPT} --vanilla test-calcs.R
fi

cd ${PARAM_EXEC_PATH}


#####################################################################
# REPORT GENERATION
#####################################################################
echo "Generating report..."

cd ${PARAM_KNITR_SCRATCH}

${RSCRIPT} --vanilla -e "library(knitr); knit('report.Rnw', output = 'report.tex')"

# Remove the report.pdf that may be present in the scratch directory,
# so we can later check whether pdflatex successfully built a report
# or not.
rm -f ${PARAM_KNITR_SCRATCH}/report.pdf

# Run pdflatex
# Latex often 'fails' (returns a nonzero exit status), but still 
# generates a report.  Keep going when this happens, and test for
# failure explicitly later.
pdflatex -interaction nonstopmode report.tex || true
pdflatex -interaction nonstopmode report.tex || true

# Check whether the report.pdf was generated
if [ ! -e ${PARAM_KNITR_SCRATCH}/report.pdf ]; then
	echo >&2 "  Error: pdflatex did not successfully generate report.pdf."
	echo >&2 "  Check ${PARAM_KNITR_SCRATCH}/report.tex and the latex log for errors."
	exit 12
fi

cp ${PARAM_KNITR_SCRATCH}/report.pdf ${PARAM_OUTPUT_PDF_PATH}

echo "Done."
