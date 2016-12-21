#!/bin/bash
set -e -x -u -o pipefail
set +x
IFS=$'\n\t'

function message {
  echo >&2 "$(date -u)" "$1"
}


#####################################################################
# VERSION
#####################################################################
export CONST_VERSION_SCRIPT="2.1.0"


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

  RSCRIPT="/usr/bin/Rscript"
  R="/usr/bin/R"
  JAVA=`which java`
  PYTHON=`which python`
  BEDTOOLS=`which bedtools`
  TABIX=`which tabix`
  BGZIP=`which bgzip`
  BCFTOOLS=`which bcftools`
  PARALLEL="`which parallel` --gnu --halt now,fail=1"
  SAMTOOLS=`which samtools`

  RTG_TOOLS="${PATH_RESOURCES_HEAD}/rtg-tools/RTG.jar"
  RTG_THREADS=`nproc`
  mem_in_mb=`head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}'`                 # Calculate 80% of memory size, for java
  RTG_VCFEVAL="${JAVA} -Xmx${mem_in_mb}m -jar ${RTG_TOOLS} vcfeval -T ${RTG_THREADS}"
else
  # Wolfpack settings (marpin only for now)
  PATH_RESOURCES_HEAD="/directflow/ClinicalGenomicsPipeline/projects/vcfscope-reporter/vcfscope_reporter_resources_bundle-2.0"
  PATH_SCRATCH_DEFAULT="/directflow/ClinicalGenomicsPipeline/tmp"

  RSCRIPT="/home/marpin/bin/Rscript"
  R="/home/marpin/bin/R"
  JAVA="/usr/java/latest/bin/java"
  PYTHON=`which python`
  BEDTOOLS="/home/marpin/software/bedtools2/bin/bedtools"
  TABIX="/home/marpin/software/htslib/tabix"
  BGZIP="/home/marpin/software/htslib/bgzip"
  BCFTOOLS="/home/marpin/software/bcftools/bcftools"
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
export CONST_REFERENCE_BSGENOME="BSgenome.HSapiens.1000g.37d5"    # This is a custom package, available at /share/ClusterShare/biodata/contrib/marpin/reference/hs37d5/build/BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz

# Script location
export PARAM_SCRIPT_PATH=$(readlink -f $(dirname $0))
export PARAM_EXEC_PATH=$(pwd)

# Scratch space
mkdir -p ${PATH_SCRATCH_DEFAULT}
export PARAM_SCRATCH=$(mktemp -d --tmpdir=${PATH_SCRATCH_DEFAULT} vcfscope.XXXXXXXXXX)
export PARAM_INPUT_SCRATCH="${PARAM_SCRATCH}/input"
export PARAM_RTG_OVERLAP_SCRATCH="${PARAM_SCRATCH}/overlap"
export PARAM_R_SCRATCH="${PARAM_SCRATCH}/R"
export PARAM_DEPTH_SCRATCH="${PARAM_SCRATCH}/depth"

# Program parameters
export PARAM_INPUT_VCFGZ_PATH
export PARAM_INPUT_BAM_PATH
export PARAM_REGION_BED_SUPPLIED
export PARAM_REGION_BED_PATH
export PARAM_OUTPUT_RDS_PATH
export PARAM_VERSION_EXEC_HOST
export PARAM_VERSION_RTG
export PARAM_VERSION_JAVA
export PARAM_VERSION_BEDTOOLS
export PARAM_VERSION_SAMTOOLS
export PARAM_STORE_ALL_VARIANTS

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
Usage: ${0##*/} [-r BEDFILE] <VCF> <BAM> <RDSOUT>

Perform WGS performance calculations on calls for a single sample.

    VCF          Input NA12878 genotype calls, in vcf.gz format.
    BAM          Mapped NA12878 reads used to generate the VCF calls.
                 An accompanying index (<BAM>.bai) must be present in 
                 the same directory.
    RDSOUT       Write the output performance data to this file.
    -r BEDFILE   Restrict analysis to the regions in BEDFILE only.
                 Default: the full genome is considered.
    -a           Store all region variants in the output RDS.
    -h           Display this help and exit.

Version ${CONST_VERSION_SCRIPT}

Mark Pinese
EOF
}

OPTIND=1
PARAM_INPUT_VCFGZ_PATH=""
PARAM_REGION_BED_SUPPLIED=0
PARAM_REGION_BED_PATH="NA"
PARAM_STORE_ALL_VARIANTS=0

while getopts "r:ha" opt; do
  case "$opt" in
    h)
      print_usage
      exit 0
      ;;
    a)
      PARAM_STORE_ALL_VARIANTS=1
      ;;
    r)
      PARAM_REGION_BED_SUPPLIED=1
      PARAM_REGION_BED_PATH=$(readlink -f "${OPTARG}")
      ;;
    '?')
      print_usage >&2
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))

if [ $# -ne 3 ]; then
  print_usage >&2
  exit 1
fi

PARAM_INPUT_VCFGZ_PATH=$(readlink -f $1)
PARAM_INPUT_BAM_PATH=$(readlink -f $2)
PARAM_OUTPUT_RDS_PATH=$(readlink -f $3)

#####################################################################
# PARAMETER CHECKING
#####################################################################

if [ ! -d "${PATH_RESOURCES_HEAD}" ]; then
  message "Error: Resources path ${PATH_RESOURCES_HEAD} not found."
  exit 2
fi

if [ ! -e ${PARAM_INPUT_VCFGZ_PATH} ]; then
  message "Error: Input file ${PARAM_INPUT_VCFGZ_PATH} not found."
  exit 3
fi

if [ ! -e ${PARAM_INPUT_BAM_PATH} ]; then
  message "Error: Input file ${PARAM_INPUT_BAM_PATH} not found."
  exit 3
fi

if [ ! -e ${PARAM_INPUT_BAM_PATH}.bai ]; then
  message "Error: Index not found for input BAM file (expected ${PARAM_INPUT_BAM_PATH}.bai)"
  exit 3
fi

if [ -e ${PARAM_OUTPUT_RDS_PATH} ]; then
  message "Error: Output file ${PARAM_OUTPUT_RDS_PATH} already exists."
  exit 4
fi

if [ ${PARAM_REGION_BED_SUPPLIED} -eq 1 ] && [ ! -e ${PARAM_REGION_BED_PATH} ]; then
  message "Error: Region file ${PARAM_REGION_BED_PATH} not found."
  exit 5
fi


#####################################################################
# SOFTWARE CHECKING
#####################################################################
if [ ! -e ${RSCRIPT} ]; then
  message "Error: Rscript executable not found at ${RSCRIPT}."
  exit 6
fi

if [ ! -e ${JAVA} ]; then
  message "Error: Java executable not found at ${JAVA}."
  exit 7
fi

if [ ! -e ${RTG_TOOLS} ]; then
  message "Error: RTG.jar not found at ${RTG_TOOLS}."
  exit 8
fi

if [ ! -f ${PARAM_SCRIPT_PATH}/measure_functions.R -o ! -f ${PARAM_SCRIPT_PATH}/measure_calculations.R ]; then
  message "Error: Missing at least one required R source file."
  exit 9
fi


#####################################################################
# CHECK THE INPUT FILE IS A SINGLE-SAMPLE VCF
#####################################################################
if [ $(gzip -dc ${PARAM_INPUT_VCFGZ_PATH} | grep -m1 '^#[^#]' | tr '\t' '\n' | wc -l) -ne 10 ]; then
  message "Error: Multi-sample VCF found.  Input VCF must contain variants for only a single sample."
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
PARAM_VERSION_SAMTOOLS=$(${SAMTOOLS} --version | head -n 1 | cut -d' ' -f 2)


#####################################################################
# PARAMETER LOGGING
#####################################################################
message "performance_report.sh parameters:"
message "  PATH_RESOURCES_HEAD=${PATH_RESOURCES_HEAD}"
message "  PATH_SCRATCH_DEFAULT=${PATH_SCRATCH_DEFAULT}"
message "  BEDTOOLS=${BEDTOOLS}"
message "  BCFTOOLS=${BCFTOOLS}"
message "  BGZIP=${BGZIP}"
message "  JAVA=${JAVA}"
message "  PARALLEL=${PARALLEL}"
message "  PYTHON=${PYTHON}"
message "  R=${R}"
message "  RSCRIPT=${RSCRIPT}"
message "  RTG_TOOLS=${RTG_TOOLS}"
message "  RTG_VCFEVAL=${RTG_VCFEVAL}"
message "  SAMTOOLS=${SAMTOOLS}"
message "  TABIX=${TABIX}"
message "  RTG_THREADS=${RTG_THREADS}"
message "  mem_in_mb=${mem_in_mb}"
message "  CONST_GOLD_CALLS_VCFGZ=${CONST_GOLD_CALLS_VCFGZ}"
message "  CONST_GOLD_CALLS_VCFGZTBI=${CONST_GOLD_CALLS_VCFGZTBI}"
message "  CONST_GOLD_HARDMASK_VALID_REGIONS_BEDGZ=${CONST_GOLD_HARDMASK_VALID_REGIONS_BEDGZ}"
message "  CONST_REFERENCE_SDF=${CONST_REFERENCE_SDF}"
message "  CONST_GENOME_BEDGZ=${CONST_GENOME_BEDGZ}"
message "  CONST_REFERENCE_BSGENOME=${CONST_REFERENCE_BSGENOME}"
message "  PARAM_SCRIPT_PATH=${PARAM_SCRIPT_PATH}"
message "  PARAM_SCRATCH=${PARAM_SCRATCH}"
message "  PARAM_INPUT_SCRATCH=${PARAM_INPUT_SCRATCH}"
message "  PARAM_RTG_OVERLAP_SCRATCH=${PARAM_RTG_OVERLAP_SCRATCH}"
message "  PARAM_R_SCRATCH=${PARAM_R_SCRATCH}"
message "  PARAM_INPUT_VCFGZ_PATH=${PARAM_INPUT_VCFGZ_PATH}"
message "  PARAM_INPUT_BAM_PATH=${PARAM_INPUT_BAM_PATH}"
message "  PARAM_REGION_BED_SUPPLIED=${PARAM_REGION_BED_SUPPLIED}"
message "  PARAM_REGION_BED_PATH=${PARAM_REGION_BED_PATH}"
message "  PARAM_OUTPUT_RDS_PATH=${PARAM_OUTPUT_RDS_PATH}"
message "  PARAM_STORE_ALL_VARIANTS=${PARAM_STORE_ALL_VARIANTS}"
message "  PARAM_VERSION_EXEC_HOST=${PARAM_VERSION_EXEC_HOST}"
message "  PARAM_VERSION_RTG=${PARAM_VERSION_RTG}"
message "  PARAM_VERSION_JAVA=${PARAM_VERSION_JAVA}"
message "  PARAM_VERSION_BEDTOOLS=${PARAM_VERSION_BEDTOOLS}"
message "  PARAM_VERSION_SAMTOOLS=${PARAM_VERSION_SAMTOOLS}"
message "  PATH_TEST_VARIANTS=${PATH_TEST_VARIANTS}"
message "  PATH_TEST_VARIANTS_INDEX=${PATH_TEST_VARIANTS_INDEX}"
message "  PATH_TEST_READS=${PATH_TEST_READS}"
message "  PATH_TEST_READS_INDEX=${PATH_TEST_READS_INDEX}"
message "  PATH_GOLD_VARIANTS=${PATH_GOLD_VARIANTS}"
message "  PATH_GOLD_VARIANTS_INDEX=${PATH_GOLD_VARIANTS_INDEX}"
message "  PATH_GOLD_REGIONS=${PATH_GOLD_REGIONS}"


#####################################################################
# CREATE TEMPORARY DIRECTORIES
#####################################################################
mkdir -p ${PARAM_SCRATCH}
mkdir -p ${PARAM_INPUT_SCRATCH}
mkdir -p ${PARAM_DEPTH_SCRATCH}


#####################################################################
# CREATE OUTPUT DIRECTORIES
#####################################################################
mkdir -p $(dirname ${PARAM_OUTPUT_RDS_PATH})


#####################################################################
# SUBSET TO REGION BED
#####################################################################
if [ ${PARAM_REGION_BED_SUPPLIED} -eq 1 ]; then
  message "Subsetting input files to supplied BED..."
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
message "Computing VCF overlaps..."
if [ -e ${PARAM_RTG_OVERLAP_SCRATCH} ]; then
  message "PARAM_RTG_OVERLAP_SCRATCH directory ${PARAM_RTG_OVERLAP_SCRATCH} already exists.  Clearing scratch directory and continuing..."
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
  chown -R dnanexus:dnanexus $(dirname ${PARAM_OUTPUT_RDS_PATH})
fi
# DIRTY DIRTY DIRTY
#####################################################################

# TODO: Parse ${PARAM_RTG_OVERLAP_SCRATCH}/vcfeval.log to identify regions to exclude
# eg Evaluation too complex (5001 unresolved paths, 18033 iterations) at reference region 2:105849275-105849281. Variants in this region will not be included in results.
# This will be required to remove the TNs in this region, but it's polish.


#####################################################################
# CALCULATE BAM DEPTH AT VARIANT LOCATIONS
#####################################################################
message "Calculating depth at variant loci..."
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
message "  Sharding BAM..."
# Process chromosomes in descending order of size, to optimize core usage
${SAMTOOLS} view -H ${PATH_TEST_READS} | grep '^@SQ' | cut -f 2,3 | sed 's/[SL]N://g' | sort -k2,2rn | cut -f 1 > ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt ${SAMTOOLS} view -1 -o ${PARAM_DEPTH_SCRATCH}/reads_{1}.bam ${PATH_TEST_READS} {1}
message "  Indexing shards..."
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt ${SAMTOOLS} index ${PARAM_DEPTH_SCRATCH}/reads_{1}.bam

# Fill in empty beds for chromosomes in the BAM but not the beds
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt touch ${PARAM_DEPTH_SCRATCH}/tp_{1}.bed
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt touch ${PARAM_DEPTH_SCRATCH}/fp_{1}.bed
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt touch ${PARAM_DEPTH_SCRATCH}/fn_{1}.bed

# For each shard (chromosome), calculate depth in the BAM in the VCF regions
message "  True positives..."
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt ${PYTHON} ${PARAM_SCRIPT_PATH}/depth.py ${PARAM_DEPTH_SCRATCH}/reads_{1}.bam ${PARAM_DEPTH_SCRATCH}/tp_{1}.bed {1} 20 20 ${PARAM_DEPTH_SCRATCH}/tp_{1}.depth
message "  False positives..."
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt ${PYTHON} ${PARAM_SCRIPT_PATH}/depth.py ${PARAM_DEPTH_SCRATCH}/reads_{1}.bam ${PARAM_DEPTH_SCRATCH}/fp_{1}.bed {1} 20 20 ${PARAM_DEPTH_SCRATCH}/fp_{1}.depth
message "  False negatives..."
eval ${PARALLEL} -a ${PARAM_DEPTH_SCRATCH}/bam_chroms.txt ${PYTHON} ${PARAM_SCRIPT_PATH}/depth.py ${PARAM_DEPTH_SCRATCH}/reads_{1}.bam ${PARAM_DEPTH_SCRATCH}/fn_{1}.bed {1} 20 20 ${PARAM_DEPTH_SCRATCH}/fn_{1}.depth


#####################################################################
# AUGMENT OVERLAP VCFs WITH BAM DEPTH DATA
#####################################################################
PATH_SAMPLE_OVERLAP_TP="${PARAM_RTG_OVERLAP_SCRATCH}/tp-augmented.vcf.gz"
PATH_SAMPLE_OVERLAP_FP="${PARAM_RTG_OVERLAP_SCRATCH}/fp-augmented.vcf.gz"
PATH_SAMPLE_OVERLAP_FN="${PARAM_RTG_OVERLAP_SCRATCH}/fn-augmented.vcf.gz"

message "  Augmenting VCFs..."
${PYTHON} ${PARAM_SCRIPT_PATH}/add-depth.py ${PATH_SAMPLE_OVERLAP_TP_RAW} ${PARAM_DEPTH_SCRATCH}/tp_ | ${BGZIP} -c > ${PATH_SAMPLE_OVERLAP_TP}
${PYTHON} ${PARAM_SCRIPT_PATH}/add-depth.py ${PATH_SAMPLE_OVERLAP_FP_RAW} ${PARAM_DEPTH_SCRATCH}/fp_ | ${BGZIP} -c > ${PATH_SAMPLE_OVERLAP_FP}
${PYTHON} ${PARAM_SCRIPT_PATH}/add-depth.py ${PATH_SAMPLE_OVERLAP_FN_RAW} ${PARAM_DEPTH_SCRATCH}/fn_ | ${BGZIP} -c > ${PATH_SAMPLE_OVERLAP_FN}

${TABIX} -p vcf ${PATH_SAMPLE_OVERLAP_TP}
${TABIX} -p vcf ${PATH_SAMPLE_OVERLAP_FP}
${TABIX} -p vcf ${PATH_SAMPLE_OVERLAP_FN}


#####################################################################
# CALCULATIONS FOR REPORT
#####################################################################
message "Performing calculations for report..."

mkdir -p ${PARAM_R_SCRATCH}
cp -f ${PARAM_SCRIPT_PATH}/measure_functions.R ${PARAM_R_SCRATCH}
cp -f ${PARAM_SCRIPT_PATH}/measure_calculations.R ${PARAM_R_SCRATCH}
cd ${PARAM_R_SCRATCH}

# Run the script.  All options are passed via exported environment 
# variables.  Also save these variables to a file for later source-ing,
# to ease debugging.
export > environment
${RSCRIPT} --vanilla measure_calculations.R

message "Done."
