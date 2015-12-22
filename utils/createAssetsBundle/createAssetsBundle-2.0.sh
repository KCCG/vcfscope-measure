#!/bin/bash
#set -e -x -o pipefail


############################################
## LOCATIONS

export RESOURCES_DIR="/directflow/ClinicalGenomicsPipeline/projects/performance-reporter/assets_bundle_generator_resources"
export JAVA="java"
BUNDLE_DIR="kccg_performance_reporter_resources_bundle-2.0"


############################################
## PATH MANIPULATION

export EXEC_DIR=$(pwd)
export BUNDLE_DIR=$(readlink -f ${BUNDLE_DIR})
export TEMP_DIR=$(mktemp -d)


############################################
## PREPARE BUNDLE DIRECTORY

mkdir -p ${BUNDLE_DIR}

mkdir -p ${BUNDLE_DIR}/gold_standard
mkdir -p ${BUNDLE_DIR}/reportable_range
mkdir -p ${BUNDLE_DIR}/redundant_regions
mkdir -p ${BUNDLE_DIR}/reference
mkdir -p ${BUNDLE_DIR}/rtg-tools


############################################
## GENOME IN A BOTTLE STANDARD

# Copy over the GiaB standard
(
	cp ${RESOURCES_DIR}/GiaB/2.19/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz ${BUNDLE_DIR}/gold_standard/calls-2.19.vcf.gz && \
	cp ${RESOURCES_DIR}/GiaB/2.19/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz.tbi ${BUNDLE_DIR}/gold_standard/calls-2.19.vcf.gz.tbi && \
	cp ${RESOURCES_DIR}/GiaB/2.19/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed.gz ${BUNDLE_DIR}/gold_standard/valid_regions-2.19.bed.gz
)&
PID_GIABCOPY=$!


############################################
## REAL TIME GENOMICS RTG-TOOLS

# Copy over RTG
cp -r ${RESOURCES_DIR}/RTG/* ${BUNDLE_DIR}/rtg-tools

# Create an RTG-format genome
(
	if [[ ! -e ${BUNDLE_DIR}/reference/hs37d5.sdf ]]; then
		${JAVA} -Xmx4G -jar ${BUNDLE_DIR}/rtg-tools/RTG.jar format -f fasta -o ${BUNDLE_DIR}/reference/hs37d5.sdf ${RESOURCES_DIR}/reference/hs37d5.fa.gz
	fi
)&
PID_RTGREF=$!


############################################
## REDUNDANT REGIONS

# Currently only LCRs from mDust.
if [[ ! -e ${BUNDLE_DIR}/redundant_regions/mdust.bed.gz ]]; then
	gzip -dc ${RESOURCES_DIR}/reference/hs37d5.fa.gz | ${RESOURCES_DIR}/tools/mdust -c | awk 'BEGIN {OFS="\t"} {print $1, $3-1, $4}' | grep -E '^([0-9]+|[XY]|MT)' | sort -k1,1 -k2,2n | bedtools merge | bgzip > ${BUNDLE_DIR}/redundant_regions/mdust.bed.gz
fi


############################################
## REPORTABLE RANGE

# Call createGenomicRegionBeds.sh to fetch Ensembl 75 genomic regions 
# from BioMart and combine into a target reportable range bed, less
# the mdust LCRs.
if [[ ! -e ${BUNDLE_DIR}/reportable_range/genome.bed.gz ]]; then
	cd scripts
	source createGenomicRegionBeds.sh ${BUNDLE_DIR}/redundant_regions/mdust.bed.gz ${BUNDLE_DIR} && \
	gzip -dc ${RESOURCES_DIR}/reference/hs37d5.fa.gz | awk -f calcChromLengths.awk | grep -E '^([0-9]+|[XY]|MT)' | sed 's/ /\t/g' | sort -k1,1 -k2,2n | bgzip > ${BUNDLE_DIR}/reportable_range/genome.bed.gz
fi

############################################
## CREATE BUNDLE ARCHIVE
wait $PID_GIABCOPY $PID_RTGREF
tar cvf ${BUNDLE_DIR}.tar -C ${BUNDLE_DIR} .


rm -rf ${TEMP_DIR}
