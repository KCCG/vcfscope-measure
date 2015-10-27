#!/bin/bash
set -e -x -o pipefail


############################################
## LOCATIONS

RESOURCES_DIR="/directflow/ClinicalGenomicsPipeline/projects/performance-reporter/assets_bundle_generator_resources"
JAVA="java"
BUNDLE_DIR="kccg_performance_reporter_resources_bundle-2.0"


############################################
## PATH MANIPULATION

EXEC_DIR=$(pwd)
BUNDLE_DIR=$(readlink -f ${BUNDLE_DIR})
TEMP_DIR=$(mktemp -d)


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
	cp ${RESOURCES_DIR}/GiaB/2.19/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz ${BUNDLE_DIR}/gold_standard/calls-2.19.vcf.gz
	cp ${RESOURCES_DIR}/GiaB/2.19/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz.tbi ${BUNDLE_DIR}/gold_standard/calls-2.19.vcf.gz.tbi
	cp ${RESOURCES_DIR}/GiaB/2.19/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed.gz ${BUNDLE_DIR}/gold_standard/valid_regions-2.19.bed.gz
)&


############################################
## REAL TIME GENOMICS RTG-TOOLS

# Copy over RTG
cp -r ${RESOURCES_DIR}/RTG ${BUNDLE_DIR}/rtg-tools

# Create an RTG-format genome
(
	${JAVA} -Xmx4G -jar ${BUNDLE_DIR}/orig/RTG/RTG.jar format -f fasta -o ${BUNDLE_DIR}/reference/hs37d5.sdf ${RESOURCES_DIR}/reference/hs37d5.fa.gz
)&


############################################
## REPORTABLE RANGE

# Call createGenomicRegionBeds.sh to fetch Ensembl 75 genomic regions 
# from BioMart and combine into a target reportable range bed.
(
	cd scripts
	source createGenomicRegionBeds.sh ${BUNDLE_DIR}
	gzip -dc ${RESOURCES_DIR}/reference/hs37d5.fa.gz | awk -f scripts/calcChromLengths.awk | grep -E '^([0-9]+|[XY]|MT)' | sed 's/ /\t/g' | sort -k1,1 -k2,2n | bgzip > ${BUNDLE_DIR}/reportable_range/genome.bed.gz
)&


############################################
## REDUNDANT REGIONS

# Currently a mix of UCSC's rmsk track, and LCRs from mDust.  Break these two into
# disjoint subsets (R and L, R and notL, notR and L, notR and notL), to allow
# marginalization later in the performance reporter analysis script.
(
	gzip -dc ${RESOURCES_DIR}/ucsc/rmsk.txt.gz | awk 'BEGIN {OFS="\t"} {print $6, $7, $8}' | grep -E '^chr([0-9]+|[XYM])' | sed 's/^chrM/MT/; s/^chr//' | sort -k1,1 -k2,2n | bedtools merge | bgzip > ${BUNDLE_DIR}/redundant_regions/rmsk.bed.gz
	gzip -dc ${RESOURCES_DIR}/reference/hs37d5.fa.gz | ${RESOURCES_DIR}/tools/mdust -c | awk 'BEGIN {OFS="\t"} {print $1, $3-1, $4}' | grep -E '^([0-9]+|[XY]|MT)' | sort -k1,1 -k2,2n | bedtools merge | bgzip > ${BUNDLE_DIR}/redundant_regions/mdust.bed.gz
)&


############################################
## CREATE BUNDLE ARCHIVE
tar cvf ${BUNDLE_DIR}.tar -C ${BUNDLE_DIR} .


rm -rf ${TEMP_DIR}