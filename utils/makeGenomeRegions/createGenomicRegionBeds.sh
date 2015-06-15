#!/bin/sh
set -e

RSCRIPT="/usr/bin/env Rscript"
BGZIP="/home/marpin/software/htslib/bgzip"
TABIX="/home/marpin/software/htslib/tabix"
BEDTOOLS="/home/marpin/software/bedtools2/bin/bedtools"
SORT="sort -k1,1 -k2,2n"
RESULT_DIR="./ensembl_75_regions"

${RSCRIPT} getRegionsFromBiomart.R

${SORT} grch37_ensembl.exonic_coding.bed | ${BGZIP} -c > grch37_ensembl.exonic_coding.bed.gz
${SORT} grch37_ensembl.exonic_utr.bed | ${BGZIP} -c > grch37_ensembl.exonic_utr.bed.gz
${SORT} grch37_ensembl.splice.bed | ${BGZIP} -c > grch37_ensembl.splice.bed.gz
${SORT} grch37_ensembl.intronic.bed | ${BGZIP} -c > grch37_ensembl.intronic.bed.gz
${SORT} grch37_ensembl.genes.bed | ${BGZIP} -c > grch37_ensembl.genes.bed.gz

${TABIX} grch37_ensembl.exonic_coding.bed.gz
${TABIX} grch37_ensembl.exonic_utr.bed.gz
${TABIX} grch37_ensembl.splice.bed.gz
${TABIX} grch37_ensembl.intronic.bed.gz
${TABIX} grch37_ensembl.genes.bed.gz

mkdir -p ${RESULT_DIR}

if [ ! -e hs37d5.fa.gz.fai ]; then
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai
fi
# if [ ! -e hs37d5.fa.gz ]; then
# 	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
# fi
# gzip -dc hs37d5.fa.gz | awk -f getChromLengths.awk

awk 'BEGIN { FS="\t"; OFS="\t" } { print $1, 0, $2-1 }' < hs37d5.fa.gz.fai | ${SORT} | ${BGZIP} -c > hs37d5_genome.bed.gz
${TABIX} hs37d5_genome.bed.gz

mv *.bed.gz ${RESULT_DIR}
mv *.bed.gz.tbi ${RESULT_DIR}
rm *.bed

${BEDTOOLS} subtract -a ${RESULT_DIR}/hs37d5_genome.bed.gz -b ${RESULT_DIR}/grch37_ensembl.genes.bed.gz | ${BGZIP} -c > ${RESULT_DIR}/grch37_ensembl.intergenic.bed.gz
${TABIX} ${RESULT_DIR}/grch37_ensembl.intergenic.bed.gz
