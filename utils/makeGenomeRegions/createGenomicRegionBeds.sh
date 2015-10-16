#!/bin/sh
set -e -x -o pipefail

RSCRIPT="/usr/bin/env Rscript"
BGZIP="bgzip"
TABIX="tabix"
BEDTOOLS="bedtools"
SORT="sort -k1,1 -k2,2n"
RESULT_DIR="./ensembl_75_regions"
UCSC_TEMP_DIR="./ucsc_temp"

# Access Biomart to get the core regions, and save them as BEDs
${RSCRIPT} getRegionsFromBiomart.R

mkdir -p ${RESULT_DIR}
mkdir -p ${UCSC_TEMP_DIR}

# Sort and compress all BEDs
for file in *.unsorted.bed; do
	${SORT} ${file} | ${BGZIP} > ${RESULT_DIR}/${file%.unsorted.bed}.bed.gz
	${TABIX} ${RESULT_DIR}/${file%.unsorted.bed}.bed.gz
	grep -P '^[0-9]+\t' ${file} | sed 's/^/chr/' | ${SORT} | ${BGZIP} > ${UCSC_TEMP_DIR}/hg19_${file%.unsorted.bed}.bed.gz
	rm ${file}
done

# Generate a hg19 track bed for display on UCSC
rm -f ${RESULT_DIR}/hg19_tracks.bed ${RESULT_DIR}/hg19_tracks.bed.gz
echo 'track name=Genome description=Genome visibility=1' >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
gzip -dc ${UCSC_TEMP_DIR}/hg19_grch37_ensembl.targetgenome.bed.gz >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
echo 'track name=Intergenic description=Intergenic visibility=1' >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
gzip -dc ${UCSC_TEMP_DIR}/hg19_grch37_ensembl.intergenic.bed.gz >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
echo 'track name=Genic description="Coding genes" visibility=1' >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
gzip -dc ${UCSC_TEMP_DIR}/hg19_grch37_ensembl.genes.bed.gz >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
echo 'track name=CDS description="Coding sequence" visibility=1' >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
gzip -dc ${UCSC_TEMP_DIR}/hg19_grch37_ensembl.exonic_coding.bed.gz >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
echo 'track name=3-UTR description="3-UTR" visibility=1' >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
gzip -dc ${UCSC_TEMP_DIR}/hg19_grch37_ensembl.exonic_3utr.bed.gz >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
echo 'track name=5-UTR description="5-UTR" visibility=1' >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
gzip -dc ${UCSC_TEMP_DIR}/hg19_grch37_ensembl.exonic_5utr.bed.gz >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
echo 'track name=Splice description="Essential splice bases" visibility=1' >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
gzip -dc ${UCSC_TEMP_DIR}/hg19_grch37_ensembl.splice.bed.gz >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
echo 'track name=Introns description=Introns visibility=1' >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
gzip -dc ${UCSC_TEMP_DIR}/hg19_grch37_ensembl.intronic.bed.gz >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
echo 'track name=CodingSlop10 description="Coding exons + UTRs with slop 10 into introns" visibility=1' >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
gzip -dc ${UCSC_TEMP_DIR}/hg19_grch37_ensembl.cds_slop10.bed.gz >> ${UCSC_TEMP_DIR}/hg19_tracks.bed
gzip ${UCSC_TEMP_DIR}/hg19_tracks.bed
mv ${UCSC_TEMP_DIR}/hg19_tracks.bed.gz ${RESULT_DIR}
rm -rf ${UCSC_TEMP_DIR}
