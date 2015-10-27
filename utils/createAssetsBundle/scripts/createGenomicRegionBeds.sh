#!/bin/sh
#set -e -x -o pipefail
set -x

RSCRIPT="/usr/bin/env Rscript"
BGZIP="bgzip"
TABIX="tabix"
BEDTOOLS="bedtools"
SORT="sort -k1,1 -k2,2n"
ASSETS_DIR="$1"
TEMP_DIR=$(mktemp -d)

# Access Biomart to get the core regions, and save them as BEDs
${RSCRIPT} getRegionsFromBiomart.R

mkdir -p ${TEMP_DIR}/ucsc
mkdir -p ${TEMP_DIR}/grch37


# Sort and compress all BEDs
for file in *.unsorted.bed; do
	${SORT} ${file} | ${BGZIP} > ${TEMP_DIR}/grch37/${file%.unsorted.bed}.bed.gz
	${TABIX} ${TEMP_DIR}/grch37/${file%.unsorted.bed}.bed.gz
	grep -P '^[0-9]+\t' ${file} | sed 's/^/chr/' | ${SORT} | ${BGZIP} > ${TEMP_DIR}/ucsc/hg19_${file%.unsorted.bed}.bed.gz
	rm ${file}
done


# Generate a hg19 track bed for display on UCSC
echo 'track name=Genome description=Genome visibility=1' >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
gzip -dc ${TEMP_DIR}/ucsc/hg19_grch37_ensembl.targetgenome.bed.gz >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
echo 'track name=Intergenic description=Intergenic visibility=1' >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
gzip -dc ${TEMP_DIR}/ucsc/hg19_grch37_ensembl.intergenic.bed.gz >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
echo 'track name=Genic description="Coding genes" visibility=1' >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
gzip -dc ${TEMP_DIR}/ucsc/hg19_grch37_ensembl.genes.bed.gz >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
echo 'track name=CDS description="Coding sequence" visibility=1' >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
gzip -dc ${TEMP_DIR}/ucsc/hg19_grch37_ensembl.exonic_coding.bed.gz >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
echo 'track name=3-UTR description="3-UTR" visibility=1' >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
gzip -dc ${TEMP_DIR}/ucsc/hg19_grch37_ensembl.exonic_3utr.bed.gz >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
echo 'track name=5-UTR description="5-UTR" visibility=1' >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
gzip -dc ${TEMP_DIR}/ucsc/hg19_grch37_ensembl.exonic_5utr.bed.gz >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
echo 'track name=Splice description="Essential splice bases" visibility=1' >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
gzip -dc ${TEMP_DIR}/ucsc/hg19_grch37_ensembl.splice.bed.gz >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
echo 'track name=Introns description=Introns visibility=1' >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
gzip -dc ${TEMP_DIR}/ucsc/hg19_grch37_ensembl.intronic.bed.gz >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
echo 'track name=CodingSlop10 description="Coding exons + UTRs with slop 10 into introns" visibility=1' >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
gzip -dc ${TEMP_DIR}/ucsc/hg19_grch37_ensembl.cds_slop10.bed.gz >> ${TEMP_DIR}/ucsc/hg19_tracks.bed
gzip ${TEMP_DIR}/ucsc/hg19_tracks.bed


# Place final beds into locations within the assets bundle
mkdir -p ${ASSETS_DIR}/orig/reportable_range
mkdir -p ${ASSETS_DIR}/reportable_range

cp ${TEMP_DIR}/grch37/grch37_ensembl.cds_slop10.bed.gz ${ASSETS_DIR}/reportable_range/target_range.bed.gz
cp ${TEMP_DIR}/grch37/grch37_ensembl.cds_slop10.bed.gz.tbi ${ASSETS_DIR}/reportable_range/target_range.bed.gz.tbi
cp ${TEMP_DIR}/ucsc/hg19_tracks.bed.gz ${ASSETS_DIR}/reportable_range/target_range.bed.gz.tbi

rm -rf ${TEMP_DIR}