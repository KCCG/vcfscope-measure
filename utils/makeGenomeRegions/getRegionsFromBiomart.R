library(biomaRt)
library(GenomicRanges)

# listMarts(host = "grch37.ensembl.org")

# attributePages(ensembl)
# listFilters(ensembl)
# listAttributes(ensembl)
# listAttributes(ensembl, page = "sequences")


fetchEnsemblData = function(ensembl)
{
	getBM(
		attributes = c(
			"ensembl_gene_id", 
			"ensembl_transcript_id", 
			"strand", 
			"rank", 
			"chromosome_name", 
			"start_position", 
			"end_position", 
			"5_utr_start", 
			"5_utr_end", 
			"3_utr_start", 
			"3_utr_end", 
			"transcript_start", 
			"transcript_end", 
			"exon_chrom_start", 
			"exon_chrom_end", 
			"genomic_coding_start", 
			"genomic_coding_end",
			"gene_biotype",
			"transcript_biotype"), 
		filters = c("biotype", "transcript_biotype"), values = list("protein_coding", "protein_coding"), 
		mart = ensembl
	)
}


createRegions = function(data)
{
	data.genes = data[,c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand")]
	data.genes = data.genes[!duplicated(data.genes),]
	ranges.genes = GRanges(
		seqnames = Rle(data.genes$chromosome_name), 
		ranges = IRanges(start = data.genes$start_position, end = data.genes$end_position), 
		strand = Rle(data.genes$strand), 
		ensembl_gene_id = data.genes$ensembl_gene_id)

	ranges.exons = GRanges(
		seqnames = Rle(data$chromosome_name), 
		ranges = IRanges(start = data$exon_chrom_start, end = data$exon_chrom_end), 
		strand = Rle(data$strand), 
		ensembl_gene_id = data$ensembl_gene_id, 
		ensembl_transcript_id = data$ensembl_transcript_id, 
		exon_rank = data$rank)

	data.transcripts = data[,c("ensembl_transcript_id", "chromosome_name", "transcript_start", "transcript_end", "strand")]
	data.transcripts = data.transcripts[!duplicated(data.transcripts),]
	ranges.transcripts = GRanges(
		seqnames = Rle(data.transcripts$chromosome_name), 
		ranges = IRanges(start = data.transcripts$transcript_start, end = data.transcripts$transcript_end), 
		strand = Rle(data.transcripts$strand),
		ensembl_transcript_id = data.transcripts$ensembl_transcript_id)

	data.5utr = data[!is.na(data[,"5_utr_start"]),]
	ranges.5utr = GRanges(
		seqnames = Rle(data.5utr$chromosome_name), 
		ranges = IRanges(start = data.5utr[,"5_utr_start"], end = data.5utr[,"5_utr_end"]), 
		strand = Rle(data.5utr$strand),
		ensembl_transcript_id = data.5utr$ensembl_transcript_id)

	data.3utr = data[!is.na(data[,"3_utr_start"]),]
	ranges.3utr = GRanges(
		seqnames = Rle(data.3utr$chromosome_name), 
		ranges = IRanges(start = data.3utr[,"3_utr_start"], end = data.3utr[,"3_utr_end"]), 
		strand = Rle(data.3utr$strand),
		ensembl_transcript_id = data.3utr$ensembl_transcript_id)

	ranges.exon_starts = GRanges(
		seqnames = Rle(data$chromosome_name), 
		ranges = IRanges(start = data$exon_chrom_start, end = data$exon_chrom_start), 
		strand = Rle(data$strand), 
		ensembl_gene_id = data$ensembl_gene_id, 
		ensembl_transcript_id = data$ensembl_transcript_id, 
		exon_rank = data$rank)
	ranges.exon_ends = GRanges(
		seqnames = Rle(data$chromosome_name), 
		ranges = IRanges(start = data$exon_chrom_end, end = data$exon_chrom_end), 
		strand = Rle(data$strand), 
		ensembl_gene_id = data$ensembl_gene_id, 
		ensembl_transcript_id = data$ensembl_transcript_id, 
		exon_rank = data$rank)
	ranges.transcript_starts = GRanges(
		seqnames = Rle(data.transcripts$chromosome_name), 
		ranges = IRanges(start = data.transcripts$transcript_start, end = data.transcripts$transcript_start), 
		strand = Rle(data.transcripts$strand),
		ensembl_transcript_id = data.transcripts$ensembl_transcript_id)
	ranges.transcript_ends = GRanges(
		seqnames = Rle(data.transcripts$chromosome_name), 
		ranges = IRanges(start = data.transcripts$transcript_end, end = data.transcripts$transcript_end), 
		strand = Rle(data.transcripts$strand),
		ensembl_transcript_id = data.transcripts$ensembl_transcript_id)
	ranges.exon_starts.internal = setdiff(ranges.exon_starts, ranges.transcript_starts)
	ranges.exon_ends.internal = setdiff(ranges.exon_ends, ranges.transcript_ends)
	ranges.splice_site_starts = GRanges(
		seqnames = seqnames(ranges.exon_starts.internal),
		ranges = IRanges(start = start(ranges.exon_starts.internal) - 2, end = start(ranges.exon_starts.internal) - 1),
		strand = strand(ranges.exon_starts.internal))
	ranges.splice_site_ends = GRanges(
		seqnames = seqnames(ranges.exon_ends.internal),
		ranges = IRanges(start = start(ranges.exon_ends.internal) + 1, end = start(ranges.exon_ends.internal) + 2),
		strand = strand(ranges.exon_ends.internal))
	ranges.splice_sites = union(ranges.splice_site_starts, ranges.splice_site_ends)

	data.coding = data[!is.na(data$genomic_coding_start),]
	ranges.coding = GRanges(
		seqnames = Rle(data.coding$chromosome_name), 
		ranges = IRanges(start = data.coding$genomic_coding_start, end = data.coding$genomic_coding_end), 
		strand = Rle(data.coding$strand), 
		ensembl_gene_id = data.coding$ensembl_gene_id, 
		ensembl_transcript_id = data.coding$ensembl_transcript_id, 
		exon_rank = data.coding$rank)

	ranges.introns = setdiff(ranges.transcripts, ranges.exons)

	list(
		genes = ranges.genes,
		transcripts = ranges.transcripts,
		exons = ranges.exons,
		introns = ranges.introns,
		coding = ranges.coding,
		utr5 = ranges.5utr,
		utr3 = ranges.3utr,
		splice_sites = ranges.splice_sites
		)
}


writeBed = function(ranges, path)
{
	data = data.frame(
		chrom = as.vector(seqnames(ranges)),
		chromStart = as.integer(as.vector(start(ranges))) - 1,
		chromEnd = as.integer(as.vector(end(ranges))))
	scipen = options()$scipen
	options(scipen = 999)
	write.table(data, file = path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
	options(scipen = scipen)
}


ensembl = useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

data = fetchEnsemblData(ensembl)

regions = createRegions(data)

regions.exonic.coding = reduce(regions$coding)
regions.exonic.utr = reduce(union(regions$utr5, regions$utr3))
regions.splice = reduce(regions$splice_sites)
regions.intronic = reduce(regions$introns)
regions.genes = reduce(regions$genes)

writeBed(regions.exonic.coding, "grch37_ensembl.exonic_coding.bed")
writeBed(regions.exonic.utr, "grch37_ensembl.exonic_utr.bed")
writeBed(regions.splice, "grch37_ensembl.splice.bed")
writeBed(regions.intronic, "grch37_ensembl.intronic.bed")
writeBed(regions.genes, "grch37_ensembl.genes.bed")
