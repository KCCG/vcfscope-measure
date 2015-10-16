library(biomaRt)
library(GenomicRanges)
library(BSgenome)

# listMarts(host = "grch37.ensembl.org")

# attributePages(ensembl)
# listFilters(ensembl)
# listAttributes(ensembl)
# listAttributes(ensembl, page = "sequences")


createAspirationalReportableRangeRegions = function(data, seqinfo, slop_delta = 10)
{
	# Create the following ranges: 
	#   Ensembl 75 coding exons + UTRs + 10 bp into introns - mitochondrial genes - non-canonical chromosomes (ie patch contigs etc)
	# These are our aspirational reportable range for version 1.0 of the clinical pipeline.
	# The actual reportable range will be somewhat smaller, due to some regions not being covered
	# by empirical sequencing data, but the logic to reduce these regions appropriately based on
	# observed depth is in another file.

	# Rank is the exon number in the transcript, in mRNA order.
	# Use this to define the maximum rank (ie the number of exons
	# in each mRNA)
	data$maxrank = tapply(data$rank, data$ensembl_transcript_id, max)[data$ensembl_transcript_id]

	# Using rank (exon number) and maxrank (number of exons), assign
	# each exon a type, depending on whether it is first (the first
	# exon, in genomic coordinates, of a multi-exon mRNA), last 
	# (the last exon, in genomic coordinates, of a multi-exon mRNA),
	# internal (an internal exon of a multi-exon mRNA), or only 
	# (the only exon of an unspliced mRNA).
	data$exon_type = NA
	data$exon_type[data$rank == 1 & data$maxrank > 1 & data$strand == 1] = "first"
	data$exon_type[data$rank == 1 & data$maxrank > 1 & data$strand == -1] = "last"
	data$exon_type[data$rank == data$maxrank & data$maxrank > 1 & data$strand == 1] = "last"
	data$exon_type[data$rank == data$maxrank & data$maxrank > 1 & data$strand == -1] = "first"
	data$exon_type[data$maxrank == 1] = "only"
	data$exon_type[data$rank > 1 & data$rank < data$maxrank] = "internal"
	data$exon_type = factor(data$exon_type)

	# Define exon coordinates with internal slop (into the introns).
	# The slop applies to start and end coordinates depending on the
	# exon type determined above, as:
	#   Type        Start slop    End slop
	#   first       0             +delta
	#   last        -delta        0
	#   only        0             0
	#   internal    -delta        +delta
	# Where delta is the number of slop bases (here defined as 10)

	data$exon_chrom_start_slop = NA
	data$exon_chrom_start_slop[data$exon_type == "first"] =    data$exon_chrom_start[data$exon_type == "first"]
	data$exon_chrom_start_slop[data$exon_type == "last"] =     data$exon_chrom_start[data$exon_type == "last"] - slop_delta
	data$exon_chrom_start_slop[data$exon_type == "only"] =     data$exon_chrom_start[data$exon_type == "only"]
	data$exon_chrom_start_slop[data$exon_type == "internal"] = data$exon_chrom_start[data$exon_type == "internal"] - slop_delta

	data$exon_chrom_end_slop = NA
	data$exon_chrom_end_slop[data$exon_type == "first"] =      data$exon_chrom_end[data$exon_type == "first"] + slop_delta
	data$exon_chrom_end_slop[data$exon_type == "last"] =       data$exon_chrom_end[data$exon_type == "last"]
	data$exon_chrom_end_slop[data$exon_type == "only"] =       data$exon_chrom_end[data$exon_type == "only"]
	data$exon_chrom_end_slop[data$exon_type == "internal"] =   data$exon_chrom_end[data$exon_type == "internal"] + slop_delta

	# Keep only entries in the autosomes and sex chromosomes.  In particular,
	# exclude all mtDNA and patch contigs.
	data = data[as.character(data$chromosome_name) %in% as.character(c(1:22, "X", "Y")),]

	# Transfer the slopped exon ranges into a GRanges object, for reduction and trimming.
	ranges.exonic_slopped = trim(GRanges(
		seqnames = Rle(data$chromosome_name), 
		ranges = IRanges(start = data$exon_chrom_start_slop, end = data$exon_chrom_end_slop), 
		strand = Rle(data$strand), 
		ensembl_gene_id = data$ensembl_gene_id, 
		ensembl_transcript_id = data$ensembl_transcript_id, 
		exon_rank = data$rank,
		seqinfo = seqinfo))

	# Simplify overlapping and abutting ranges
	ranges.exonic_slopped_reduced = reduce(ranges.exonic_slopped)

	# Return the final reportable range object
	return(ranges.exonic_slopped_reduced)
}


createFunctionalRegions = function(data, seqinfo)
{
	# Keep only entries in the autosomes, sex chromosomes, and MT.  This ensures
	# overlap with the hs37d5 reference in seqinfo.
	data = data[as.character(data$chromosome_name) %in% as.character(c(1:22, "X", "Y", "MT")),]

	data.genes = data[,c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand")]
	data.genes = data.genes[!duplicated(data.genes),]
	ranges.genes = GRanges(
		seqnames = Rle(data.genes$chromosome_name), 
		ranges = IRanges(start = data.genes$start_position, end = data.genes$end_position), 
		strand = Rle(data.genes$strand), 
		ensembl_gene_id = data.genes$ensembl_gene_id,
		seqinfo = seqinfo)

	ranges.exons = GRanges(
		seqnames = Rle(data$chromosome_name), 
		ranges = IRanges(start = data$exon_chrom_start, end = data$exon_chrom_end), 
		strand = Rle(data$strand), 
		ensembl_gene_id = data$ensembl_gene_id, 
		ensembl_transcript_id = data$ensembl_transcript_id, 
		exon_rank = data$rank,
		seqinfo = seqinfo)

	data.transcripts = data[,c("ensembl_transcript_id", "chromosome_name", "transcript_start", "transcript_end", "strand")]
	data.transcripts = data.transcripts[!duplicated(data.transcripts),]
	ranges.transcripts = GRanges(
		seqnames = Rle(data.transcripts$chromosome_name), 
		ranges = IRanges(start = data.transcripts$transcript_start, end = data.transcripts$transcript_end), 
		strand = Rle(data.transcripts$strand),
		ensembl_transcript_id = data.transcripts$ensembl_transcript_id,
		seqinfo = seqinfo)

	data.5utr = data[!is.na(data[,"5_utr_start"]),]
	ranges.5utr = GRanges(
		seqnames = Rle(data.5utr$chromosome_name), 
		ranges = IRanges(start = data.5utr[,"5_utr_start"], end = data.5utr[,"5_utr_end"]), 
		strand = Rle(data.5utr$strand),
		ensembl_transcript_id = data.5utr$ensembl_transcript_id,
		seqinfo = seqinfo)

	data.3utr = data[!is.na(data[,"3_utr_start"]),]
	ranges.3utr = GRanges(
		seqnames = Rle(data.3utr$chromosome_name), 
		ranges = IRanges(start = data.3utr[,"3_utr_start"], end = data.3utr[,"3_utr_end"]), 
		strand = Rle(data.3utr$strand),
		ensembl_transcript_id = data.3utr$ensembl_transcript_id,
		seqinfo = seqinfo)

	ranges.exon_starts = GRanges(
		seqnames = Rle(data$chromosome_name), 
		ranges = IRanges(start = data$exon_chrom_start, end = data$exon_chrom_start), 
		strand = Rle(data$strand), 
		ensembl_gene_id = data$ensembl_gene_id, 
		ensembl_transcript_id = data$ensembl_transcript_id, 
		exon_rank = data$rank,
		seqinfo = seqinfo)

	ranges.exon_ends = GRanges(
		seqnames = Rle(data$chromosome_name), 
		ranges = IRanges(start = data$exon_chrom_end, end = data$exon_chrom_end), 
		strand = Rle(data$strand), 
		ensembl_gene_id = data$ensembl_gene_id, 
		ensembl_transcript_id = data$ensembl_transcript_id, 
		exon_rank = data$rank,
		seqinfo = seqinfo)

	ranges.transcript_starts = GRanges(
		seqnames = Rle(data.transcripts$chromosome_name), 
		ranges = IRanges(start = data.transcripts$transcript_start, end = data.transcripts$transcript_start), 
		strand = Rle(data.transcripts$strand),
		ensembl_transcript_id = data.transcripts$ensembl_transcript_id,
		seqinfo = seqinfo)

	ranges.transcript_ends = GRanges(
		seqnames = Rle(data.transcripts$chromosome_name), 
		ranges = IRanges(start = data.transcripts$transcript_end, end = data.transcripts$transcript_end), 
		strand = Rle(data.transcripts$strand),
		ensembl_transcript_id = data.transcripts$ensembl_transcript_id,
		seqinfo = seqinfo)

	ranges.exon_starts.internal = setdiff(ranges.exon_starts, ranges.transcript_starts)
	ranges.exon_ends.internal = setdiff(ranges.exon_ends, ranges.transcript_ends)

	ranges.splice_site_starts = GRanges(
		seqnames = seqnames(ranges.exon_starts.internal),
		ranges = IRanges(start = start(ranges.exon_starts.internal) - 2, end = start(ranges.exon_starts.internal) - 1),
		strand = strand(ranges.exon_starts.internal),
		seqinfo = seqinfo)

	ranges.splice_site_ends = GRanges(
		seqnames = seqnames(ranges.exon_ends.internal),
		ranges = IRanges(start = start(ranges.exon_ends.internal) + 1, end = start(ranges.exon_ends.internal) + 2),
		strand = strand(ranges.exon_ends.internal),
		seqinfo = seqinfo)

	ranges.splice_sites = union(ranges.splice_site_starts, ranges.splice_site_ends, ignore.strand = TRUE)

	data.coding = data[!is.na(data$genomic_coding_start),]
	ranges.coding = GRanges(
		seqnames = Rle(data.coding$chromosome_name), 
		ranges = IRanges(start = data.coding$genomic_coding_start, end = data.coding$genomic_coding_end), 
		strand = Rle(data.coding$strand), 
		ensembl_gene_id = data.coding$ensembl_gene_id, 
		ensembl_transcript_id = data.coding$ensembl_transcript_id, 
		exon_rank = data.coding$rank,
		seqinfo = seqinfo)

	ranges.introns = setdiff(ranges.transcripts, ranges.exons)

	ranges.genome = GRanges(
		seqnames = Rle(seqnames(seqinfo)),
		ranges = IRanges(start = 1, end = seqlengths(seqinfo)),
		strand = "*",
		seqinfo = seqinfo)

	ranges.intergenic = setdiff(ranges.genome, ranges.genes, ignore.strand = TRUE)

	list(
		genome = ranges.genome,
		genes = ranges.genes,
		intergenic = ranges.intergenic,
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


# Define the 'validation target genome' -- for validation, we disregard all variants outside of 
# this set of ranges.  At this point, the target genome is autosomes and sex chromosomes 
# only, with mtDNA.  No patch contigs, or decoys etc.
# The reportable range will additionally be restricted to nuclear DNA only -- no mtDNA.
genome.bsgenome = getBSgenome("BSgenome.HSapiens.1000g.37d5")
genome.seqinfo = seqinfo(genome.bsgenome)
genome.seqinfo = genome.seqinfo[c(as.character(1:22), "X", "Y", "MT")]

ensembl = useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

data = getBM(
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


regions.functional = createFunctionalRegions(data, genome.seqinfo)
regions.aspirational_rr = createAspirationalReportableRangeRegions(data, genome.seqinfo, slop_delta = 10)

regions.functional.exonic.coding = reduce(regions.functional$coding)
regions.functional.exonic.5utr = reduce(regions.functional$utr5)
regions.functional.exonic.3utr = reduce(regions.functional$utr3)
regions.functional.exonic.utr = reduce(union(regions.functional$utr5, regions.functional$utr3))
regions.functional.splice = reduce(regions.functional$splice_sites)
regions.functional.intronic = reduce(regions.functional$introns)
regions.functional.genes = reduce(regions.functional$genes)
regions.functional.intergenic = reduce(regions.functional$intergenic)
regions.functional.target_genome = reduce(regions.functional$genome)

regions.aspirational_rr = reduce(intersect(regions.aspirational_rr, regions.functional.target_genome, ignore.strand = TRUE))

writeBed(regions.functional.exonic.coding, "grch37_ensembl.exonic_coding.unsorted.bed")
writeBed(regions.functional.exonic.utr, "grch37_ensembl.exonic_utr.unsorted.bed")
writeBed(regions.functional.exonic.5utr, "grch37_ensembl.exonic_5utr.unsorted.bed")
writeBed(regions.functional.exonic.3utr, "grch37_ensembl.exonic_3utr.unsorted.bed")
writeBed(regions.functional.splice, "grch37_ensembl.splice.unsorted.bed")
writeBed(regions.functional.intronic, "grch37_ensembl.intronic.unsorted.bed")
writeBed(regions.functional.genes, "grch37_ensembl.genes.unsorted.bed")
writeBed(regions.functional.intergenic, "grch37_ensembl.intergenic.unsorted.bed")
writeBed(regions.functional.target_genome, "grch37_ensembl.targetgenome.unsorted.bed")

writeBed(regions.aspirational_rr, "grch37_ensembl.cds_slop10.unsorted.bed")
