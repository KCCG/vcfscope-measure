####################################################################
#  
#  KCCG WGS Validation Reporter -- Report calculations
#  
#  Usage: 
#    Rscript report_calculations.R [-d] <infile> <tp> <fp> <fn> 
#      <genome> <gold_regions> <call_regions>
#  
#  Flags:
#    -d  			Turn on debug mode
#  
#  Positional parameters:
#    infile			Path to the input vcf.gz
#    tp, fp, fn		Paths to the overlap files output by vcfeval
#    genome			Genome label, as used by the VariantAnnotation
#  					package (eg. "hg19")
#    gold_regions	A bed or bed.gz of regions that are considered to
#    				be callable in the gold standard NA12878 data.
#    call_regions 	A bed or bed.gz of regions that are considered to
#    				be callable on the platform to be tested.
#  
#  
#  Mark Pinese
#  26 May 2015	MP 	Split off from report.Rnw
#  
####################################################################

options(stringsAsFactors = FALSE, warn = 1)

source("report_functions.R")


####################################################################
# COMMAND LINE PARSING
####################################################################
argv = commandArgs(TRUE)
if (length(argv) != 8)
	stop("Usage: Rscript report_calculations.R <debugflag> <debugchrom> <wd> <tp> <fp> <fn> <genome> <gold_regions> <call_regions>")

DEBUG = argv[1] == "1"
DEBUG.chrom = argv[2]
path.input = argv[3]
path.tp = argv[4]
path.fp = argv[5]
path.fn = argv[6]
genome = argv[7]
path.gold_regions = argv[8]
path.call_regions = argv[9]



#####################################################################
# DEBUG SETTINGS -- REMOVE FOR PRODUCTION
	# argv = "None -- debug settings used"
	# DEBUG = TRUE
	# path.input = "../../data/test_data/HiSeqX_v2_TKCC/calls.vcf.gz"
	# path.tp = "../overlap/tp.vcf.gz"
	# path.fp = "../overlap/fp.vcf.gz"
	# path.fn = "../overlap/fn.vcf.gz"
	# genome = "BSgenome.HSapiens.1000g.37d5"
	# path.gold_regions = "../../data/gold_standard/valid_regions.bed.gz"
	# path.call_regions = "../../data/kccg/not_hardmasked.bed.gz"
# END DEBUG SETTINGS
#####################################################################



if (DEBUG)
{
	message(sprintf("Command line: %s", paste(argv, collapse = " ")))
	message(sprintf("  DEBUG:             %s", DEBUG))
	message(sprintf("  DEBUG.chrom:       %s", DEBUG.chrom))
	message(sprintf("  path.tp:           %s", path.tp))
	message(sprintf("  path.fp:           %s", path.fp))
	message(sprintf("  path.fn:           %s", path.fn))
	message(sprintf("  genome:            %s", genome))
	message(sprintf("  path.gold_regions: %s", path.gold_regions))
	message(sprintf("  path.call_regions: %s", path.call_regions))
}


#####################################################################
# LIBRARIES
#####################################################################
library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome)

library(bit64)

library(ROCR)


#####################################################################
# LOAD DATA
#####################################################################
# The reference genome
genome.bsgenome = getBSgenome(genome)
genome.seqinfo = seqinfo(genome.bsgenome)
genome(genome.seqinfo) = genome 	# To get around disagreement 
		# between the vcfs (which by default use the genome name 
		# abbreviation), and the beds (which use the full name).

# The call overlaps.  If in debug mode, subset to chromosome DEBUG.chrom
if (DEBUG)
{
	temp = as.data.frame(seqinfo(genome.bsgenome))
	debug.region = GRanges(seqnames = DEBUG.chrom, IRanges(1, temp[DEBUG.chrom,]$seqlengths), seqinfo = genome.seqinfo)
	data.tp = readVcf(TabixFile(path.tp), genome, ScanVcfParam(which = debug.region))
	data.fp = readVcf(TabixFile(path.fp), genome, ScanVcfParam(which = debug.region))
	data.fn = readVcf(TabixFile(path.fn), genome, ScanVcfParam(which = debug.region))
} else {
	data.tp = readVcf(TabixFile(path.tp), genome)
	data.fp = readVcf(TabixFile(path.fp), genome)
	data.fn = readVcf(TabixFile(path.fn), genome)
}

# Simple data sanity check
temp.sample.tp = header(data.tp)@samples
temp.sample.fp = header(data.fp)@samples
stopifnot(temp.sample.tp == temp.sample.fp)
stopifnot(length(temp.sample.tp) == 1)

data.sampleid = header(data.tp)@samples

# The region beds
regions.gold = bed2GRanges(path.gold_regions, genome.seqinfo)
regions.call = bed2GRanges(path.call_regions, genome.seqinfo)

# If we're debugging (chr DEBUG.chrom only), subset these beds to just this area
if (DEBUG)
{
	regions.gold = intersect(regions.gold, debug.region)
	regions.call = intersect(regions.call, debug.region)
}



#####################################################################
# SUBSET CALLS
#####################################################################
# Subset all calls to the gold standard valid regions only.
data.tp = data.tp[queryHits(findOverlaps(rowData(data.tp), regions.gold, type = "within", maxgap = 0))]
data.fp = data.fp[queryHits(findOverlaps(rowData(data.fp), regions.gold, type = "within", maxgap = 0))]
data.fn = data.fn[queryHits(findOverlaps(rowData(data.fn), regions.gold, type = "within", maxgap = 0))]



#####################################################################
# SNV PERFORMANCE
#####################################################################
data.tp.snv = data.tp[isSNV(data.tp)]
data.fp.snv = data.fp[isSNV(data.fp)]
data.fn.snv = data.fn[isSNV(data.fn)]

stopifnot(sum(width(data.tp.snv)) == nrow(data.tp.snv))
stopifnot(sum(width(data.fp.snv)) == nrow(data.fp.snv))
stopifnot(sum(width(data.fn.snv)) == nrow(data.fn.snv))

snv.tn = setdiff(regions.gold, union(rowData(data.tp.snv), rowData(data.fn.snv), rowData(data.fp.snv), ignore.strand = TRUE))

snv.tp.count = nrow(data.tp.snv)
snv.fp.count = nrow(data.fp.snv)
snv.fn.count = nrow(data.fn.snv)
snv.tn.count = sum(as.integer64(width(snv.tn)))

snv.truth = makeTruthVector(snv.tp.count, snv.fp.count, snv.fn.count, snv.tn.count)
snv.score.VQSLOD = makeScoreVector(data.tp.snv, data.fp.snv, snv.fn.count, snv.tn.count, function(x) info(x)$VQSLOD)
snv.score.QUAL =   makeScoreVector(data.tp.snv, data.fp.snv, snv.fn.count, snv.tn.count, function(x) rowData(x)$QUAL)
snv.score.GQ =     makeScoreVector(data.tp.snv, data.fp.snv, snv.fn.count, snv.tn.count, function(x) geno(x)$GQ)
snv.score.DP =     makeScoreVector(data.tp.snv, data.fp.snv, snv.fn.count, snv.tn.count, function(x) info(x)$DP)
snv.score.FILTER = makeScoreVector(data.tp.snv, data.fp.snv, snv.fn.count, snv.tn.count, function(x) (rowData(x)$FILTER == "PASS")*1)

snv.pred.VQSLOD = prediction(snv.score.VQSLOD, snv.truth)
snv.pred.QUAL =   prediction(snv.score.QUAL,   snv.truth)
snv.pred.GQ =     prediction(snv.score.GQ,     snv.truth)
snv.pred.DP =     prediction(snv.score.DP,     snv.truth)
snv.pred.FILTER = prediction(snv.score.FILTER, snv.truth)


#####################################################################
# SAVE RESULTS
#####################################################################
save.image("report_data.rda")
