####################################################################
#  
#  KCCG WGS Validation Reporter -- Report calculations
#  
#  Usage: 
#    Rscript report_calculations.R <debugflag> <debugchrom>
#      <input_vcf> <tp> <fp> <fn> <gold_vcf> <genome> 
#      <gold_regions> <call_regions>
#  
#  Positional parameters:
#    debugflag		Debug mode flag.  1 if this is a debug run,
#    				any other value otherwise.
#    debugchrom		If debugging, limit analysis to this chromosome.
#    				If not in debug mode, this value is ignored.
#    input_vcf		Path the the input VCF (can be .vcf or .vcf.gz)
#    tp, fp, fn		Paths to the overlap files output by vcfeval
#    gold_vcf		Path to the gold-standard VCF (can be .vcf or 
#                   .vcf.gz)
#    genome			Genome label, as used by the VariantAnnotation
#    				package (eg. "hg19")
#    gold_regions	A bed or bed.gz of regions that are considered to
#    				be callable in the gold standard NA12878 data.
#    call_regions 	A bed or bed.gz of regions that are considered to
#    				be callable on the platform to be tested.
#  
#  
#  Mark Pinese, 2015
#  
####################################################################

options(stringsAsFactors = FALSE, warn = 1)

source("report_functions.R")


####################################################################
# COMMAND LINE PARSING
####################################################################
argv = commandArgs(TRUE)
if (length(argv) != 11)
	stop("Usage: Rscript report_calculations.R <debugflag> <debugchrom> <input_vcf> <tp> <fp> <fn> <tp_baseline> <gold_vcf> <genome> <gold_regions> <call_regions>")

DEBUG = argv[1] == "1"
DEBUG.chrom = argv[2]
path.input = argv[3]
path.tp = argv[4]
path.fp = argv[5]
path.fn = argv[6]
path.tp.baseline = argv[7]
path.gold = argv[8]
genome = argv[9]
path.gold_regions = argv[10]
path.call_regions = argv[11]



#####################################################################
# DEBUG SETTINGS -- REMOVE FOR PRODUCTION
	# argv = "None -- debug settings used"
	# DEBUG = TRUE
	# DEBUG.chrom = "X"
	# path.input = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/test_data/HiSeqX_v2_TKCC/calls.vcf.gz"
	# path.tp = "/directflow/ClinicalGenomicsPipeline/tmp/valrept.bwifPHtQ1x/overlap/tp.vcf.gz"
	# path.fp = "/directflow/ClinicalGenomicsPipeline/tmp/valrept.bwifPHtQ1x/overlap/fp.vcf.gz"
	# path.fn = "/directflow/ClinicalGenomicsPipeline/tmp/valrept.bwifPHtQ1x/overlap/fn.vcf.gz"
	# path.tp.baseline = "/directflow/ClinicalGenomicsPipeline/tmp/valrept.bwifPHtQ1x/overlap/tp-baseline.vcf.gz"
	# path.gold = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/gold_standard/calls.vcf.gz"
	# genome = "BSgenome.HSapiens.1000g.37d5"
	# path.gold_regions = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/gold_standard/valid_regions.bed.gz"
	# path.call_regions = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/kccg/not_hardmasked.bed.gz"
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
	message(sprintf("  path.gold:         %s", path.gold))
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
	data.tp.baseline = readVcf(TabixFile(path.tp.baseline), genome, ScanVcfParam(which = debug.region))
} else {
	data.tp = readVcf(TabixFile(path.tp), genome)
	data.fp = readVcf(TabixFile(path.fp), genome)
	data.fn = readVcf(TabixFile(path.fn), genome)
	data.tp.baseline = readVcf(TabixFile(path.tp.baseline), genome)
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
data.tp.baseline = data.tp.baseline[queryHits(findOverlaps(rowData(data.tp.baseline), regions.gold, type = "within", maxgap = 0))]


#####################################################################
# CLASSIFY CALLS
#####################################################################
# By gold standard zygosity
class.tp.zyg = classifyZygosity(data.tp.baseline)		# This works as data.tp and data.tp.baseline have identical entries, in the same order
class.fn.zyg = classifyZygosity(data.fn)
# class.fp.zyg and class.tn.zyg are not useful: there is no such 
# thing as gold standard zygosity for FP samples, and TN samples are 
# always homozygous reference.

# By region (coding exonic, splice, intronic, NC exonic, intergenic)
# (Let splice override NC exonic)
# TP and FP we have for free from the VEP calls, but TN and FN will
# be much harder.  And if we need to recalculate TN and FN, we may
# as well do everything, for consistency.
# It's tempting to put the gold standard baseline VCF through VEP, but
# the issue of VEP version matching could be serious.  I'd rather
# not recode VEP, so perhaps just a very simple breakdown of location,
# based on the UCSC tables or similar, will do.
#class.tp.region = classifyRegion()
#class.fn.region = 
# SO Splice region: Within 1-3 bases of exon, or 3-8 bases of intron

# By type: SNV, insertion, deletion, other
# Note a nuance here: tp and fn are based on the baseline alt allele, 
# whereas fp is based on the called alt allele.  This is only relevant
# in the fn vs fp case, as the tps will be the same between both 
# baseline and the test calls.
class.tp.muttype = classifyMutationType(data.tp.baseline)		# Based on baseline
class.fn.muttype = classifyMutationType(data.fn)				# Based on baseline
class.fp.muttype = classifyMutationType(data.fp)				# Based on calls

# The following may be an interesting diagnostic:
# table(data.frame(
# 	CallType = rep(c("TP", "FN", "FP"), c(length(class.tp.muttype), length(class.fn.muttype), length(class.fp.muttype))), 
# 	MutType = c(as.character(class.tp.muttype), as.character(class.fn.muttype), as.character(class.fp.muttype))))
# It indicates that TPs are predominantly SNVs.  FP is enriched for
# indels (particularly deletions), and FN for both classes of indels.

# The 'size' of mutation (see getMutationSize for a definition of size
# in this context).  As for muttype, measures slightly different 
# quantities in the fn and fp groups.
class.tp.mutsize = getMutationSize(data.tp.baseline)
class.fn.mutsize = getMutationSize(data.fn)
class.fp.mutsize = getMutationSize(data.fp)


#####################################################################
# SNV PERFORMANCE
#####################################################################
data.tp.snv = data.tp[isSNV(data.tp)]
data.fp.snv = data.fp[isSNV(data.fp)]
data.fn.snv = data.fn[isSNV(data.fn)]
data.gold.snv = data.gold[isSNV(data.gold)]

stopifnot(sum(width(data.tp.snv)) == nrow(data.tp.snv))
stopifnot(sum(width(data.fp.snv)) == nrow(data.fp.snv))
stopifnot(sum(width(data.fn.snv)) == nrow(data.fn.snv))
stopifnot(sum(width(data.gold.snv)) == nrow(data.gold.snv))

snv.tn = setdiff(regions.gold, union(rowData(data.tp.snv), rowData(data.fn.snv), rowData(data.fp.snv), ignore.strand = TRUE))

snv.tp.count = nrow(data.tp.snv)
snv.fp.count = nrow(data.fp.snv)
snv.fn.count = nrow(data.fn.snv)
snv.tn.count = sum(as.numeric(width(snv.tn)))

snv.perf.VQSLOD = vcfPerf(data.tp.snv, data.fp.snv, snv.fn.count, snv.tn.count, function(x) info(x)$VQSLOD)
snv.perf.QUAL =   vcfPerf(data.tp.snv, data.fp.snv, snv.fn.count, snv.tn.count, function(x) rowData(x)$QUAL)
snv.perf.GQ =     vcfPerf(data.tp.snv, data.fp.snv, snv.fn.count, snv.tn.count, function(x) geno(x)$GQ)
snv.perf.DP =     vcfPerf(data.tp.snv, data.fp.snv, snv.fn.count, snv.tn.count, function(x) info(x)$DP)
snv.perf.FILTER = vcfPerf(data.tp.snv, data.fp.snv, snv.fn.count, snv.tn.count, function(x) (rowData(x)$FILTER == "PASS")*1)


#####################################################################
# SAVE RESULTS
#####################################################################
save.image("report_data.rda")
