####################################################################
#  
#  KCCG WGS Validation Reporter -- Report calculations
#  
#  Usage: 
#    Rscript report_calculations.R <debugflag> <debugchrom>
#      <input_vcf> <tp> <tp_baseline> <fp> <fn> <gold_vcf> <genome> 
#      <gold_regions> <func_regions> <mask_regions>
#  
#  Positional parameters:
#    debugflag      Debug mode flag.  1 if this is a debug run,
#                   any other value otherwise.
#    debugchrom     If debugging, limit analysis to this chromosome.
#                   If not in debug mode, this value is ignored.
#    input_vcf      Path the the input VCF (can be .vcf or .vcf.gz)
#    tp, fp, fn     Paths to the overlap files output by vcfeval
#    tp_baseline    Path to the tp overlap file with baseline 
#                   annotations.
#    gold_vcf       Path to the gold-standard VCF (can be .vcf or 
#                   .vcf.gz)
#    genome         Genome label, as used by the VariantAnnotation
#                   package (eg. "hg19")
#    gold_regions   A bed or bed.gz of regions that are considered to
#                   be callable in the gold standard NA12878 data.
#    func_regions   A path prefix for the functional region BEDs. 
#    mask_regions   A path prefix for the masking region BEDs. 
#  
#  
#  Mark Pinese, 2015
#  
####################################################################

options(stringsAsFactors = FALSE, warn = 2)     # Treat warnings as errors -- this will be production code

source("report_functions.R")


####################################################################
# COMMAND LINE PARSING
####################################################################
argv = commandArgs(TRUE)
if (length(argv) != 12)
    stop("Usage: Rscript report_calculations.R <debugflag> <debugchrom> <input_vcf> <tp> <tp_baseline> <fp> <fn> <gold_vcf> <genome> <gold_regions> <func_regions> <mask_regions>")

DEBUG = argv[1] == "1"
DEBUG.chrom = argv[2]
path.input = argv[3]
path.tp = argv[4]
path.tp.baseline = argv[5]
path.fp = argv[6]
path.fn = argv[7]
path.gold = argv[8]
genome = argv[9]
path.gold_regions = argv[10]
path.function_regions_prefix = argv[11]
path.mask_regions_prefix = argv[12]


#####################################################################
# DEBUG SETTINGS -- REMOVE FOR PRODUCTION
    source("~/repos/validation-reporter/report_functions.R")
    argv = "None -- debug settings used"
    DEBUG = TRUE
    DEBUG.chrom = "1"
    path.input = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/test_data/HiSeqX_v2_TKCC/calls.vcf.gz"
    path.tp = "/directflow/ClinicalGenomicsPipeline/tmp/valrept.0aVRkXV53Q/overlap/tp.vcf.gz"
    path.fp = "/directflow/ClinicalGenomicsPipeline/tmp/valrept.0aVRkXV53Q/overlap/fp.vcf.gz"
    path.fn = "/directflow/ClinicalGenomicsPipeline/tmp/valrept.0aVRkXV53Q/overlap/fn.vcf.gz"
    path.tp.baseline = "/directflow/ClinicalGenomicsPipeline/tmp/valrept.0aVRkXV53Q/overlap/tp-baseline.vcf.gz"
    path.gold = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/gold_standard/calls.vcf.gz"
    genome = "BSgenome.HSapiens.1000g.37d5"
    path.gold_regions = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/gold_standard/valid_regions.bed.gz"
    path.function_regions_prefix = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/functional_regions/"
    path.mask_regions_prefix = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/mask_regions/"
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
    message(sprintf("  path.function_regions_prefix: %s", path.function_regions_prefix))
    message(sprintf("  path.mask_regions_prefix:     %s", path.mask_regions_prefix))
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
genome(genome.seqinfo) = genome     # To get around disagreement 
        # between the vcfs (which by default use the genome name 
        # abbreviation), and the beds (which use the full name).

# The call overlaps.  If in debug mode, subset to chromosome DEBUG.chrom
# Note the use of suppressWarnings.  The following readVcf calls result
# in warnings of the form:
#   Warning in FUN(X[[5L]], ...) :
#     duplicate ID's in header will be forced to unique rownames
#   Calls: readVcf ... .bcfHeaderAsSimpleList -> tapply -> tapply -> lapply -> FUN
# These can be safely ignored -- for now -- they are a 
# consequence of the VCF containing duplicate 
#   ##GATKCommandLine=<ID=ApplyRecalibration, ...
# lines.  TODO: Consider tweaking the front-end so that these entries
# are de-duplicated.  Then, the suppressWarnings calls can be 
# removed here, and genuine file read errors will be caught.
if (DEBUG)
{
    temp = as.data.frame(seqinfo(genome.bsgenome))
    DEBUG.region = GRanges(seqnames = DEBUG.chrom, IRanges(1, temp[DEBUG.chrom,]$seqlengths), seqinfo = genome.seqinfo)
    data.tp = suppressWarnings(readVcf(TabixFile(path.tp), genome, ScanVcfParam(which = DEBUG.region)))
    data.fp = suppressWarnings(readVcf(TabixFile(path.fp), genome, ScanVcfParam(which = DEBUG.region)))
    data.fn = suppressWarnings(readVcf(TabixFile(path.fn), genome, ScanVcfParam(which = DEBUG.region)))
    data.tp.baseline = suppressWarnings(readVcf(TabixFile(path.tp.baseline), genome, ScanVcfParam(which = DEBUG.region)))
} else {
    data.tp = suppressWarnings(readVcf(TabixFile(path.tp), genome))
    data.fp = suppressWarnings(readVcf(TabixFile(path.fp), genome))
    data.fn = suppressWarnings(readVcf(TabixFile(path.fn), genome))
    data.tp.baseline = suppressWarnings(readVcf(TabixFile(path.tp.baseline), genome))
}

# Simple data sanity check
temp.sample.tp = header(data.tp)@samples
temp.sample.fp = header(data.fp)@samples
stopifnot(temp.sample.tp == temp.sample.fp)
stopifnot(length(temp.sample.tp) == 1)

data.sampleid = header(data.tp)@samples

# The gold standard calls, subsetting under debug as before
if (DEBUG)
{
    data.gold = suppressWarnings(readVcf(TabixFile(path.gold), genome, ScanVcfParam(which = DEBUG.region)))
} else {
    data.gold = suppressWarnings(readVcf(TabixFile(path.gold), genome))
}

# The gold standard valid call regions
regions.gold = bed2GRanges(path.gold_regions, genome.seqinfo)

# The masking beds
regions.mask = readMaskRegions(path.mask_regions_prefix, genome.seqinfo)

# The 'function classes' of the genome
regions.function = readFunctionalRegions(path.function_regions_prefix, genome.seqinfo)

# If we're debugging (chr DEBUG.chrom only), subset these beds to just this area
if (DEBUG)
{
    regions.gold = intersect(regions.gold, DEBUG.region)
    regions.mask = lapply(regions.mask, function(x) intersect(x, DEBUG.region))
    regions.function = lapply(regions.function, function(x) intersect(x, DEBUG.region))
}



#####################################################################
# SUBSET CALLS
#####################################################################
# Subset all calls to the gold standard valid regions only.
data.tp = data.tp[queryHits(findOverlaps(rowData(data.tp), regions.gold, type = "within", maxgap = 0))]
data.fp = data.fp[queryHits(findOverlaps(rowData(data.fp), regions.gold, type = "within", maxgap = 0))]
data.fn = data.fn[queryHits(findOverlaps(rowData(data.fn), regions.gold, type = "within", maxgap = 0))]
data.tp.baseline = data.tp.baseline[queryHits(findOverlaps(rowData(data.tp.baseline), regions.gold, type = "within", maxgap = 0))]
data.gold = data.gold[queryHits(findOverlaps(rowData(data.gold), regions.gold, type = "within", maxgap = 0))]       # Just in case

#####################################################################
# SUBSET REGIONS
#####################################################################
# Perform this subsetting on the genome regions also
regions.mask = lapply(regions.mask, function(x) intersect(x, regions.gold))
regions.function = lapply(regions.function, function(x) intersect(x, regions.gold))


#####################################################################
# CLASSIFY CALLS
#####################################################################
# By gold standard zygosity
class.tp.zyg = simplifyZygosityClass(classifyZygosity(data.tp.baseline))        # This works as data.tp and data.tp.baseline have identical entries, in the same order
class.fn.zyg = simplifyZygosityClass(classifyZygosity(data.fn))
class.fp.zyg = data.frame("R" = rep(TRUE, nrow(data.fp)))       #
class.fp.zyg$A = FALSE                                          # FP will always be hom ref in the baseline
class.fp.zyg$RA = FALSE                                         #
class.fp.zyg$AB = FALSE                                         #


# By somy
class.tp.somy = classifySomy(data.tp.baseline)
class.fn.somy = classifySomy(data.fn)
class.fp.somy = classifySomy(data.fp)


# By sequence function (coding exonic, splice, intronic, UTR, intergenic)
# Although we have VEP calls for TP and FP, this will make some of the
# validation dependent on the pipeline classifications, which is taboo.
# To get around this problem, use fixed genome region BEDs, which have
# been independently derived using the utils/makeGenomeRegions scripts.
regions.function.min_overlap_levels = c("coding" = 1, "genic" = 1, "splice" = 1, "intronic" = 2, "utr" = 2, "intergenic" = 2)
class.tp.function = simplifyRegionClass(classifyRegion(data.tp, regions.function), regions.function.min_overlap_levels)
class.fp.function = simplifyRegionClass(classifyRegion(data.fp, regions.function), regions.function.min_overlap_levels)
class.fn.function = simplifyRegionClass(classifyRegion(data.fn, regions.function), regions.function.min_overlap_levels)


# By masking status
regions.mask.min_overlap_levels = c("ambiguous" = 1, "low_complexity" = 1, "repetitive" = 1)
class.tp.mask = simplifyRegionClass(classifyRegion(data.tp, regions.mask), regions.mask.min_overlap_levels)
class.fp.mask = simplifyRegionClass(classifyRegion(data.fp, regions.mask), regions.mask.min_overlap_levels)
class.fn.mask = simplifyRegionClass(classifyRegion(data.fn, regions.mask), regions.mask.min_overlap_levels)
class.tp.mask = cbind(class.tp.mask, "unmasked" = !apply(class.tp.mask, 1, any))
class.fp.mask = cbind(class.fp.mask, "unmasked" = !apply(class.fp.mask, 1, any))
class.fn.mask = cbind(class.fn.mask, "unmasked" = !apply(class.fn.mask, 1, any))


# By type: SNV, insertion, deletion, other
# Note a nuance here: tp and fn are based on the baseline alt allele, 
# whereas fp is based on the called alt allele.  This is only relevant
# in the fn vs fp case, as the tps will be the same between both 
# baseline and the test calls.
class.tp.muttype = classifyMutationType(data.tp.baseline)       # Based on baseline
class.fn.muttype = classifyMutationType(data.fn)                # Based on baseline
class.fp.muttype = classifyMutationType(data.fp)                # Based on calls

# The following may be an interesting diagnostic:
# table(data.frame(
#   CallType = rep(c("TP", "FN", "FP"), c(length(class.tp.muttype), length(class.fn.muttype), length(class.fp.muttype))), 
#   MutType = c(as.character(class.tp.muttype), as.character(class.fn.muttype), as.character(class.fp.muttype))))
# It indicates that TPs are predominantly SNVs.  FP is enriched for
# indels (particularly deletions), and FN for both classes of indels.


# By the 'size' of mutation (see getMutationSize for a definition of 
# size in this context).  As for muttype, measures slightly different 
# quantities in the fn and fp groups.
class.tp.mutsize = classifyMutationSize(getMutationSize(data.tp.baseline))
class.fn.mutsize = classifyMutationSize(getMutationSize(data.fn))
class.fp.mutsize = classifyMutationSize(getMutationSize(data.fp))


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

data.snv = list(vcf.tp = data.tp.snv, vcf.fp = data.fp.snv, n.fn = snv.fn.count, n.tn = snv.tn.count)
snv.perf.all.VQSLOD = vcfPerf(data.snv, function(x) info(x)$VQSLOD)
snv.perf.all.QUAL =   vcfPerf(data.snv, function(x) rowData(x)$QUAL)
snv.perf.all.GQ =     vcfPerf(data.snv, function(x) geno(x)$GQ)
snv.perf.all.DP =     vcfPerf(data.snv, function(x) info(x)$DP)
snv.perf.all.FILTER = vcfPerf(data.snv, function(x) (rowData(x)$FILTER == "PASS")*1)
snv.perf.all.VQSLOD_FILTER = vcfPerf(data.snv, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"))
snv.perf.all.QUAL_FILTER =   vcfPerf(data.snv, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"))
snv.perf.all.GQ_FILTER =     vcfPerf(data.snv, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"))
snv.perf.all.DP_FILTER =     vcfPerf(data.snv, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"))


class.snv.zyg = list(tp = class.tp.zyg[isSNV(data.tp),], fp = class.fp.zyg[isSNV(data.fp),], fn = class.fn.zyg[isSNV(data.fn),])
snv.perf.zyg.VQSLOD = vcfPerfGrouped(data.snv, function(x) info(x)$VQSLOD, class.snv.zyg)

class.snv.somy = list(tp = class.tp.somy[isSNV(data.tp),], fp = class.fp.somy[isSNV(data.fp),], fn = class.fn.somy[isSNV(data.fn),])
snv.perf.somy.VQSLOD = vcfPerfGrouped(data.snv, function(x) info(x)$VQSLOD, class.snv.somy)

class.snv.function = list(tp = class.tp.function[isSNV(data.tp),], fp = class.fp.function[isSNV(data.fp),], fn = class.fn.function[isSNV(data.fn),])
snv.perf.function.VQSLOD = vcfPerfGrouped(data.snv, function(x) info(x)$VQSLOD, class.snv.function)

class.snv.mask = list(tp = class.tp.mask[isSNV(data.tp),], fp = class.fp.mask[isSNV(data.fp),], fn = class.fn.mask[isSNV(data.fn),])
snv.perf.mask.VQSLOD = vcfPerfGrouped(data.snv, function(x) info(x)$VQSLOD, class.snv.mask)




#####################################################################
# COMBINED PERFORMANCE
#####################################################################

data.all = list(vcf.tp = data.tp, vcf.fp = data.fp, n.fn = nrow(data.fn), n.tn = 0)

class.all.zyg = list(tp = class.tp.zyg, fp = class.fp.zyg, fn = class.fn.zyg)
all.perf.zyg.VQSLOD = vcfPerfGrouped(data.all, function(x) info(x)$VQSLOD, class.all.zyg)

class.all.somy = list(tp = class.tp.somy, fp = class.fp.somy, fn = class.fn.somy)
all.perf.somy.VQSLOD = vcfPerfGrouped(data.all, function(x) info(x)$VQSLOD, class.all.somy)

class.all.function = list(tp = class.tp.function, fp = class.fp.function, fn = class.fn.function)
all.perf.function.VQSLOD = vcfPerfGrouped(data.all, function(x) info(x)$VQSLOD, class.all.function)

class.all.mask = list(tp = class.tp.mask, fp = class.fp.mask, fn = class.fn.mask)
all.perf.mask.VQSLOD = vcfPerfGrouped(data.all, function(x) info(x)$VQSLOD, class.all.mask)

class.all.muttype = list(tp = class.tp.muttype, fp = class.fp.muttype, fn = class.fn.muttype)
all.perf.muttype.VQSLOD = vcfPerfGrouped(data.all, function(x) info(x)$VQSLOD, class.all.muttype)

class.all.mutsize = list(tp = class.tp.mutsize, fp = class.fp.mutsize, fn = class.fn.mutsize)
all.perf.mutsize.VQSLOD = vcfPerfGrouped(data.all, function(x) info(x)$VQSLOD, class.all.mutsize)




#####################################################################
# SAVE RESULTS
#####################################################################
save.image("report_data.rda")
