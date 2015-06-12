####################################################################
#  
#  KCCG WGS Validation Reporter -- Report calculations
#  
#  Usage: 
#    Rscript report_calculations.R <debugflag> <debugchrom>
#      <input_vcf> <tp> <tp_baseline> <fp> <fn> <gold_vcf> <genome> 
#      <gold_regions> <func_regions> <mask_regions> <ver_branch> 
#      <ver_commit> <ver_host>
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
#    ver_*          Version strings
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
if (length(argv) != 15)
{
    stop(sprintf("Usage: Rscript report_calculations.R <debugflag> <debugchrom> <input_vcf> <tp> <tp_baseline> <fp> <fn> <gold_vcf> <genome> <gold_regions> <func_regions> <mask_regions> <ver_branch> <ver_commit> <ver_host>\nargv = %s", paste(argv, sep = " ")))
}

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
versions = list(branch = argv[13], commit = argv[14], host = argv[15])


#####################################################################
# DEBUG SETTINGS -- REMOVE FOR PRODUCTION
    # source("~/repos/validation-reporter/report_functions.R")
    # argv = "None -- debug settings used"
    # DEBUG = TRUE
    # DEBUG.chrom = "13"
    # path.input = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/test_data/HiSeqX_v2_TKCC/calls.vcf.gz"
    # path.tp = "/directflow/ClinicalGenomicsPipeline/tmp/valrept.UONTBAfcgH/overlap/tp.vcf.gz"
    # path.fp = "/directflow/ClinicalGenomicsPipeline/tmp/valrept.UONTBAfcgH/overlap/fp.vcf.gz"
    # path.fn = "/directflow/ClinicalGenomicsPipeline/tmp/valrept.UONTBAfcgH/overlap/fn.vcf.gz"
    # path.tp.baseline = "/directflow/ClinicalGenomicsPipeline/tmp/valrept.UONTBAfcgH/overlap/tp-baseline.vcf.gz"
    # path.gold = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/gold_standard/calls.vcf.gz"
    # genome = "BSgenome.HSapiens.1000g.37d5"
    # path.gold_regions = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/gold_standard/valid_regions.bed.gz"
    # path.function_regions_prefix = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/functional_regions/"
    # path.mask_regions_prefix = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/mask_regions/"
    # versions = list(branch = "feat-subsetvariants", commit = "90e89bba08e299dec789b91b2b44f340935d2e69", host = "Linux gamma02.local 2.6.32-504.8.1.el6.x86_64 #1 SMP Wed Jan 28 21:11:36 UTC 2015 x86_64 x86_64 x86_64 GNU/Linux")
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
    message(        "  versions:")
    message(sprintf("    branch:          %s", versions$branch))
    message(sprintf("    commit:          %s", versions$commit))
    message(sprintf("    host:            %s", versions$host))
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
    vcf.scan_param.called = ScanVcfParam(info = c("DP", "GQ_MEAN", "QD", "VQSLOD"), geno = c("GT", "DP", "GQ"), which = DEBUG.region)
    vcf.scan_param.uncalled = ScanVcfParam(info = c("DP", "TYPE"), geno = c("GT", "DP", "GQ"), which = DEBUG.region)
} else {
    vcf.scan_param.called = ScanVcfParam(info = c("DP", "GQ_MEAN", "QD", "VQSLOD"), geno = c("GT", "DP", "GQ"))
    vcf.scan_param.uncalled = ScanVcfParam(info = c("DP", "TYPE"), geno = c("GT", "DP", "GQ"))
}
calls = list(
    tp = suppressWarnings(readVcf(TabixFile(path.tp), genome, vcf.scan_param.called)),
    fp = suppressWarnings(readVcf(TabixFile(path.fp), genome, vcf.scan_param.called)),
    fn = suppressWarnings(readVcf(TabixFile(path.fn), genome, vcf.scan_param.uncalled)),
    tp.baseline = suppressWarnings(readVcf(TabixFile(path.tp.baseline), genome, vcf.scan_param.uncalled)))

# Simple data sanity check
temp.sample.tp = header(calls$tp)@samples
temp.sample.fp = header(calls$fp)@samples
stopifnot(temp.sample.tp == temp.sample.fp)
stopifnot(length(temp.sample.tp) == 1)

calls.sampleid = header(calls$tp)@samples


# The gold standard calls, subsetting under debug as before
calls$gold = suppressWarnings(readVcf(TabixFile(path.gold), genome, vcf.scan_param.uncalled))


# Various genomic regions, for later subsetting of performance measures
regions = list(
    gold = list(callable = bed2GRanges(path.gold_regions, genome.seqinfo)),             # Gold standard valid call regions
    mask = readMaskRegions(path.mask_regions_prefix, genome.seqinfo),                   # Masking beds
    functional = readFunctionalRegions(path.function_regions_prefix, genome.seqinfo)    # 'Function classes' of the genome
    genome = GRanges(                                                                   # The whole genome, for set ops.
        seqnames = seqnames(genome.seqinfo), 
        ranges = IRanges(1, seqlengths(genome.seqinfo)), 
        strand = "*", seqinfo = genome.seqinfo)
)

# Create a new function class, of coding +/- 10 bp
regions$functional$coding_10 = suppressWarnings(trim(reduce(
    union(union(
        regions$functional$coding, 
        trim(flank(regions$functional$coding, 10, start = TRUE)), ignore.strand = TRUE), 
        trim(flank(regions$functional$coding, 10, start = FALSE)), ignore.strand = TRUE))))

# And a new mask class, of "unmasked"
regions$mask$unmasked = setdiff(regions$genome, reduce(union(union(regions$mask$ambiguous, regions$mask$low_complexity), regions$mask$repetitive)))


# If we're debugging (chr DEBUG.chrom only), subset all regions to just this area
if (DEBUG)
{
    regions$gold = lapply(regions$gold, function(x) intersect(x, DEBUG.region))
    regions$mask = lapply(regions$mask, function(x) intersect(x, DEBUG.region))
    regions$functional = lapply(regions$functional, function(x) intersect(x, DEBUG.region))
}



#####################################################################
# CLASSIFY CALLS
#####################################################################
class = list()

# By gold standard zygosity
class$zyg = list(
    tp = simplifyZygosityClass(classifyZygosity(calls$tp.baseline)),        # This works as data.tp and data.tp.baseline have identical entries, in the same order
    fn = simplifyZygosityClass(classifyZygosity(calls$fn)),
    fp = simplifyZygosityClass(classifyZygosity(calls$fp)))


# By somy
class$somy = list(
    tp = classifySomy(calls$tp.baseline),   # Note the use of baseline here
    fn = classifySomy(calls$fn),
    fp = classifySomy(calls$fp))


# By presence / absence in the gold standard callable regions
class$goldcall = lapply(calls, function(.calls) simplifyRegionClass(classifyRegion(.calls, regions$gold), 
    c("callable" = 2)))
class$goldcall = lapply(class$goldcall, function(x) cbind(x, "not_callable" = !x[,"callable"]))

# table(class$goldcall$tp[,1])
# table(class$goldcall$fp[,1])
# table(class$goldcall$fn[,1])


# By sequence function (coding exonic, splice, intronic, UTR, intergenic)
# Although we have VEP calls for TP and FP, this will make some of the
# validation dependent on the pipeline classifications, which is taboo.
# To get around this problem, use fixed genome region BEDs, which have
# been independently derived using the utils/makeGenomeRegions scripts.
class$functional = lapply(calls, function(.calls) simplifyRegionClass(classifyRegion(.calls, regions$functional), 
    c("coding" = 1, "coding_10" = 1, "genic" = 1, "splice" = 1, "intronic" = 2, "utr" = 2, "intergenic" = 2)))


# By masking status
class$mask = lapply(calls, function(.calls) simplifyRegionClass(classifyRegion(.calls, regions$mask), c("ambiguous" = 1, "low_complexity" = 1, "repetitive" = 1)))


# By type: SNV, insertion, deletion, other
# Note a nuance here: tp and fn are based on the baseline alt allele, 
# whereas fp is based on the called alt allele.  This is only relevant
# in the fn vs fp case, as the tps will be the same between both 
# baseline and the test calls.
class$muttype = list(
    tp = classifyMutationType(calls$tp.baseline),       # Based on baseline (though *should* be the same for actual calls)
    fp = classifyMutationType(calls$fp),                # Based on calls
    fn = classifyMutationType(calls$fn))                # Based on baseline

# The following may be an interesting diagnostic:
# table(data.frame(
#   CallType = rep(c("TP", "FN", "FP"), c(length(class.muttype$tp), length(class.muttype$fn), length(class.muttype$fp))), 
#   MutType = c(as.character(class.muttype$tp), as.character(class.muttype$fn), as.character(class.muttype$fp))))
# It indicates that TPs are predominantly SNVs.  FP is enriched for
# indels (particularly deletions), and FN for both classes of indels.


# By the 'size' of mutation (see getMutationSize for a definition of 
# size in this context).  As for muttype, measures slightly different 
# quantities in the fn and fp groups.
class$mutsize = list(
    tp = classifyMutationSize(getMutationSize(calls$tp.baseline)),
    fp = classifyMutationSize(getMutationSize(calls$fp)),
    fn = classifyMutationSize(getMutationSize(calls$fn)))



#####################################################################
# COMPRESS CALL CLASSES
#####################################################################
# The class list gets pretty big.  Compress it.  As part of this,
# unfortunately, the type of the class data changes, from
#   list(list(data.frame(logical vector)))
# to
#   list(list(list(Rle)))
# TODO: Re-engineer the above classification functions to return
# their results in list(list(list(vector))) or list(list(list(Rle)))
# form, so that pre-compressed and post-compressed data have a 
# similar method of access.
class = lapply(class, function(type) lapply(type, function(perf_categ) apply(perf_categ, 2, Rle)))


# FYI for the following section: VariantAnnotation mutation classes:
# 
# • isSNV: Reference and alternate alleles are both a single
#   nucleotide long.
# 
# • isInsertion: Reference allele is a single nucleotide and the
#   alternate allele is greater (longer) than a single nucleotide
#   and the first nucleotide of the alternate allele matches the
#   reference.
# 
# • isDeletion: Alternate allele is a single nucleotide and the
#   reference allele is greater (longer) than a single nucleotide
#   and the first nucleotide of the reference allele matches the
#   alternate.
# 
# • isIndel: The variant is either a deletion or insertion as
#   determined by ‘isDeletion’ and ‘isInsertion’.
# 
# • isSubstition: Reference and alternate alleles are the same
#   length (1 or more nucleotides long).
# 
# • isTransition: Reference and alternate alleles are both a
#   single nucleotide long.  The reference-alternate pair
#   interchange is of either two-ring purines (A <-> G) or
#   one-ring pyrimidines (C <-> T).


#####################################################################
# SNV METRIC COMPARISON AND CUTOFF DETERMINATION: CODING_10 REGIONS
#####################################################################
# Subset to SNVs in gold-callable regions
# NOTE: subset.snv.regions, and subset.snv, MUST MATCH
# subset.snv.regions = intersect(regions$gold$callable, regions$functional$coding_10)
# subset.snv = sapply(names(calls), function(name) isSNV(calls[[name]]) & class$goldcall[[name]][["callable"]] & class$functional[[name]][["coding_10"]], simplify = FALSE, USE.NAMES = TRUE)
subset.snv.regions = regions$gold$callable
subset.snv = sapply(names(calls), function(name) isSNV(calls[[name]]) & class$goldcall[[name]][["callable"]], simplify = FALSE, USE.NAMES = TRUE)

calls.snv = sapply(names(calls), function(name) calls[[name]][subset.snv[[name]]], simplify = FALSE, USE.NAMES = TRUE)

# Sanity check -- are these in fact SNVs?
stopifnot(sum(width(calls.snv$tp)) == nrow(calls.snv$tp))
stopifnot(sum(width(calls.snv$fp)) == nrow(calls.snv$fp))
stopifnot(sum(width(calls.snv$fn)) == nrow(calls.snv$fn))
stopifnot(sum(width(calls.snv$gold)) == nrow(calls.snv$gold))

# We can define TNs for SNVs
calls.snv$tn = setdiff(subset.snv.regions, union(rowData(calls.snv$tp), rowData(calls.snv$fp), rowData(calls.snv$fn), ignore.strand = TRUE))

count.snv.fn = nrow(calls.snv$fn)
count.snv.tn = sum(as.numeric(width(calls.snv$tn)))

# Calculate performance for various scores, across all
# SNVs.  Later we will also subset by zygosity, and 
# sequence context, to see if certain scores perform
# better in certain contexts.  Globally, we are also
# going to perform this analysis on indels, and complex
# variants, again examining the performance of different
# scores on different variant types.
perf.snv = list()
perfdata.snv = list(vcf.tp = calls.snv$tp, vcf.fp = calls.snv$fp, n.fn = count.snv.fn, n.tn = count.snv.tn)

perf.snv$all = list(
    VQSLOD = vcfPerf(perfdata.snv, function(x) info(x)$VQSLOD),
    QUAL = vcfPerf(perfdata.snv, function(x) rowData(x)$QUAL),
    GQ = vcfPerf(perfdata.snv, function(x) geno(x)$GQ),
    DP = vcfPerf(perfdata.snv, function(x) info(x)$DP),
    FILTER = vcfPerf(perfdata.snv, function(x) (rowData(x)$FILTER == "PASS")*1),
    "VQSLOD:FILTER" = vcfPerf(perfdata.snv, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS")),
    "QUAL:FILTER" = vcfPerf(perfdata.snv, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS")),
    "GQ:FILTER" = vcfPerf(perfdata.snv, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS")),
    "DP:FILTER" = vcfPerf(perfdata.snv, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS")))

perf.snv$all2 = list(
    VQSLOD = list(
        no_filter = vcfPerf(perfdata.snv, function(x) info(x)$VQSLOD), 
        with_filter = vcfPerf(perfdata.snv, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"))),
    QUAL = list(
        no_filter = vcfPerf(perfdata.snv, function(x) rowData(x)$QUAL), 
        with_filter = vcfPerf(perfdata.snv, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"))),
    GQ = list(
        no_filter = vcfPerf(perfdata.snv, function(x) geno(x)$GQ), 
        with_filter = vcfPerf(perfdata.snv, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"))),
    DP = list(
        no_filter = vcfPerf(perfdata.snv, function(x) info(x)$DP), 
        with_filter = vcfPerf(perfdata.snv, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"))),
    FILTER = list(
        with_filter = vcfPerf(perfdata.snv, function(x) (rowData(x)$FILTER == "PASS")*1)))

# Subset by zygosity
# TNs belong in all three zygosity classes, as they are always of genotype RR.
class.snv.zyg = subsetClass(class$zyg, subset.snv, tn = list(RRvsAA = count.snv.tn, RRvsRA = count.snv.tn, RRvsAB = count.snv.tn))
perf.snv$zyg = list(
    VQSLOD = vcfPerfGrouped(perfdata.snv, function(x) info(x)$VQSLOD, class.snv.zyg),
    QUAL = vcfPerfGrouped(perfdata.snv, function(x) rowData(x)$QUAL, class.snv.zyg),
    GQ = vcfPerfGrouped(perfdata.snv, function(x) geno(x)$GQ, class.snv.zyg),
    DP = vcfPerfGrouped(perfdata.snv, function(x) info(x)$DP, class.snv.zyg),
    FILTER = vcfPerfGrouped(perfdata.snv, function(x) (rowData(x)$FILTER == "PASS")*1, class.snv.zyg),
    "VQSLOD:FILTER" = vcfPerfGrouped(perfdata.snv, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"), class.snv.zyg),
    "QUAL:FILTER" = vcfPerfGrouped(perfdata.snv, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"), class.snv.zyg),
    "GQ:FILTER" = vcfPerfGrouped(perfdata.snv, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"), class.snv.zyg),
    "DP:FILTER" = vcfPerfGrouped(perfdata.snv, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"), class.snv.zyg))

# And sequence context
# The nasty code for counting TNs just counts, for each class in 
# regions$mask, the number of TN bases (from calls.snv$tn) that
# overlap this mask class.
class.snv.mask = subsetClass(class$mask, subset.snv, tn = lapply(regions$mask, function(mask) sum(as.numeric(width(intersect(mask, calls.snv$tn))))))
perf.snv$mask = list(
    VQSLOD = vcfPerfGrouped(perfdata.snv, function(x) info(x)$VQSLOD, class.snv.mask),
    QUAL = vcfPerfGrouped(perfdata.snv, function(x) rowData(x)$QUAL, class.snv.mask),
    GQ = vcfPerfGrouped(perfdata.snv, function(x) geno(x)$GQ, class.snv.mask),
    DP = vcfPerfGrouped(perfdata.snv, function(x) info(x)$DP, class.snv.mask),
    FILTER = vcfPerfGrouped(perfdata.snv, function(x) (rowData(x)$FILTER == "PASS")*1, class.snv.mask),
    "VQSLOD:FILTER" = vcfPerfGrouped(perfdata.snv, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"), class.snv.mask),
    "QUAL:FILTER" = vcfPerfGrouped(perfdata.snv, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"), class.snv.mask),
    "GQ:FILTER" = vcfPerfGrouped(perfdata.snv, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"), class.snv.mask),
    "DP:FILTER" = vcfPerfGrouped(perfdata.snv, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"), class.snv.mask))



#####################################################################
# INDEL AND SUBSTITUTION METRIC COMPARISON AND CUTOFF DETERMINATION: CODING_10 REGIONS
#####################################################################

# Repeat the analysis performed in the SNV case, except this time subset to indels.
# This time, there is no known TN background (and no practical universe of possible
# mutations).
#subset.indelsubst = sapply(names(calls), function(name) (isIndel(calls[[name]]) | isSubstitution(calls[[name]])) & class$goldcall[[name]][["callable"]] & class$functional[[name]][["coding_10"]], simplify = FALSE, USE.NAMES = TRUE)
subset.indelsubst = sapply(names(calls), function(name) (isIndel(calls[[name]]) | isSubstitution(calls[[name]])) & class$goldcall[[name]][["callable"]], simplify = FALSE, USE.NAMES = TRUE)

calls.indelsubst = sapply(names(calls), function(name) calls[[name]][subset.indelsubst[[name]]], simplify = FALSE, USE.NAMES = TRUE)

count.indelsubst.fn = nrow(calls.indelsubst$fn)

# Score performance across all indels.  Note that n.tn = 0, as there
# is an infinite number of potential TN cases; setting n.tn to zero 
# here effectively removes these events from consideration entirely,
# and we will still have valid counts for TP, FP, and FN cases.
perf.indelsubst = list()
perfdata.indelsubst = list(vcf.tp = calls.indelsubst$tp, vcf.fp = calls.indelsubst$fp, n.fn = count.indelsubst.fn, n.tn = 0)


perf.indelsubst$all = list(
    VQSLOD = vcfPerf(perfdata.indelsubst, function(x) info(x)$VQSLOD),
    QUAL = vcfPerf(perfdata.indelsubst, function(x) rowData(x)$QUAL),
    GQ = vcfPerf(perfdata.indelsubst, function(x) geno(x)$GQ),
    DP = vcfPerf(perfdata.indelsubst, function(x) info(x)$DP),
    FILTER = vcfPerf(perfdata.indelsubst, function(x) (rowData(x)$FILTER == "PASS")*1),
    "VQSLOD:FILTER" = vcfPerf(perfdata.indelsubst, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS")),
    "QUAL:FILTER" = vcfPerf(perfdata.indelsubst, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS")),
    "GQ:FILTER" = vcfPerf(perfdata.indelsubst, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS")),
    "DP:FILTER" = vcfPerf(perfdata.indelsubst, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS")))

perf.indelsubst$all2 = list(
    VQSLOD = list(
        no_filter = vcfPerf(perfdata.indelsubst, function(x) info(x)$VQSLOD), 
        with_filter = vcfPerf(perfdata.indelsubst, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"))),
    QUAL = list(
        no_filter = vcfPerf(perfdata.indelsubst, function(x) rowData(x)$QUAL), 
        with_filter = vcfPerf(perfdata.indelsubst, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"))),
    GQ = list(
        no_filter = vcfPerf(perfdata.indelsubst, function(x) geno(x)$GQ), 
        with_filter = vcfPerf(perfdata.indelsubst, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"))),
    DP = list(
        no_filter = vcfPerf(perfdata.indelsubst, function(x) info(x)$DP), 
        with_filter = vcfPerf(perfdata.indelsubst, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"))),
    FILTER = list(
        with_filter = vcfPerf(perfdata.indelsubst, function(x) (rowData(x)$FILTER == "PASS")*1)))

# Subset by zygosity
# We can't count TNs in an indel context; just set to a null value
class.indelsubst.zyg = subsetClass(class$zyg, subset.indelsubst, tn = NULL)
perf.indelsubst$zyg = list(
    VQSLOD = vcfPerfGrouped(perfdata.indelsubst, function(x) info(x)$VQSLOD, class.indelsubst.zyg),
    QUAL = vcfPerfGrouped(perfdata.indelsubst, function(x) rowData(x)$QUAL, class.indelsubst.zyg),
    GQ = vcfPerfGrouped(perfdata.indelsubst, function(x) geno(x)$GQ, class.indelsubst.zyg),
    DP = vcfPerfGrouped(perfdata.indelsubst, function(x) info(x)$DP, class.indelsubst.zyg),
    FILTER = vcfPerfGrouped(perfdata.indelsubst, function(x) (rowData(x)$FILTER == "PASS")*1, class.indelsubst.zyg),
    "VQSLOD:FILTER" = vcfPerfGrouped(perfdata.indelsubst, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"), class.indelsubst.zyg),
    "QUAL:FILTER" = vcfPerfGrouped(perfdata.indelsubst, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"), class.indelsubst.zyg),
    "GQ:FILTER" = vcfPerfGrouped(perfdata.indelsubst, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"), class.indelsubst.zyg),
    "DP:FILTER" = vcfPerfGrouped(perfdata.indelsubst, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"), class.indelsubst.zyg))

# And sequence context
# Again, set tn = NULL
class.indelsubst.mask = subsetClass(class$mask, subset.indelsubst, tn = NULL)
perf.indelsubst$mask = list(
    VQSLOD = vcfPerfGrouped(perfdata.indelsubst, function(x) info(x)$VQSLOD, class.indelsubst.mask),
    QUAL = vcfPerfGrouped(perfdata.indelsubst, function(x) rowData(x)$QUAL, class.indelsubst.mask),
    GQ = vcfPerfGrouped(perfdata.indelsubst, function(x) geno(x)$GQ, class.indelsubst.mask),
    DP = vcfPerfGrouped(perfdata.indelsubst, function(x) info(x)$DP, class.indelsubst.mask),
    FILTER = vcfPerfGrouped(perfdata.indelsubst, function(x) (rowData(x)$FILTER == "PASS")*1, class.indelsubst.mask),
    "VQSLOD:FILTER" = vcfPerfGrouped(perfdata.indelsubst, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"), class.indelsubst.mask),
    "QUAL:FILTER" = vcfPerfGrouped(perfdata.indelsubst, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"), class.indelsubst.mask),
    "GQ:FILTER" = vcfPerfGrouped(perfdata.indelsubst, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"), class.indelsubst.mask),
    "DP:FILTER" = vcfPerfGrouped(perfdata.indelsubst, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"), class.indelsubst.mask))

# And mutation size
# Again, set tn = NULL
class.indelsubst.size = subsetClass(class$mutsize, subset.indelsubst, tn = NULL)
perf.indelsubst$mutsize = list(
    VQSLOD = vcfPerfGrouped(perfdata.indelsubst, function(x) info(x)$VQSLOD, class.indelsubst.size),
    QUAL = vcfPerfGrouped(perfdata.indelsubst, function(x) rowData(x)$QUAL, class.indelsubst.size),
    GQ = vcfPerfGrouped(perfdata.indelsubst, function(x) geno(x)$GQ, class.indelsubst.size),
    DP = vcfPerfGrouped(perfdata.indelsubst, function(x) info(x)$DP, class.indelsubst.size),
    FILTER = vcfPerfGrouped(perfdata.indelsubst, function(x) (rowData(x)$FILTER == "PASS")*1, class.indelsubst.size),
    "VQSLOD:FILTER" = vcfPerfGrouped(perfdata.indelsubst, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"), class.indelsubst.size),
    "QUAL:FILTER" = vcfPerfGrouped(perfdata.indelsubst, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"), class.indelsubst.size),
    "GQ:FILTER" = vcfPerfGrouped(perfdata.indelsubst, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"), class.indelsubst.size),
    "DP:FILTER" = vcfPerfGrouped(perfdata.indelsubst, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"), class.indelsubst.size))



#####################################################################
# COMBINED PERFORMANCE
#####################################################################

# data.all = list(vcf.tp = data.tp, vcf.fp = data.fp, n.fn = nrow(data.fn), n.tn = 0)

# class.all.zyg = list(tp = class.tp.zyg, fp = class.fp.zyg, fn = class.fn.zyg)
# all.perf.zyg.VQSLOD = vcfPerfGrouped(data.all, function(x) info(x)$VQSLOD, class.all.zyg)

# class.all.somy = list(tp = class.tp.somy, fp = class.fp.somy, fn = class.fn.somy)
# all.perf.somy.VQSLOD = vcfPerfGrouped(data.all, function(x) info(x)$VQSLOD, class.all.somy)

# class.all.function = list(tp = class.tp.function, fp = class.fp.function, fn = class.fn.function)
# all.perf.function.VQSLOD = vcfPerfGrouped(data.all, function(x) info(x)$VQSLOD, class.all.function)

# class.all.mask = list(tp = class.tp.mask, fp = class.fp.mask, fn = class.fn.mask)
# all.perf.mask.VQSLOD = vcfPerfGrouped(data.all, function(x) info(x)$VQSLOD, class.all.mask)

# class.all.muttype = list(tp = class.tp.muttype, fp = class.fp.muttype, fn = class.fn.muttype)
# all.perf.muttype.VQSLOD = vcfPerfGrouped(data.all, function(x) info(x)$VQSLOD, class.all.muttype)

# class.all.mutsize = list(tp = class.tp.mutsize, fp = class.fp.mutsize, fn = class.fn.mutsize)
# all.perf.mutsize.VQSLOD = vcfPerfGrouped(data.all, function(x) info(x)$VQSLOD, class.all.mutsize)





#####################################################################
# SAVE RESULTS
#####################################################################
# For debugging: object sizes in GB
sort(sapply(ls(), function(id) object.size(get(id)))) / 1024^3

save.image("report_data.rda")
