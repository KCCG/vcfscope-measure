####################################################################
#  
#  KCCG WGS Validation Reporter -- Report calculations
#  
#  Usage: 
#    Rscript report_calculations.R <debugflag> <debugchrom>
#      <extendedflag> <input_vcf> <tp> <fp> <fn> <gold_vcf> 
#      <genome> <gold_regions> <func_regions> <mask_regions> 
#      <ver_branch> <ver_commit> <ver_host>
#  
#  Positional parameters:
#    debugflag      Debug mode flag.  1 if this is a debug run,
#                   any other value otherwise.
#    debugchrom     If debugging, limit analysis to this chromosome.
#                   If not in debug mode, this value is ignored.
#    extendedflag   Generate an extended report, with additional
#                   plots for threshold and score diagnostics.
#    input_vcf      Path the the input VCF (can be .vcf or .vcf.gz)
#    tp, fp, fn     Paths to the overlap files output by vcfeval
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
if (length(argv) != 16)
{
    stop(sprintf("Usage: Rscript report_calculations.R <debugflag> <debugchrom> <extendedflag>\n    <input_vcf> <tp> <fp> <fn> <gold_vcf> <genome> <gold_regions> <func_regions>\n     <mask_regions> <ver_script> <ver_branch> <ver_commit> <ver_host>\n\nargv = %s", paste(argv, sep = " ")))
}

DEBUG = argv[1] == "1"
DEBUG.chrom = argv[2]
extendedflag = argv[3] == "1"
path.input = argv[4]
path.tp = argv[5]
path.fp = argv[6]
path.fn = argv[7]
path.gold = argv[8]
genome = argv[9]
path.gold_regions = argv[10]
path.function_regions_prefix = argv[11]
path.mask_regions_prefix = argv[12]
versions = list(script = argv[13], branch = argv[14], commit = argv[15], host = argv[16])


if (DEBUG)
{
    message(sprintf("Command line: %s", paste(argv, collapse = " ")))
    message(sprintf("  DEBUG:             %s", DEBUG))
    message(sprintf("  DEBUG.chrom:       %s", DEBUG.chrom))
    message(sprintf("  extendedflag:      %s", extendedflag))
    message(sprintf("  path.tp:           %s", path.tp))
    message(sprintf("  path.fp:           %s", path.fp))
    message(sprintf("  path.fn:           %s", path.fn))
    message(sprintf("  path.gold:         %s", path.gold))
    message(sprintf("  genome:            %s", genome))
    message(sprintf("  path.gold_regions: %s", path.gold_regions))
    message(sprintf("  path.function_regions_prefix: %s", path.function_regions_prefix))
    message(sprintf("  path.mask_regions_prefix:     %s", path.mask_regions_prefix))
    message(        "  versions:")
    message(sprintf("    script:          %s", versions$script))
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
# removed here, and genuine file read warnings will be caught.
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
    fn = suppressWarnings(readVcf(TabixFile(path.fn), genome, vcf.scan_param.uncalled)))


# Simple data sanity check: do the vcfs 
# have the same sample name?
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
    functional = readFunctionalRegions(path.function_regions_prefix, genome.seqinfo),   # 'Function classes' of the genome
    genome = list(genome = GRanges(                                                     # The whole genome, for set ops.
        seqnames = seqnames(genome.seqinfo), 
        ranges = IRanges(1, seqlengths(genome.seqinfo)), 
        strand = "*", seqinfo = genome.seqinfo))
)

# Create a new function class, of coding +/- 10 bp
regions$functional$coding_10 = suppressWarnings(trim(reduce(
    union(union(
        regions$functional$coding, 
        trim(flank(regions$functional$coding, 10, start = TRUE)), ignore.strand = TRUE), 
        trim(flank(regions$functional$coding, 10, start = FALSE)), ignore.strand = TRUE))))

# And new mask classes, of "unmasked" and "not gold standard callable"
regions$mask$unmasked = setdiff(regions$genome$genome, reduce(union(union(regions$mask$ambiguous, regions$mask$low_complexity), regions$mask$repetitive)))
regions$gold$notcallable = setdiff(regions$genome$genome, regions$gold$callable)

# If we're debugging (chr DEBUG.chrom only), subset all regions to just this area
if (DEBUG)
{
    regions = lapply(regions, function(region_class) lapply(region_class, function(region) intersect(region, DEBUG.region)))
}


#####################################################################
# CLASSIFY CALLS
#####################################################################

class = list(
    zyg = classifyZygosity(calls),          # By zygosity
    somy = classifySomy(calls),             # By somy
    muttype = classifyMutationType(calls),  # By mutation type: SNV, InsDelSubst, other
    mutsize = classifyMutationSize(calls),  # By mutation 'size' (see getMutationSizeVcf for the definition of size in this context)

    # By presence / absence in the gold standard callable regions
    goldcall = classifyRegionOverlap(calls, regions$gold, c("callable" = "All", "notcallable" = "Any")),
    
    # By sequence function (coding exonic, splice, intronic, UTR, intergenic).
    # Use region BEDs that have been independently derived using the 
    # utils/makeGenomeRegions scripts.
    functional = classifyRegionOverlap(calls, regions$functional, c("coding" = "Any", "coding_10" = "Any", "genic" = "Any", "splice" = "Any", "intronic" = "All", "utr" = "All", "intergenic" = "All")),

    # By masking status
    mask = classifyRegionOverlap(calls, regions$mask, c("ambiguous" = "Any", "low_complexity" = "Any", "repetitive" = "Any", "unmasked" = "All"))
)


# There is a bit of subtlety in this.  Some classes (zygosity, muttype,
# and mutsize), are based on the properties of the mutation, as given
# in the gold standard.  However, in the case of FP or FN calls, 
# these properties may be wrong, as they are based on a known 
# incorrect call.  To resolve this, we manually go back and set these
# calls to the correct values from the gold standard.  For the majority 
# of classes -- anything based on regions or reference sequence 
# properties -- this correction is not necessary, as these classes are 
# based on data that are not affected by an incorrect call.

# FP is always genotype RR, in truth, and so belongs in all three classes:
class$zyg$RRvsRA$fp = class$zyg$RRvsRA$fp | TRUE
class$zyg$RRvsAA$fp = class$zyg$RRvsAA$fp | TRUE
class$zyg$RRvsAB$fp = class$zyg$RRvsAB$fp | TRUE

# Mutation type fixes do not need to be performed, provided we remain
# aware of what we're doing:
# As for zygosity, TPs, TNs, and FNs need no correction.  FPs we leave 
# alone, so that they reflect the mutation type of the *called* mutation, 
# even though it was incorrect.  This gives an indication on the tendency 
# of FP calls to be one type of mutation over another.

# Likewise, we leave the mutation size classes alone.  These decisions
# become important when interpreting the final ROCs.


# Basic consistency check: ensure all class vectors match the length 
# of the corresponding call vcf.
for (temp.i in class)
{
    for (temp.j in temp.i)
    {
        for (temp.k in names(temp.j))
            stopifnot(length(temp.j[[temp.k]]) == nrow(calls[[temp.k]]))
    }
}


#####################################################################
# SNV METRIC COMPARISON AND CUTOFF DETERMINATION: GOLD CALLABLE REGIONS
#####################################################################
# Subset to SNVs in gold-callable regions
# NOTE: subset.snv.regions, and subset.snv, MUST MATCH
subset.snv.regions = regions$gold$callable
subset.snv = sapply(names(calls), function(name) class$muttype$SNV[[name]] & class$goldcall$callable[[name]], simplify = FALSE, USE.NAMES = TRUE)

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
# better in certain contexts.  Later, we will also
# perform this analysis on indels, and complex variants, 
# again examining the performance of different scores 
# on different variant types.
perf.snv = list()
perfdata.snv = list(vcf.tp = calls.snv$tp, vcf.fp = calls.snv$fp, n.fn = count.snv.fn, n.tn = count.snv.tn)

# Subset by zygosity
# TNs belong in all three zygosity classes, as they are always of genotype RR.
class.snv.zyg = subsetClass(class$zyg, subset.snv, tn = list(RRvsAA = count.snv.tn, RRvsRA = count.snv.tn, RRvsAB = count.snv.tn))
# Basic checks for subsetClass
stopifnot(all(class$zyg$RRvsRA$tp[subset.snv$tp] == class.snv.zyg$RRvsRA$tp))
stopifnot(all(class$zyg$RRvsRA$fp[subset.snv$fp] == class.snv.zyg$RRvsRA$fp))
stopifnot(all(class$zyg$RRvsAA$tp[subset.snv$tp] == class.snv.zyg$RRvsAA$tp))

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

if (extendedflag) {
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
}



#####################################################################
# INDEL AND SUBSTITUTION METRIC COMPARISON AND CUTOFF DETERMINATION: GOLD CALLABLE REGIONS
#####################################################################

# Repeat the analysis performed in the SNV case, except this time subset to indels.
# This time, there is no known TN background (and no practical universe of possible
# mutations).
subset.indelsubst = sapply(names(calls), function(name) class$muttype$InsDelSubst[[name]] & class$goldcall$callable[[name]], simplify = FALSE, USE.NAMES = TRUE)

calls.indelsubst = sapply(names(calls), function(name) calls[[name]][subset.indelsubst[[name]]], simplify = FALSE, USE.NAMES = TRUE)

count.indelsubst.fn = nrow(calls.indelsubst$fn)

# Score performance across all indels.  Note that n.tn = 0, as there
# is an infinite number of potential TN cases; setting n.tn to zero 
# here effectively removes these events from consideration entirely,
# and we will still have valid counts for TP, FP, and FN cases.
perf.indelsubst = list()
perfdata.indelsubst = list(vcf.tp = calls.indelsubst$tp, vcf.fp = calls.indelsubst$fp, n.fn = count.indelsubst.fn, n.tn = 0)

# Subset by zygosity
# We can't count TNs in an indel context; just set to a null value
class.indelsubst.zyg = subsetClass(class$zyg, subset.indelsubst, tn = NULL)
# Basic checks for subsetClass
stopifnot(all(class$zyg$RRvsRA$tp[subset.indelsubst$tp] == class.indelsubst.zyg$RRvsRA$tp))
stopifnot(all(class$zyg$RRvsRA$fp[subset.indelsubst$fp] == class.indelsubst.zyg$RRvsRA$fp))
stopifnot(all(class$zyg$RRvsAA$tp[subset.indelsubst$tp] == class.indelsubst.zyg$RRvsAA$tp))

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

if (extendedflag) {
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
}


#####################################################################
# ALL VARIANT METRIC COMPARISON AND CUTOFF DETERMINATION: GOLD CALLABLE REGIONS
#####################################################################
if (extendedflag) {
# Repeat the analysis performed in the SNV case, except this time 
# consider all variants.
# As for indels, there is no known TN background (and no practical universe of possible
# mutations).
subset.all = sapply(names(calls), function(name) class$goldcall$callable[[name]], simplify = FALSE, USE.NAMES = TRUE)

calls.all = sapply(names(calls), function(name) calls[[name]][subset.all[[name]]], simplify = FALSE, USE.NAMES = TRUE)

count.all.fn = nrow(calls.all$fn)

# Score performance.  Note that n.tn = 0, as there
# is an infinite number of potential TN cases; setting n.tn to zero 
# here effectively removes these events from consideration entirely,
# and we will still have valid counts for TP, FP, and FN cases.
perf.all = list()
perfdata.all = list(vcf.tp = calls.all$tp, vcf.fp = calls.all$fp, n.fn = count.all.fn, n.tn = 0)


perf.all$all = list(
    VQSLOD = vcfPerf(perfdata.all, function(x) info(x)$VQSLOD),
    QUAL = vcfPerf(perfdata.all, function(x) rowData(x)$QUAL),
    GQ = vcfPerf(perfdata.all, function(x) geno(x)$GQ),
    DP = vcfPerf(perfdata.all, function(x) info(x)$DP),
    FILTER = vcfPerf(perfdata.all, function(x) (rowData(x)$FILTER == "PASS")*1),
    "VQSLOD:FILTER" = vcfPerf(perfdata.all, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS")),
    "QUAL:FILTER" = vcfPerf(perfdata.all, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS")),
    "GQ:FILTER" = vcfPerf(perfdata.all, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS")),
    "DP:FILTER" = vcfPerf(perfdata.all, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS")))

perf.all$all2 = list(
    VQSLOD = list(
        no_filter = vcfPerf(perfdata.all, function(x) info(x)$VQSLOD), 
        with_filter = vcfPerf(perfdata.all, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"))),
    QUAL = list(
        no_filter = vcfPerf(perfdata.all, function(x) rowData(x)$QUAL), 
        with_filter = vcfPerf(perfdata.all, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"))),
    GQ = list(
        no_filter = vcfPerf(perfdata.all, function(x) geno(x)$GQ), 
        with_filter = vcfPerf(perfdata.all, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"))),
    DP = list(
        no_filter = vcfPerf(perfdata.all, function(x) info(x)$DP), 
        with_filter = vcfPerf(perfdata.all, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"))),
    FILTER = list(
        with_filter = vcfPerf(perfdata.all, function(x) (rowData(x)$FILTER == "PASS")*1)))

# Subset by zygosity
# We can't count TNs in an all variant context; just set to a null value
class.all.zyg = subsetClass(class$zyg, subset.all, tn = NULL)
# Basic checks for subsetClass
stopifnot(all(class$zyg$RRvsRA$tp[subset.all$tp] == class.all.zyg$RRvsRA$tp))
stopifnot(all(class$zyg$RRvsRA$fp[subset.all$fp] == class.all.zyg$RRvsRA$fp))
stopifnot(all(class$zyg$RRvsAA$tp[subset.all$tp] == class.all.zyg$RRvsAA$tp))

perf.all$zyg = list(
    VQSLOD = vcfPerfGrouped(perfdata.all, function(x) info(x)$VQSLOD, class.all.zyg),
    QUAL = vcfPerfGrouped(perfdata.all, function(x) rowData(x)$QUAL, class.all.zyg),
    GQ = vcfPerfGrouped(perfdata.all, function(x) geno(x)$GQ, class.all.zyg),
    DP = vcfPerfGrouped(perfdata.all, function(x) info(x)$DP, class.all.zyg),
    FILTER = vcfPerfGrouped(perfdata.all, function(x) (rowData(x)$FILTER == "PASS")*1, class.all.zyg),
    "VQSLOD:FILTER" = vcfPerfGrouped(perfdata.all, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"), class.all.zyg),
    "QUAL:FILTER" = vcfPerfGrouped(perfdata.all, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"), class.all.zyg),
    "GQ:FILTER" = vcfPerfGrouped(perfdata.all, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"), class.all.zyg),
    "DP:FILTER" = vcfPerfGrouped(perfdata.all, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"), class.all.zyg))

# And sequence context
# Again, set tn = NULL
class.all.mask = subsetClass(class$mask, subset.all, tn = NULL)
perf.all$mask = list(
    VQSLOD = vcfPerfGrouped(perfdata.all, function(x) info(x)$VQSLOD, class.all.mask),
    QUAL = vcfPerfGrouped(perfdata.all, function(x) rowData(x)$QUAL, class.all.mask),
    GQ = vcfPerfGrouped(perfdata.all, function(x) geno(x)$GQ, class.all.mask),
    DP = vcfPerfGrouped(perfdata.all, function(x) info(x)$DP, class.all.mask),
    FILTER = vcfPerfGrouped(perfdata.all, function(x) (rowData(x)$FILTER == "PASS")*1, class.all.mask),
    "VQSLOD:FILTER" = vcfPerfGrouped(perfdata.all, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"), class.all.mask),
    "QUAL:FILTER" = vcfPerfGrouped(perfdata.all, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"), class.all.mask),
    "GQ:FILTER" = vcfPerfGrouped(perfdata.all, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"), class.all.mask),
    "DP:FILTER" = vcfPerfGrouped(perfdata.all, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"), class.all.mask))

# And mutation type
# Again, set tn = NULL
class.all.type = subsetClass(class$muttype, subset.all, tn = NULL)
perf.all$muttype = list(
    VQSLOD = vcfPerfGrouped(perfdata.all, function(x) info(x)$VQSLOD, class.all.type),
    QUAL = vcfPerfGrouped(perfdata.all, function(x) rowData(x)$QUAL, class.all.type),
    GQ = vcfPerfGrouped(perfdata.all, function(x) geno(x)$GQ, class.all.type),
    DP = vcfPerfGrouped(perfdata.all, function(x) info(x)$DP, class.all.type),
    FILTER = vcfPerfGrouped(perfdata.all, function(x) (rowData(x)$FILTER == "PASS")*1, class.all.type),
    "VQSLOD:FILTER" = vcfPerfGrouped(perfdata.all, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"), class.all.type),
    "QUAL:FILTER" = vcfPerfGrouped(perfdata.all, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"), class.all.type),
    "GQ:FILTER" = vcfPerfGrouped(perfdata.all, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"), class.all.type),
    "DP:FILTER" = vcfPerfGrouped(perfdata.all, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"), class.all.type))

# And mutation size
# Again, set tn = NULL
class.all.size = subsetClass(class$mutsize, subset.all, tn = NULL)
perf.all$mutsize = list(
    VQSLOD = vcfPerfGrouped(perfdata.all, function(x) info(x)$VQSLOD, class.all.size),
    QUAL = vcfPerfGrouped(perfdata.all, function(x) rowData(x)$QUAL, class.all.size),
    GQ = vcfPerfGrouped(perfdata.all, function(x) geno(x)$GQ, class.all.size),
    DP = vcfPerfGrouped(perfdata.all, function(x) info(x)$DP, class.all.size),
    FILTER = vcfPerfGrouped(perfdata.all, function(x) (rowData(x)$FILTER == "PASS")*1, class.all.size),
    "VQSLOD:FILTER" = vcfPerfGrouped(perfdata.all, function(x) info(x)$VQSLOD*(rowData(x)$FILTER == "PASS"), class.all.size),
    "QUAL:FILTER" = vcfPerfGrouped(perfdata.all, function(x) rowData(x)$QUAL*(rowData(x)$FILTER == "PASS"), class.all.size),
    "GQ:FILTER" = vcfPerfGrouped(perfdata.all, function(x) geno(x)$GQ*(rowData(x)$FILTER == "PASS"), class.all.size),
    "DP:FILTER" = vcfPerfGrouped(perfdata.all, function(x) info(x)$DP*(rowData(x)$FILTER == "PASS"), class.all.size))
}


#####################################################################
# SAVE RESULTS
#####################################################################
# For debugging: object sizes in GB
# sort(sapply(ls(), function(id) object.size(get(id)))) / 1024^3

# Remove everything but the essential data required by the report.
temp = NA
temp = ls()
temp = temp[!(grepl("^path\\.", temp) | temp %in% c("calls.sampleid", "versions", "DEBUG", "DEBUG.chrom", "genome") | grepl("^perf\\.", temp))]
rm(list = temp)

save.image("report_data.rda")
