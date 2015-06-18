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


# 1000 genomes SNV / Indel filter, as described by Brad Chapman.
# SNV:      Fail if (QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0)
# Indel:    Fail if (QD < 2.0 ||              FS > 200.0 ||                     ReadPosRankSum < -20.0)
filter_1000G = function(vcf)
{
    # Make sure the vcf has all the metrics used by the filter.
    stopifnot(length(setdiff("QD", "MQ", "FS", "MQRankSum", "ReadPosRankSum"), colnames(info(vcf))) == 0)

    snv = isSNV(vcf)

    # Calculate PASS logicals for each metric, splitting by snv / indel as appropriate.
    QD.pass = info(vcf)$QD >= 2.0
    MQ.pass = info(vcf)$MQ >= ifelse(snv, 40.0, -Inf)
    FS.pass = info(vcf)$FS <= ifelse(snv, 60.0, 200.0)
    MQRankSum.pass = info(vcf)$MQRankSum >= ifelse(snv, -12.5, -Inf)
    ReadPosRankSum.pass = info(vcf)$ReadPosRankSum >= ifelse(snv, -8.0, -20.0)

    # If the desired values are missing, assume they pass
    QD.pass[is.na(info(vcf)$QD)] = TRUE
    MQ.pass[is.na(info(vcf)$MQ)] = TRUE
    FS.pass[is.na(info(vcf)$FS)] = TRUE
    MQRankSum.pass[is.na(info(vcf)$MQRankSum)] = TRUE
    ReadPosRankSum.pass[is.na(info(vcf)$ReadPosRankSum)] = TRUE

    # Final pass is the AND of all individual filters.
    QD.pass & MQ.pass & FS.pass & MQRankSum.pass & ReadPosRankSum.pass
}


criteria = list(
    "VQSLOD" =          list(scoreFunc = function(x) info(x)$VQSLOD,                                    callFunc = function(x) info(x)$VQSLOD > 2.7),
    "QUAL" =            list(scoreFunc = function(x) rowData(x)$QUAL,                                   callFunc = function(x) rowData(x)$QUAL > 200),
    "FILTER" =          list(scoreFunc = function(x) (rowData(x)$FILTER == "PASS")*1,                   callFunc = function(x) rowData(x)$FILTER == "PASS"),
    "VQSLOD:1000G" =    list(scoreFunc = function(x) info(x)$VQSLOD * filter_1000G(x),                  callFunc = function(x) (info(x)$VQSLOD > 2.7) * filter_1000G(x)),
    "QUAL:1000G" =      list(scoreFunc = function(x) rowData(x)$QUAL * filter_1000G(x),                 callFunc = function(x) (rowData(x)$QUAL > 200) * filter_1000G(x)),
    "FILTER:1000G" =    list(scoreFunc = function(x) (rowData(x)$FILTER == "PASS") * filter_1000G(x),   callFunc = function(x) (rowData(x)$FILTER == "PASS") & filter_1000G(x))
)



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
    stopifnot(length(setdiff("QD", ), colnames(info(vcf))) == 0)
    vcf.scan_param.called = ScanVcfParam(info = c("DP", "GQ_MEAN", "QD", "VQSLOD", "MQ", "FS", "MQRankSum", "ReadPosRankSum"), geno = c("GT", "DP", "GQ"), which = DEBUG.region)
    vcf.scan_param.uncalled = ScanVcfParam(info = c("DP", "TYPE"), geno = c("GT", "DP", "GQ"), which = DEBUG.region)
} else {
    vcf.scan_param.called = ScanVcfParam(info = c("DP", "GQ_MEAN", "QD", "VQSLOD", "MQ", "FS", "MQRankSum", "ReadPosRankSum"), geno = c("GT", "DP", "GQ"))
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

perf.snv$zyg = lapply(criteria, function(crit) vcfPerfGrouped(perfdata.snv, crit$scoreFunc, class.snv.zyg))

if (extendedflag) {
perf.snv$all = lapply(criteria, function(crit) vcfPerf(perfdata.snv, crit$scoreFunc))

# And sequence context
# The nasty code for counting TNs just counts, for each class in 
# regions$mask, the number of TN bases (from calls.snv$tn) that
# overlap this mask class.
class.snv.mask = subsetClass(class$mask, subset.snv, tn = lapply(regions$mask, function(mask) sum(as.numeric(width(intersect(mask, calls.snv$tn))))))
perf.snv$mask = lapply(criteria, function(crit) vcfPerfGrouped(perfdata.snv, crit$scoreFunc, class.snv.mask))
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

perf.indelsubst$zyg = lapply(criteria, function(crit) vcfPerfGrouped(perfdata.indelsubst, crit$scoreFunc, class.indelsubst.zyg))

if (extendedflag) {
perf.indelsubst$all = lapply(criteria, function(crit) vcfPerf(perfdata.indelsubst, crit$scoreFunc))

# And sequence context
# Again, set tn = NULL
class.indelsubst.mask = subsetClass(class$mask, subset.indelsubst, tn = NULL)
perf.indelsubst$mask = lapply(criteria, function(crit) vcfPerfGrouped(perfdata.indelsubst, crit$scoreFunc, class.indelsubst.mask))

# And mutation size
# Again, set tn = NULL
class.indelsubst.size = subsetClass(class$mutsize, subset.indelsubst, tn = NULL)
perf.indelsubst$mutsize = lapply(criteria, function(crit) vcfPerfGrouped(perfdata.indelsubst, crit$scoreFunc, class.indelsubst.size))
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
perf.all$all = lapply(criteria, function(crit) vcfPerf(perfdata.all, crit$scoreFunc))

# Subset by zygosity
# We can't count TNs in an all variant context; just set to a null value
class.all.zyg = subsetClass(class$zyg, subset.all, tn = NULL)
# Basic checks for subsetClass
stopifnot(all(class$zyg$RRvsRA$tp[subset.all$tp] == class.all.zyg$RRvsRA$tp))
stopifnot(all(class$zyg$RRvsRA$fp[subset.all$fp] == class.all.zyg$RRvsRA$fp))
stopifnot(all(class$zyg$RRvsAA$tp[subset.all$tp] == class.all.zyg$RRvsAA$tp))

perf.all$zyg = lapply(criteria, function(crit) vcfPerfGrouped(perfdata.all, crit$scoreFunc, class.all.zyg))

# And sequence context
# Again, set tn = NULL
class.all.mask = subsetClass(class$mask, subset.all, tn = NULL)
perf.all$mask = lapply(criteria, function(crit) vcfPerfGrouped(perfdata.all, crit$scoreFunc, class.all.mask))

# And mutation type
# Again, set tn = NULL
class.all.type = subsetClass(class$muttype, subset.all, tn = NULL)
perf.all$muttype = lapply(criteria, function(crit) vcfPerfGrouped(perfdata.all, crit$scoreFunc, class.all.type))

# And mutation size
# Again, set tn = NULL
class.all.size = subsetClass(class$mutsize, subset.all, tn = NULL)
perf.all$mutsize = lapply(criteria, function(crit) vcfPerfGrouped(perfdata.all, crit$scoreFunc, class.all.size))
}


#####################################################################
# SNV METRIC PERFORMANCE: GOLD CALLABLE, CODING +/- 10 REGIONS
#####################################################################
# Subset to SNVs in gold-callable regions within 10 bp of CDs
# NOTE: subset.snv.coding10.regions, and subset.snv.coding10, MUST MATCH
subset.snv.coding10.regions = intersect(regions$gold$callable, regions$functional$coding_10)
subset.snv.coding10 = sapply(names(calls), function(name) class$muttype$SNV[[name]] & class$goldcall$callable[[name]] & class$functional$coding_10[[name]], simplify = FALSE, USE.NAMES = TRUE)

calls.snv.coding10 = sapply(names(calls), function(name) calls[[name]][subset.snv.coding10[[name]]], simplify = FALSE, USE.NAMES = TRUE)
calls.snv.coding10$tn = setdiff(subset.snv.coding10.regions, union(rowData(calls.snv.coding10$tp), rowData(calls.snv.coding10$fp), rowData(calls.snv.coding10$fn), ignore.strand = TRUE))

count.snv.coding10.fn = nrow(calls.snv.coding10$fn)
count.snv.coding10.tn = sum(as.numeric(width(calls.snv.coding10$tn)))

perfdata.snv.coding10 = list(vcf.tp = calls.snv.coding10$tp, vcf.fp = calls.snv.coding10$fp, n.fn = count.snv.coding10.fn, n.tn = count.snv.coding10.tn)

class.snv.coding10.zyg = subsetClass(class$zyg, subset.snv.coding10, tn = list(RRvsAA = count.snv.coding10.tn, RRvsRA = count.snv.coding10.tn, RRvsAB = count.snv.coding10.tn))
perf.snv.coding10 = list(zyg = lapply(criteria, function(crit) vcfPerfGrouped(perfdata.snv.coding10, crit$scoreFunc, class.snv.coding10.zyg)))


#####################################################################
# INDEL AND SUBSTITUTION METRIC PERFORMANCE: GOLD CALLABLE, CODING +/- 10 REGIONS
#####################################################################
# Subset to SNVs in gold-callable regions within 10 bp of CDs
# NOTE: subset.indelsubst.coding10.regions, and subset.indelsubst.coding10, MUST MATCH
subset.indelsubst.coding10 = sapply(names(calls), function(name) class$muttype$InsDelSubst[[name]] & class$goldcall$callable[[name]] & class$functional$coding_10[[name]], simplify = FALSE, USE.NAMES = TRUE)

calls.indelsubst.coding10 = sapply(names(calls), function(name) calls[[name]][subset.indelsubst.coding10[[name]]], simplify = FALSE, USE.NAMES = TRUE)

count.indelsubst.coding10.fn = nrow(calls.indelsubst.coding10$fn)

perfdata.indelsubst.coding10 = list(vcf.tp = calls.indelsubst.coding10$tp, vcf.fp = calls.indelsubst.coding10$fp, n.fn = count.indelsubst.coding10.fn, n.tn = 0)

class.indelsubst.coding10.zyg = subsetClass(class$zyg, subset.indelsubst.coding10, tn = NULL)
perf.indelsubst.coding10 = list(zyg = lapply(criteria, function(crit) vcfPerfGrouped(perfdata.indelsubst.coding10, crit$scoreFunc, class.indelsubst.coding10.zyg)))


#####################################################################
# SAVE RESULTS
#####################################################################
# For debugging: object sizes in GB
# sort(sapply(ls(), function(id) object.size(get(id)))) / 1024^3

# Remove everything but the essential data required by the report.
temp = NA
temp = ls()
temp = temp[!(grepl("^path\\.", temp) | temp %in% c("calls.sampleid", "versions", "DEBUG", "DEBUG.chrom", "genome", "extendedflag") | grepl("^perf\\.", temp))]
rm(list = temp)

save.image("report_data.rda")
