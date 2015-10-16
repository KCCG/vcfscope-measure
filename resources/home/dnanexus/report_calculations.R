####################################################################
#  
#  KCCG WGS Validation Reporter -- Report calculations
#  
#  Usage: 
#    Rscript --vanilla report_calculations.R
#  
#  All parameters are passed via environment variables.  See below
#  for information on these.
#  
#  
#  Mark Pinese, 2015
#  
####################################################################

options(
    stringsAsFactors = FALSE, 
    warn = 2,                   # Treat warnings as errors -- this will be production code
    echo = TRUE)                # Emit lines to aid debugging

source("report_functions.R")



####################################################################
# LOAD ENVIRONMENT VARIABLES
####################################################################
env = as.list(Sys.getenv())

param = list()

param$path.test.orig = env$PARAM_INPUT_VCFGZ_PATH                               # Path to the original variant vcf under test

param$region.subset = env$PARAM_REGION_BED_SUPPLIED == "1"                      # Subset analysis to a set of genomic regions?
param$region.subset.path = env$PARAM_REGION_BED_PATH                            # If so, the path to the region bed file

param$path.test.subset = env$PATH_TEST_VARIANTS                                 # The test variant vcf, subset to the bed in param$region.subset.path
param$path.gold.variants.orig = env$CONST_GOLD_CALLS_VCFGZ                      # Original gold standard variant vcf
param$path.gold.regions.orig = env$CONST_GOLD_HARDMASK_VALID_REGIONS_BEDGZ      # Regions in which the gold standard vcf is reliable, as a bed.gz
param$path.gold.variants.subset = env$PATH_GOLD_VARIANTS                        # Gold standard variant vcf, subset to the bed in param$region.subset.path
param$path.gold.regions.subset = env$PATH_GOLD_REGIONS                          # Gold standard valid regions, subset to the bed in param$region.subset.path

param$sample.ids = env$PARAM_INPUT_VCF_SAMPLES
param$sample.index = as.integer(env$LOOP_SAMPLE_INDEX)
param$sample.count = as.integer(env$LOOP_NUM_SAMPLES)
param$sample.id = env$LOOP_THIS_SAMPLE_ID

param$path.tp = env$LOOP_PATH_SAMPLE_OVERLAP_TP     # Output from RTG's vcfeval, overlapping
param$path.fp = env$LOOP_PATH_SAMPLE_OVERLAP_FP     # param$path.test.subset and 
param$path.fn = env$LOOP_PATH_SAMPLE_OVERLAP_FN     # param$path.gold.variants.subset

param$extended = env$PARAM_EXTENDED == "1"          # Create an extended report?
param$genome = env$CONST_REFERENCE_BSGENOME         # The reference genome R package
param$version = list()
param$version$script = env$CONST_VERSION_SCRIPT     # The validation reporter version string
param$version$host = env$PARAM_VERSION_EXEC_HOST    # The host system for this report.
param$version$rtg = env$PARAM_VERSION_RTG           # Software versions.
param$version$java = env$PARAM_VERSION_JAVA         #
param$version$bedtools = env$PARAM_VERSION_BEDTOOLS #

param$path.function.regions.prefix = env$CONST_FUNCTIONAL_REGIONS_BEDGZ_PREFIX  # Path prefixes for functional region BEDs, and
param$path.mask.regions.prefix = env$CONST_MASK_REGIONS_BEDGZ_PREFIX            # genomic mask BEDs.

param$path.rds.output = env$PARAM_OUTPUT_RDS_PATH


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
genome.bsgenome = getBSgenome(param$genome)
genome.seqinfo = seqinfo(genome.bsgenome)
genome(genome.seqinfo) = param$genome     # To get around disagreement 
        # between the vcfs (which by default use the genome name 
        # abbreviation), and the beds (which use the full name).

# The call overlaps.
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

# Load just the required fields from the VCFs.  Which fields to
# fetch depends on the vcf type, as TP and FP vcf entries are sourced
# from the call VCF, whereas FN entries are sourced from the gold
# standard VCF.
vcf.scan_param.called = ScanVcfParam(info = c("DP", "GQ_MEAN", "QD", "VQSLOD", "MQ", "FS", "MQRankSum", "ReadPosRankSum"), geno = c("GT", "DP", "GQ"))
vcf.scan_param.uncalled = ScanVcfParam(info = c("DP", "TYPE"), geno = c("GT", "DP", "GQ"))

calls = list(
    tp = suppressWarnings(readVcf(TabixFile(param$path.tp), param$genome, vcf.scan_param.called)),
    fp = suppressWarnings(readVcf(TabixFile(param$path.fp), param$genome, vcf.scan_param.called)),
    fn = suppressWarnings(readVcf(TabixFile(param$path.fn), param$genome, vcf.scan_param.uncalled))
)


# Simple data sanity check: do the vcfs 
# have the same sample name?
temp.sample.tp = header(calls$tp)@samples
temp.sample.fp = header(calls$fp)@samples
stopifnot(temp.sample.tp == temp.sample.fp)
stopifnot(length(temp.sample.tp) == 1)

calls.sampleid = header(calls$tp)@samples


# Various genomic regions, for later subsetting of performance measures
regions = list(
    gold = list(callable = bed2GRanges(param$path.gold.regions.subset, genome.seqinfo)),      # Gold standard valid call regions
    mask = readMaskRegions(param$path.mask.regions.prefix, genome.seqinfo),                   # Masking beds
    functional = readFunctionalRegions(param$path.function.regions.prefix, genome.seqinfo),   # 'Function classes' of the genome
    universe = list(genome = GRanges(                                                         # The whole genome
        seqnames = seqnames(genome.seqinfo), 
        ranges = IRanges(1, seqlengths(genome.seqinfo)), 
        strand = "*", seqinfo = genome.seqinfo)))

# The universe for set operations.  If a regions subset BED was supplied, 
# the universe is this set of regions, and so load it.  Otherwise, the
# universe is the whole genome.
if (param$region.subset) {
    regions$universe$universe = intersect(regions$universe$genome, bed2GRanges(param$region.subset.path, genome.seqinfo))     # Universe is the supplied region BED file.
} else {
    regions$universe$universe = regions$universe$genome                                       # Universe is the whole genome
}

# Intersect all regions with the universe.
for (i in names(regions))
{
    if (i != "genome")
    {
        for (j in names(regions[[i]]))
            regions[[i]][[j]] = intersect(regions[[i]][[j]], regions$universe$universe)
    }
}

# Create new mask classes, of "unmasked" and "not gold standard callable"
regions$mask$unmasked = setdiff(regions$universe$universe, reduce(union(union(regions$mask$ambiguous, regions$mask$low_complexity), regions$mask$repetitive)))
regions$gold$notcallable = setdiff(regions$universe$universe, regions$gold$callable)


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
    functional = classifyRegionOverlap(calls, regions$functional, c("coding" = "Any", "genic" = "Any", "splice" = "Any", "intronic" = "All", "utr" = "All", "intergenic" = "All")),

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

# Sanity check: are the classes that should be exclusive, in 
# fact exclusive?  An example is mutation type, in which each
# vcf entry should be correspond to exactly one type of mutation.
checkClassExclusive = function(class, exact = TRUE)
{
    groups = names(class[[1]])
    for (group in groups)
    {
        if (length(class[[1]][[group]]) == 0)
            next
        indicators = sapply(class, function(subclass) as.vector(subclass[[group]]))
        if (is.vector(indicators))  { indicators = matrix(indicators, nrow = 1, ncol = length(indicators)) }
        if (exact)
            stopifnot(all(rowSums(indicators) == 1))
        else
            stopifnot(all(rowSums(indicators) <= 1))
    }
}
checkClassExclusive(class$somy, exact = FALSE)      # exact=FALSE, as 'chromosomes' such as GL000207.1 have no somy assignation by this code
checkClassExclusive(class$muttype)
checkClassExclusive(class$mutsize)
checkClassExclusive(class$goldcall, exact = FALSE)  # exact=FALSE, as a variant can be neither entirely inside a gold-callable region, or entirely outside of one (ie it may be only partially in a callable region)



#####################################################################
# SCORE VARIABLES AND CALL CUTOFF DEFINITIONS
#####################################################################

# If some of the fields required by a given criterion are missing, the return value
# will be NULL, and consequently downstream performance estimation functions (vcfPerf
# and vcfPerfGrouped) will deliberately produce a no-information performance estimate.
criteria = list(
    "VQSLOD" =          list(scoreFunc = function(x) info(x)$VQSLOD),
    "QUAL" =            list(scoreFunc = function(x) rowRanges(x)$QUAL),
    "GQ" =              list(scoreFunc = function(x) info(x)$GQ),
    "GL" =              list(scoreFunc = function(x) info(x)$GL),
    "FILTER" =          list(scoreFunc = function(x) (rowRanges(x)$FILTER == "PASS")*1),
    "DEPTH" =           list(scoreFunc = function(x) info(x)$DP)
)



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

# We can define TNs for SNVs
calls.snv$tn = setdiff(subset.snv.regions, union(rowRanges(calls.snv$tp), rowRanges(calls.snv$fp), rowRanges(calls.snv$fn), ignore.strand = TRUE))

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

if (param$extended) {
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

perf.indelsubst$zyg = lapply(criteria, function(crit) vcfPerfGrouped(perfdata.indelsubst, crit$scoreFunc, class.indelsubst.zyg))

if (param$extended) {
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
if (param$extended) {
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
} else {
    perf.all = NULL
}


#####################################################################
# THRESHOLDED PERFORMANCE SUMMARY CALCULATIONS FOR REPORT
#####################################################################
report = list()

report$gentime = date()

report$specifications = list(
    measure = "FILTER",
    cutoff = 0.5,
    snv.het.min.sens = 0.95,
    snv.het.min.spec = 0.95,
    snv.hom.min.sens = 0.99,
    snv.hom.min.spec = 0.99,
    indelsubst.het.min.sens = 0.80,
    indelsubst.hom.min.sens = 0.95)
report$specifications$label = "FILTER = PASS"
report$specifications$label_latex = "$\\mathrm{FILTER} = \\mathrm{PASS}$"

snv.het.perf = calcSensSpecAtCutoff(perf.snv$zyg[[report$specifications$measure]]$RRvsRA, report$specifications$cutoff)
indelsubst.het.perf = calcSensSpecAtCutoff(perf.indelsubst$zyg[[report$specifications$measure]]$RRvsRA, report$specifications$cutoff)
snv.hom.perf = calcSensSpecAtCutoff(perf.snv$zyg[[report$specifications$measure]]$RRvsAA, report$specifications$cutoff)
indelsubst.hom.perf = calcSensSpecAtCutoff(perf.indelsubst$zyg[[report$specifications$measure]]$RRvsAA, report$specifications$cutoff)

report$snv.het.sens.value = snv.het.perf$sens
report$snv.het.spec.value = snv.het.perf$spec
report$indelsubst.het.sens.value = indelsubst.het.perf$sens
report$snv.hom.sens.value = snv.hom.perf$sens
report$snv.hom.spec.value = snv.hom.perf$spec
report$indelsubst.hom.sens.value = indelsubst.hom.perf$sens

report$snv.het.sens.pass = report$snv.het.sens.value >= report$specifications$snv.het.min.sens
report$snv.het.spec.pass = report$snv.het.spec.value >= report$specifications$snv.het.min.spec
report$indelsubst.het.sens.pass = report$indelsubst.het.sens.value >= report$specifications$indelsubst.het.min.sens
report$snv.hom.sens.pass = report$snv.hom.sens.value >= report$specifications$snv.hom.min.sens
report$snv.hom.spec.pass = report$snv.hom.spec.value >= report$specifications$snv.hom.min.spec
report$indelsubst.hom.sens.pass = report$indelsubst.hom.sens.value >= report$specifications$indelsubst.hom.min.sens

report$overall.pass = 
    report$snv.het.sens.pass & report$snv.het.spec.pass & report$indelsubst.het.sens.pass &
    report$snv.hom.sens.pass & report$snv.hom.spec.pass & report$indelsubst.hom.sens.pass

report$snv.het.sens.call = ifelse(report$snv.het.sens.pass, "Pass", "\\textcolor{red}{\\textbf{FAIL}}")
report$snv.het.spec.call = ifelse(report$snv.het.spec.pass, "Pass", "\\textcolor{red}{\\textbf{FAIL}}")
report$indelsubst.het.sens.call = ifelse(report$indelsubst.het.sens.pass, "Pass", "\\textcolor{red}{\\textbf{FAIL}}")
report$snv.hom.sens.call = ifelse(report$snv.hom.sens.pass, "Pass", "\\textcolor{red}{\\textbf{FAIL}}")
report$snv.hom.spec.call = ifelse(report$snv.hom.spec.pass, "Pass", "\\textcolor{red}{\\textbf{FAIL}}")
report$indelsubst.hom.sens.call = ifelse(report$indelsubst.hom.sens.pass, "Pass", "\\textcolor{red}{\\textbf{FAIL}}")

report$overall.call = ifelse(report$overall.pass, "\\textcolor{blue}{PASS}", "\\textcolor{red}{\\textbf{FAIL}}")


#####################################################################
# EXPORT SUMMARY RESULTS AS RDS
#####################################################################
export = list(
    param = param,
    performance = list(
        whole_genome = list(
            combined = perf.all,
            snv = perf.snv,
            indelsubst = perf.indelsubst
        )
    ),
    report_summary = report
)

saveRDS(export, file = "report_summary.rds", version = 2, compress = "xz")


#####################################################################
# SAVE FULL RESULTS FOR THE REPORT
#####################################################################
save.image("report_data.rda")
