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

param$path.tp = env$PATH_SAMPLE_OVERLAP_TP          # Output from RTG's vcfeval, overlapping
param$path.fp = env$PATH_SAMPLE_OVERLAP_FP          # param$path.test.subset and 
param$path.fn = env$PATH_SAMPLE_OVERLAP_FN          # param$path.gold.variants.subset

param$genome = env$CONST_REFERENCE_BSGENOME         # The reference genome R package
param$version = list()
param$version$script = env$CONST_VERSION_SCRIPT     # The validation reporter version string
param$version$host = env$PARAM_VERSION_EXEC_HOST    # The host system for this report.
param$version$rtg = env$PARAM_VERSION_RTG           # Software versions.
param$version$java = env$PARAM_VERSION_JAVA         #
param$version$bedtools = env$PARAM_VERSION_BEDTOOLS #

param$path.genome = env$CONST_GENOME_BEDGZ          # The full target genome

param$path.rds.output = env$PARAM_OUTPUT_RDS_PATH


#####################################################################
# LIBRARIES
#####################################################################
library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome)
library(parallel)


#####################################################################
# LOAD GENOME
#####################################################################
# The reference genome
genome.bsgenome = getBSgenome(param$genome)
genome.seqinfo = seqinfo(genome.bsgenome)
genome(genome.seqinfo) = param$genome     # To get around disagreement 
        # between the vcfs (which by default use the genome name 
        # abbreviation), and the beds (which use the full name).


#####################################################################
# LOAD VARIANTS (SPLIT INTO ERROR CLASSES FROM RTG VCFEVAL)
#####################################################################
# The call overlaps.  Load just the required fields from the VCFs.
vcf.scan_param = ScanVcfParam(geno = c("GT", "DP", "KCCG_PERF_DP_MIN"), fixed = c("ALT", "FILTER"), info = NA)

calls = list(
    tp = suppressWarnings(readVcf(TabixFile(param$path.tp), param$genome, vcf.scan_param)),
    fp = suppressWarnings(readVcf(TabixFile(param$path.fp), param$genome, vcf.scan_param)),
    fn = suppressWarnings(readVcf(TabixFile(param$path.fn), param$genome, vcf.scan_param))
)

# Simple data sanity check: do the vcfs have the same, correct sample name?
stopifnot(length(header(calls$tp)@samples) == 1 && header(calls$tp)@samples == header(calls$fp)@samples)


# BEGIN HACK:
# Get the variant sizes.  This can be done by accessing the CollapsedVCF
# structures, but is extremely slow -- in particular, access to the
# alt allele lengths is very inefficient when using the following code:
# alt_widths = sapply(alt(vcf), function(alt_options) range(width(alt_options)))
# Get around this problem using the following hack, calculating the 
# range of alternate allele lengths on the shell, and then
# augmenting the CollapsedVCF objects with these data.  These data will
# be used by classifyMutationSize.
augmentCollapsedVCFWithAltLengthRange = function(vcf, path)
{
    vcf2 = vcf
    inpipe = pipe(sprintf("gzip -dc %s | grep -v '^#' | awk 'BEGIN {FS=\"\\t\"} { nalts = split($5, alts, \",\"); min = length(alts[1]); max = length(alts[1]); for (i = 2; i <= nalts; i++) { if (length(alts[i]) < min) { min = length(alts[i]) } else if (length(alts[i]) > max) { max = length(alts[i]) } }; print length($4), min, max}'", path))
    lengths = scan(inpipe, list(ref = integer(), minalt = integer(), maxalt = integer()))
    close(inpipe)
    stopifnot(all(width(vcf) == lengths$ref))
    rowRanges(vcf2)$alt_width_min = lengths$minalt
    rowRanges(vcf2)$alt_width_max = lengths$maxalt
    return(vcf2)
}

calls$tp = augmentCollapsedVCFWithAltLengthRange(calls$tp, param$path.tp)
calls$fp = augmentCollapsedVCFWithAltLengthRange(calls$fp, param$path.fp)
calls$fn = augmentCollapsedVCFWithAltLengthRange(calls$fn, param$path.fn)
# END HACK


#####################################################################
# LOAD UNIVERSE AND ANALYSIS SUBSET
#####################################################################
universe = list(
    genome = bed2GRanges(param$path.genome, genome.seqinfo))

universe$gold_standard = bed2GRanges(param$path.gold.regions.subset, genome.seqinfo)
universe$gold_standard = intersect(universe$gold_standard, universe$genome, ignore.strand = TRUE)

# If a regions subset BED was supplied, additionally restrict the analysis 
# region to the intersection with this BED.
if (param$region.subset) {
    universe$subset = bed2GRanges(param$region.subset.path, genome.seqinfo)
    universe$subset = intersect(universe$subset, universe$genome, ignore.strand = TRUE)
} else {
    universe$subset = universe$genome
}

# Define the analysis region for set operations.
universe$analysis = intersect(universe$subset, universe$gold_standard, ignore.strand = TRUE)

# The total number of bases in the subset of the genome used for this analysis
universe_analysis_size = sum(as.numeric(width(universe$analysis)))


#####################################################################
# LOAD GENOMIC REGIONS
#####################################################################
# Various genomic regions; variants will be labelled by their
# presence or absence in these regions.
regions = list(
)

regions = lapply(regions, intersect, y = universe$genome, ignore.strand = TRUE)
regions.orig = regions

# Intersect all regions with the analysis subset.
regions = lapply(regions.orig, intersect, y = universe$analysis, ignore.strand = TRUE)


#####################################################################
# SUBSET CALLS TO ANALYSIS REGION
#####################################################################
calls = lapply(calls, function(x) x[x %within% universe$analysis])


#####################################################################
# CLASSIFY VARIANTS
#####################################################################
class = mclapply(
    list(
        zyg = classifyZygosity,                     # By zygosity
        muttype = classifyMutationType,             # By mutation type: Subst, Ins, Del, Other
        mutsize = classifyMutationSize,             # By mutation 'size' (see getMutationSizeVcf for the definition of size)
        depth = classifyDepth                       # By depth
    ), function(func) func(calls), mc.preschedule = FALSE)

# There is a bit of subtlety in this.  Some classes (zygosity, muttype,
# and mutsize) should based on the properties of the true variant, as 
# given in the gold standard.  However, in the case of FP or FN calls, 
# these properties may be wrong, as they are based on a known 
# incorrect call.  To resolve this, we manually go back and set these
# calls to the correct values from the gold standard.  For the majority 
# of classes -- anything based on regions or reference sequence 
# properties -- this correction is not necessary, as these classes are 
# based on data that are not affected by an incorrect call.

# False positives are in truth zygosity RR (hom. ref.)
class$zyg$fp$RR = class$zyg$fp$RR | TRUE
for (i in setdiff(names(class$zyg$fp), "RR"))
    class$zyg$fp[[i]] = class$zyg$fp[[i]] & FALSE

# False positives are in truth mutation type None
class$muttype$fp$None = class$muttype$fp$None | TRUE
for (i in setdiff(names(class$muttype$fp), "None"))
    class$muttype$fp[[i]] = class$muttype$fp[[i]] & FALSE

# False positives are in truth mutation size zero
class$mutsize$fp[["0-4"]] = class$mutsize$fp[["0-4"]] | TRUE
for (i in setdiff(names(class$mutsize$fp), "0-4"))
    class$mutsize$fp[[i]] = class$mutsize$fp[[i]] & FALSE


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

# Sanity check: are classes exclusive, and do they cover all 
# variants?  This is required by the later marginalization 
# simplifications.
checkClassExclusive = function(class)
{
    for (subgroup in class)
    {
        membership = Rle(FALSE, length(subgroup[[1]]))
        for (value in subgroup)
        {
            stopifnot(any(value & membership) == FALSE)
            membership = membership | value
        }
        stopifnot(all(membership) == TRUE)
    }
}


lapply(class, checkClassExclusive)



#####################################################################
# SCORE VARIABLES AND CALL CUTOFF DEFINITION
#####################################################################

criteria = list(
    "FILTER" = list(
        scoreFunc = function(x) (rowRanges(x)$FILTER == "PASS")*1,
        threshold = 0.5)
)



#####################################################################
# CALCULATE PERFORMANCE ON EVERY CLASS VALUE COMBINATION
#####################################################################
class_subsets.values = expand.grid(lapply(class, function(x) names(x[[1]])))
for (i in 1:ncol(class_subsets.values))
        class_subsets.values[,i] = ordered(class_subsets.values[,i])
class_subsets.performance_path = list()

message(sprintf("Calculating performance on %d disjoint variant subsets...", prod(sapply(class, function(x) length(names(x[[1]]))))))
temp.start_time = Sys.time()

class_subset_performance_func = function(i)
{
    if (i %% 100 == 0)
    {
        elapsed = as.numeric(difftime(Sys.time(), temp.start_time, units = "secs"))
        rate = i / elapsed
        remaining = (nrow(class_subsets.values) - i) / rate
        message(sprintf("%d / %d (%.2f%%)\tElapsed: %.0fs, est. remaining: %.0fs", i, nrow(class_subsets.values), i / nrow(class_subsets.values) * 100, elapsed, remaining))
    }

    # Calculate an indicator variable for the variants (in each of the 
    # three major categories -- tp, fp, and fn), that match the 
    # combination in class_subsets.values[i,]
    temp.indicator = sapply(names(calls), function(call_type) Rle(TRUE, nrow(calls[[call_type]])), USE.NAMES = TRUE)
    for (class_name in colnames(class_subsets.values))
    {
        temp.indicator = sapply(names(calls), function(call_type) temp.indicator[[call_type]] = temp.indicator[[call_type]] & class[[class_name]][[call_type]][[class_subsets.values[i, class_name]]])
        if (!any(sapply(temp.indicator, any)))
            break
    }

    class_subsets.variant_counts = sapply(temp.indicator, sum)

    if (all(class_subsets.variant_counts == 0))
        return(vcfPerf(data = NULL, criteria$FILTER$scoreFunc))
    else
        return(vcfPerf(
            data = list(
                vcf.tp = calls$tp[temp.indicator$tp],
                vcf.fp = calls$fp[temp.indicator$fp],
                n.fn = class_subsets.variant_counts["fn"],
                n.tn = 0),
            criteria$FILTER$scoreFunc))
}

class_subsets.performance_path = mclapply(1:nrow(class_subsets.values), class_subset_performance_func, mc.preschedule = TRUE)

class_subsets.performance_thresholded = cbind(class_subsets.values, as.data.frame(t(sapply(class_subsets.performance_path, getPerfAtCutoff, cutoff = criteria$FILTER$threshold))))

stopifnot(sum(colSums(class_subsets.performance_thresholded[,c("ntp", "nfp", "ntn", "nfn")])) == sum(sapply(calls, nrow)))
stopifnot(sum(class_subsets.performance_thresholded$nfp) + sum(class_subsets.performance_thresholded$ntn) == nrow(calls$fp))


saveRDS(
    list(
        report = list(gentime = date(), criterion = "FILTER = PASS", criterion_latex = "$\\mathrm{FILTER} = \\mathrm{PASS}$"), 
        params = param, 
        class_subsets.performance_thresholded = class_subsets.performance_thresholded, 
        regions = regions,
        regions.orig = regions.orig,
        universe = universe), 
    file = param$path.rds.output)
