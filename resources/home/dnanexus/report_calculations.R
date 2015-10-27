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

param$genome = env$CONST_REFERENCE_BSGENOME         # The reference genome R package
param$version = list()
param$version$script = env$CONST_VERSION_SCRIPT     # The validation reporter version string
param$version$host = env$PARAM_VERSION_EXEC_HOST    # The host system for this report.
param$version$rtg = env$PARAM_VERSION_RTG           # Software versions.
param$version$java = env$PARAM_VERSION_JAVA         #
param$version$bedtools = env$PARAM_VERSION_BEDTOOLS #

param$path.rmsk = env$CONST_RMSK_REGIONS_BEDGZ      # Repeatmasker masked regions
param$path.mdust = env$CONST_MDUST_REGIONS_BEDGZ    # mdust low-complexity regions
param$path.genome = env$CONST_GENOME_BEDGZ          # The full target genome

param$path.rds.output = env$PARAM_OUTPUT_RDS_PATH


#####################################################################
# LIBRARIES
#####################################################################
library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome)


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
vcf.scan_param = ScanVcfParam(geno = c("GT"))

calls = list(
    tp = suppressWarnings(readVcf(TabixFile(param$path.tp), param$genome, vcf.scan_param)),
    fp = suppressWarnings(readVcf(TabixFile(param$path.fp), param$genome, vcf.scan_param)),
    fn = suppressWarnings(readVcf(TabixFile(param$path.fn), param$genome, vcf.scan_param))
)

# Simple data sanity check: do the vcfs have the same sample name?
temp.sample.tp = header(calls$tp)@samples
temp.sample.fp = header(calls$fp)@samples
stopifnot(temp.sample.tp == temp.sample.fp)
stopifnot(length(temp.sample.tp) == 1)

calls.sampleid = header(calls$tp)@samples


#####################################################################
# LOAD UNIVERSE AND ANALYSIS SUBSET
#####################################################################
universe = list(
    genome = bed2GRanges(param$path.genome, genome.seqinfo),
    gold_standard = bed2GRanges(param$path.gold.regions.subset, genome.seqinfo))

# If a regions subset BED was supplied, additionally restrict the analysis 
# region to the intersection with this BED.
if (param$region.subset) {
    universe$subset = bed2GRanges(param$region.subset.path, genome.seqinfo)
} else {
    universe$subset = universe$genome
}

# Define the analysis region for set operations.
universe$analysis = intersect(universe$subset, intersect(universe$genome, universe$gold_standard, ignore.strand = TRUE), ignore.strand = TRUE)

# The total number of bases in the subset of the genome used for this analysis
universe_analysis_size = sum(as.numeric(width(universe$analysis)))


#####################################################################
# LOAD GENOMIC REGIONS
#####################################################################
# Various genomic regions; variants will be labelled by their
# presence or absence in these regions.
regions = list(
    rmsk = bed2GRanges(param$path.rmsk, genome.seqinfo),
    mdust = bed2GRanges(param$path.mdust, genome.seqinfo))

# Intersect all regions with the analysis subset.
regions = lapply(regions, intersect, y = universe$analysis, ignore.strand = TRUE)


#####################################################################
# SUBSET CALLS TO ANALYSIS REGION
#####################################################################
calls = lapply(calls, function(x) x[x %within% universe$analysis])


#####################################################################
# CLASSIFY VARIANTS
#####################################################################
class = list(
    zyg = classifyZygosity(calls),                  # By zygosity
    muttype = classifyMutationType(calls),          # By mutation type: Subst, Ins, Del, Other
    mutsize = classifyMutationSize(calls),          # By mutation 'size' (see getMutationSizeVcf for the definition of size)

    # By any overlap with Repeatmasker masked regions
    rmsk = classifyRegionOverlap(calls, 
        list(masked = regions$rmsk, unmasked = setdiff(universe$analysis, regions$rmsk, ignore.strand = TRUE)), 
        c("masked" = "Any", "unmasked" = "All")),
    
    # By any overlap with mdust marked low-complexity regions
    mdust = classifyRegionOverlap(calls, 
        list(masked = regions$mdust, unmasked = setdiff(universe$analysis, regions$mdust, ignore.strand = TRUE)), 
        c("masked" = "Any", "unmasked" = "All"))
)

# TODO: Set zyg, muttype, and mutsize by reference variants.  TP is fine, but for FN need
# aligned true variants.  And for FP, need to hard code.

# There is a bit of subtlety in this.  Some classes (zygosity, muttype,
# and mutsize), should based on the properties of the true variant, as 
# given in the gold standard.  However, in the case of FP or FN calls, 
# these properties may be wrong, as they are based on a known 
# incorrect call.  To resolve this, we manually go back and set these
# calls to the correct values from the gold standard.  For the majority 
# of classes -- anything based on regions or reference sequence 
# properties -- this correction is not necessary, as these classes are 
# based on data that are not affected by an incorrect call.

# False positives are in truth zygosity RR (hom. ref.)
class$zyg$RR$fp = class$zyg$RR$fp | TRUE
for (i in setdiff(names(class$zyg), "RR"))
    class$zyg[[i]]$fp = class$zyg[[i]]$fp & FALSE

# False positives are in truth mutation type None
class$muttype$None$fp = class$muttype$None$fp | TRUE
for (i in setdiff(names(class$muttype), "None"))
    class$muttype[[i]]$fp = class$muttype[[i]]$fp & FALSE

# False positives are in truth mutation size zero
class$mutsize[["0"]]$fp = class$mutsize[["0"]]$fp | TRUE
for (i in setdiff(names(class$mutsize), "0"))
    class$mutsize[[i]]$fp = class$mutsize[[i]]$fp & FALSE


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
    ngroups = length(class)
    if (ngroups == 1)
        return()

    bins = names(class[[1]])

    for (bin in bins)
    {
        x = class[[1]][[bin]]
        for (i in 2:ngroups)
        {
            y = class[[i]][[bin]]
            stopifnot(any(x & y) == FALSE)
            x = x | y
        }

        stopifnot(all(x) == TRUE)
    }
}

checkClassExclusive(class$zyg)
checkClassExclusive(class$muttype)
checkClassExclusive(class$mutsize)
checkClassExclusive(class$rmsk)
checkClassExclusive(class$mdust)



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
class_subsets.values = expand.grid(lapply(class, names))
class_subsets.variant_counts = matrix(0, nrow = nrow(class_subsets.values), ncol = 3)
colnames(class_subsets.variant_counts) = names(calls)
class_subsets.performance_path = list()

message("Calculating performance on disjoint variant subsets...")
temp.progress = txtProgressBar(max = nrow(class_subsets.values), style = 3)
for (i in 1:nrow(class_subsets.values))
{
    # Calculate an indicator variable for the variants (in each of the 
    # three major categories -- tp, fp, and fn), that match the 
    # combination in class_subsets.values[i,]
    temp.indicator = sapply(names(calls), function(call_type) Rle(TRUE, length(calls[[call_type]])), USE.NAMES = TRUE)
    for (class_name in colnames(class_subsets.values))
    {
        temp.indicator = sapply(names(calls), function(call_type) temp.indicator[[call_type]] = temp.indicator[[call_type]] & class[[class_name]][[class_subsets.values[i, class_name]]][[call_type]])
    }

    # Tally the number of variants that fall into this subset, for
    # later sanity checking.    
    class_subsets.variant_counts[i,] = sapply(temp.indicator, sum)[colnames(class_subsets.variant_counts)]

    if (all(class_subsets.variant_counts[i,] == 0))
        class_subsets.performance_path[[i]] = vcfPerf(data = NULL, criteria$FILTER$scoreFunc)
    else
        class_subsets.performance_path[[i]] = vcfPerf(
            data = list(
                vcf.tp = calls$tp[temp.indicator$tp],
                vcf.fp = calls$fp[temp.indicator$fp],
                n.fn = class_subsets.variant_counts[i, "fn"],
                n.tn = 0),
            criteria$FILTER$scoreFunc)

    setTxtProgressBar(temp.progress, i)
}
stopifnot(colSums(class_subsets.variant_counts) == sapply(calls, nrow)[colnames(class_subsets.variant_counts)])

class_subsets.performance_thresholded = cbind(class_subsets.values, as.data.frame(t(sapply(class_subsets.performance_path, getPerfAtCutoff, cutoff = criteria$FILTER$threshold))))

stopifnot(rowSums(class_subsets.performance_thresholded[,c("ntp", "nfp", "ntn", "nfn")]) == rowSums(class_subsets.variant_counts))
stopifnot(sum(colSums(class_subsets.performance_thresholded[,c("ntp", "nfp", "ntn", "nfn")])) == sum(sapply(calls, nrow)))
stopifnot(sum(class_subsets.performance_thresholded$nfp) + sum(class_subsets.performance_thresholded$ntn) == nrow(calls$fp))


#####################################################################
# MARGINALIZE THRESHOLDED PERFORMANCE FOR PLOTS
#####################################################################
universe_analysis_size

library(ggplot2)
ggplot(class_subsets.performance_thresholded, aes(x = mutsize, y = ntp / (ntp + nfn), fill = zyg)) + geom_bar(stat = "sum") + facet_grid(muttype ~ rmsk + mdust) + theme_bw()

library(plyr)


marginalizePerformance = function(perf_data, subset, vars, ...)
{
    ddply(perf_data[eval(subset, perf_data, parent.frame()),], vars, function(x) { 
        ntp = sum(x$ntp)
        nfn = sum(x$nfn)
        nfp = sum(x$nfp)
        ntn = sum(x$ntn)
        if (ntp + nfn + nfp + ntn == 0)
            return(NULL)
        ci_test = binom.test(ntp, ntp + nfn, ...)

        sens = ci_test$estimate
        sens.lci = ci_test$conf.int[1]
        sens.uci = ci_test$conf.int[2]

        c("sens" = sens, "sens.lci" = sens.lci, "sens.uci" = sens.uci)
    })
}


ggplot(
    marginalizePerformance(class_subsets.performance_thresholded, quote(muttype == "Subst"), .(zyg, mdust, rmsk)), 
    aes(x = zyg, y = sens, fill = zyg)) + 
    geom_bar(stat = "identity", position = "dodge") + facet_grid(mdust ~ rmsk) + theme_bw()

ggplot(
    marginalizePerformance(class_subsets.performance_thresholded, quote(muttype == "Subst"), .(zyg)), 
    aes(x = zyg, y = sens, fill = zyg)) + 
    geom_bar(stat = "identity", position = "dodge") + theme_bw()

temp.perf = marginalizePerformance(class_subsets.performance_thresholded, quote(muttype == "Ins"), .(zyg, mutsize, mdust, rmsk))
temp.maxsize = max(as.numeric(gsub("+", "", unique(temp.perf$mutsize))))
temp.perf$mutsize = ordered(as.vector(temp.perf$mutsize), levels = c(as.character(0:(temp.maxsize-1)), paste(temp.maxsize, "+", sep = ""))), 
ggplot(
    temp.perf, 
    aes(x = mutsize, y = sens, fill = zyg)) + 
    geom_bar(stat = "identity", position = "dodge") + facet_grid(mdust ~ rmsk) + theme_bw()

temp.perf = marginalizePerformance(class_subsets.performance_thresholded, quote(muttype == "Ins"), .(zyg, mutsize))
temp.maxsize = max(as.numeric(gsub("+", "", unique(temp.perf$mutsize))))
temp.perf$mutsize = ordered(as.vector(temp.perf$mutsize), levels = c(as.character(0:(temp.maxsize-1)), paste(temp.maxsize, "+", sep = ""))), 
ggplot(
    temp.perf, 
    aes(x = mutsize, y = sens, fill = zyg)) + 
    geom_bar(stat = "identity", position = "dodge") + theme_bw()

temp.perf = marginalizePerformance(class_subsets.performance_thresholded, quote(muttype == "Del"), .(zyg, mutsize, mdust, rmsk)), 
temp.maxsize = max(as.numeric(gsub("+", "", unique(temp.perf$mutsize))))
temp.perf$mutsize = ordered(as.vector(temp.perf$mutsize), levels = c(as.character(0:(temp.maxsize-1)), paste(temp.maxsize, "+", sep = ""))), 
ggplot(
    temp.perf,
    aes(x = mutsize, y = sens, fill = zyg)) + 
    geom_bar(stat = "identity", position = "dodge") + facet_grid(mdust ~ rmsk) + theme_bw()

temp.perf = marginalizePerformance(class_subsets.performance_thresholded, quote(muttype == "Del"), .(zyg, mutsize)), 
temp.maxsize = max(as.numeric(gsub("+", "", unique(temp.perf$mutsize))))
temp.perf$mutsize = ordered(as.vector(temp.perf$mutsize), levels = c(as.character(0:(temp.maxsize-1)), paste(temp.maxsize, "+", sep = ""))), 
ggplot(
    temp.perf,
    aes(x = mutsize, y = sens, fill = zyg)) + 
    geom_bar(stat = "identity", position = "dodge") + theme_bw()

false_positive_counts = ddply(subset(class_subsets.performance_thresholded, muttype == "None"), .(mdust, rmsk), function(x) c(nfp = sum(x$nfp)))
false_positive_counts$subset_size = NA
false_positive_counts$subset_size[false_positive_counts$rmsk == "masked" & false_positive_counts$mdust == "masked"] = sum(as.numeric(width(intersect(regions$rmsk, regions$mdust, ignore.strand = TRUE))))
false_positive_counts$subset_size[false_positive_counts$rmsk == "masked" & false_positive_counts$mdust == "unmasked"] = sum(as.numeric(width(intersect(regions$rmsk, setdiff(universe$analysis, regions$mdust, ignore.strand = TRUE), ignore.strand = TRUE))))
false_positive_counts$subset_size[false_positive_counts$rmsk == "unmasked" & false_positive_counts$mdust == "masked"] = sum(as.numeric(width(intersect(setdiff(universe$analysis, regions$rmsk, ignore.strand = TRUE), regions$mdust, ignore.strand = TRUE))))
false_positive_counts$subset_size[false_positive_counts$rmsk == "unmasked" & false_positive_counts$mdust == "unmasked"] = sum(as.numeric(width(intersect(setdiff(universe$analysis, regions$rmsk, ignore.strand = TRUE), setdiff(universe$analysis, regions$mdust, ignore.strand = TRUE), ignore.strand = TRUE))))
stopifnot(sum(false_positive_counts$subset_size) == universe_analysis_size)
false_positive_counts$rate_per_Mb = false_positive_counts$nfp / false_positive_counts$subset_size * 1e6

sum(false_positive_counts$nfp) / universe_analysis_size * 1e6

false_positive_counts$class = NA
false_positive_counts$class[false_positive_counts$rmsk == "masked" & false_positive_counts$mdust == "masked"] = "Both"
false_positive_counts$class[false_positive_counts$rmsk == "masked" & false_positive_counts$mdust == "unmasked"] = "RMSK"
false_positive_counts$class[false_positive_counts$rmsk == "unmasked" & false_positive_counts$mdust == "masked"] = "mdust"
false_positive_counts$class[false_positive_counts$rmsk == "unmasked" & false_positive_counts$mdust == "unmasked"] = "Neither"
false_positive_counts$class = ordered(false_positive_counts$class, levels = c("Both", "mdust", "RMSK", "Neither"))

ggplot(false_positive_counts, aes(x = class, y = rate_per_Mb)) + geom_bar(stat = "identity")
sum(false_positive_counts$rate_per_Mb)
