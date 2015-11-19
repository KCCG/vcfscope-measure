library(inline)
library(Rcpp)


bed2GRanges = function(path, seqinfo)
{
    # Read first line, to get the number of columns
    firstline = read.table(path, header = FALSE, stringsAsFactors = FALSE, nrows = 1)
    ncols = length(firstline)

    if (ncols > 6)
        warning("BED file %s contains more than 6 columns; discarding columns 7 and above.")
    else if (ncols < 3)
        stop(sprintf("%s is not a valid BED file: fewer than three columns.", path))

    header = c("chr", "start", "end", "id", "score", "strand")
    column_types = list(character(), integer(), integer(), character(), double(), character())

    if (ncols < length(column_types))
    {
        column_types = column_types[1:ncols]
    }
    else if (ncols > length(column_types))
    {
        for (i in (length(column_types)+1):ncols)
            column_types[i] = NULL
    }

    data = scan(path, what = column_types, strip.white = TRUE)
    data = data.frame(matrix(unlist(data), ncol = ncols, byrow = FALSE), stringsAsFactors = FALSE)
    colnames(data) = header[1:ncols]

    data$start = as.integer(data$start)
    data$end = as.integer(data$end)
    if (ncols >= 5)
        data$score = as.numeric(data$score)

    data$start = data$start + 1

    data = data[data$chr %in% seqnames(seqinfo),]
    
    if (ncol(data) == 3)
        result = GRanges(data$chr, IRanges(data$start, data$end), seqinfo = seqinfo)
    else if (ncol(data) == 4)
        result = GRanges(data$chr, IRanges(data$start, data$end), id = data$id, seqinfo = seqinfo)
    else if (ncol(data) == 5)
        result = GRanges(data$chr, IRanges(data$start, data$end), id = data$id, score = data$score, seqinfo = seqinfo)
    else
        result = GRanges(data$chr, IRanges(data$start, data$end), id = data$id, score = data$score, strand = gsub("[^+-]+", "*", data$strand), seqinfo = seqinfo)

    result
}


readMaskRegions = function(prefix, seqinfo)
{
    rmsk0_mdust0 = bed2GRanges(paste(prefix, "rmsk0_mdust0.bed.gz", sep = ""), seqinfo)
    rmsk0_mdust1 = bed2GRanges(paste(prefix, "rmsk0_mdust1.bed.gz", sep = ""), seqinfo)
    rmsk1_mdust0 = bed2GRanges(paste(prefix, "rmsk1_mdust0.bed.gz", sep = ""), seqinfo)
    rmsk1_mdust1 = bed2GRanges(paste(prefix, "rmsk1_mdust1.bed.gz", sep = ""), seqinfo)

    list(rmsk0_mdust0 = rmsk0_mdust0, rmsk0_mdust1 = rmsk0_mdust1, rmsk1_mdust0 = rmsk1_mdust0, rmsk1_mdust1 = rmsk1_mdust1)
}


classifyZygosity = function(vcfs)
{
    lapply(vcfs, classifyZygosityVcf)
}


classifyZygosityVcf = function(vcf)
{
    genotypes = geno(vcf)$GT
    # Perform the following transformation:
    # Genotype                  Output
    # .                         NA
    # 0                         R (shouldn't occur in GIAB VCF)
    # 1                         A
    # 0/0 or 0|0                R/R (shouldn't occur in GIAB VCF)
    # 0/[^0] or 0|[^0]          R/A
    # [^0]/[^0] or [^0]|[^0]    A/A if both genotypes the same, A/B otherwise

    # Remove phasing information
    genotypes = gsub("\\|", "/", genotypes)

    # Split into pairs
    genotypes = strsplit(genotypes, "/", fixed = TRUE)

    # Convert components to numeric
    genotypes = suppressWarnings(lapply(genotypes, as.integer))

    # Perfom the zygosity classification
    result = sapply(genotypes, function(gt) {
        if (length(gt) == 1)
        {
            if (is.na(gt))      # as.numeric(".") = NA, so this condition catches the GT = "." case
                return("Unknown")
            
            if (gt == 0)
                return("R")
            return("A")
        }
        else if (length(gt) == 2)
        {
            if (is.na(gt[1]) || is.na(gt[2]))
                return("Unknown")
            else if (gt[1] == gt[2])
            {
                if (gt[1] == 0)
                    return("R/R")
                return("A/A")
            }
            else
            {
                if (gt[1] == 0 || gt[2] == 0)
                    return("R/A")
                return("A/B")
            }
        }
        else
        {
            warning(sprintf("Multiallelic genotype found (%s); I cannot yet handle this case.", gt), gt)
            return(NA)
        }
    })

    # For now, place haploid regions together with diploid.
    list(
        RR = Rle(result == "R" | result == "R/R"),
        AA = Rle(result == "A" | result == "A/A"),
        RA = Rle(result == "R/A"),
        AB = Rle(result == "A/B"),
        Unknown = Rle(result == "Unknown"))
}


classifyMutationType = function(vcfs)
{
    lapply(vcfs, classifyMutationTypeVcf)
}


classifyMutationTypeVcf = function(vcf)
{
    if (nrow(vcf) == 0)
    {
        return(list(Subst = Rle(), Ins = Rle(), Del = Rle(), Other = Rle(), None = Rle()))
    }
    result = list(
        Subst = Rle(isSubstitution(vcf)),
        Ins = Rle(isInsertion(vcf)),
        Del = Rle(isDeletion(vcf)))
    result$Other = !(result$Subst | result$Ins | result$Del)
    result$None = result$Subst & FALSE      # Always false at first; this is set by other logic (it should be manually set for FPs, and left alone otherwise)
    result
}


# Classifies the variants in vcf by their overlap with the GRanges 
# regions in the list regions.  Returns overlaps as a list of logical
# Rle objects, with each list member corresponding to a region class 
# in the regions list, and each item of each Rle vector corresponding
# to a variant in vcf.  The Rle elements are TRUE if the given variant
# overlaps any element in the given range, else FALSE.  The degree of
# overlap that is required for a TRUE overlap value is controlled by
# mode.  mode is a named vector of strings, with names corresponding
# to the region classes (names(mode) == names(regions)), and values of
# either "All" or "Any".  If "All", a variant is only called true if
# every one of its bases overlaps at least one range in regions; if
# "Any", a variant is called true if any base overlaps at least one
# range in regions.
classifyRegionOverlapVcf = function(vcf, regions, mode)
{
    stopifnot(all(sort(names(mode)) == sort(names(regions))))
    stopifnot(all(mode %in% c("Any", "All")))

    # This is required to ensure that region_mode = "All" works, as
    # it falsely fails if overlapping or abutting regions are present.
    regions = lapply(regions, reduce)

    sapply(names(regions), function(region_name) {
        this_region_granges = regions[[region_name]]
        this_region_mode = mode[region_name]
        if (this_region_mode == "Any")
            return(Rle(overlapsAny(rowRanges(vcf), regions[[region_name]], maxgap = 0L, minoverlap = 1L, type = "any", ignore.strand = TRUE)))
        else
            return(Rle(overlapsAny(rowRanges(vcf), regions[[region_name]], maxgap = 0L, minoverlap = 1L, type = "within", ignore.strand = TRUE)))
    }, simplify = FALSE, USE.NAMES = TRUE)
}


classifyRegionOverlap = function(vcfs, regions, mode)
{
    lapply(vcfs, classifyRegionOverlapVcf, regions = regions, mode = mode)
}


classifyMutationSize = function(vcfs)
{
    lapply(vcfs, classifyMutationSizeVcf)
}


classifyMutationSizeVcf = function(vcf)
{
    # breaks = c(seq(0, 50, 10), Inf)
    breaks = c(0, 5, seq(10, 50, 10), Inf)
    levels = levels(cut(0, breaks, right = FALSE))
    interval_start_inclusive = as.numeric(sub("^\\[", "", sub(",.*", "", levels)))
    interval_end_inclusive = as.numeric(sub(".*,", "", sub("\\)$", "", levels))) - 1
    interval_labels = mapply(function(start, end) { 
        if (start == end)
            return(start)
        else if (end == "Inf")
            return(paste(start, "+", sep = ""))
        else
            return(paste(start, "-", end, sep = "")) },
        interval_start_inclusive, interval_end_inclusive)

    if (nrow(vcf) == 0) 
    {
        result_list = lapply(interval_labels, function(x) Rle())
    }
    else
    {
        cut_lengths = cut(getMutationSizeVcf(vcf), breaks = breaks, labels = interval_labels, right = FALSE)
        result_list = lapply(interval_labels, function(l) Rle(cut_lengths == l))
    }

    names(result_list) = interval_labels
    result_list
}


getMutationSizeVcf = function(vcf)
{
    # Return the 'size' of the mutation.

    # For substitutions, this is the size of the alt
    # For insertions, the number of inserted bases (== width(alt) - width(ref))
    # For deletions, the number of deleted bases (== width(ref) - width(alt))
    # For delins, the largest of the number of inserted and deleted bases

    # In the case of multiple alt alleles, some variants will likely fall under multiple classes,
    # and therefore multiple size definitions.  In this situation, the largest passing size is
    # used (therefore all the pmaxes).

    subst = isSubstitution(vcf, singleAltOnly = FALSE)
    del = isDeletion(vcf, singleAltOnly = FALSE)
    ins = isInsertion(vcf, singleAltOnly = FALSE)
    delins = isDelins(vcf, singleAltOnly = FALSE)

    alts = alt(vcf)

    # The following code requires the fields alt_width_min and alt_width_max
    # to be present in the CollapsedVCF rowRanges.  These fields aren't present
    # in the original files, but are instead added at execute time by the 
    # function augmentCollapsedVCFWithAltLengthRange (see report_calculations.R).
    # This hack is to get around the canonical code to do this being extremely
    # slow.
    alt_widths_min = rowRanges(vcf)$alt_width_min
    alt_widths_max = rowRanges(vcf)$alt_width_max
    # The canonical (using CollapsedVCF / XStringSet calls) code which is 
    # replaced by the above, follows:
    # alt_widths = sapply(alts, function(alt_options) range(width(alt_options)))
    # alt_widths_min = alt_widths[1,]
    # alt_widths_max = alt_widths[2,]

    ref_widths = width(vcf)

    size = rep(0, nrow(vcf))

    size[subst] = ref_widths[subst]
    size[del] = pmax(size[del], ref_widths[del] - alt_widths_min[del])
    size[ins] = pmax(size[ins], alt_widths_max[ins] - ref_widths[ins])
    size[delins] = pmax(size[delins], pmax(alt_widths_max, ref_widths)[delins])
    size[size == 0] = NA
    size
}


classifyDepth = function(vcfs)
{
    lapply(vcfs, classifyDepthVcf)
}


classifyDepthVcf = function(vcf)
{
    breaks = c(seq(0, 50, 5), Inf)
    levels = levels(cut(0, breaks, right = FALSE))
    interval_start_inclusive = as.numeric(sub("^\\[", "", sub(",.*", "", levels)))
    interval_end_inclusive = as.numeric(sub(".*,", "", sub("\\)$", "", levels))) - 1
    interval_labels = mapply(function(start, end) { 
        if (start == end)
            return(start)
        else if (end == "Inf")
            return(paste(start, "+", sep = ""))
        else
            return(paste(start, "-", end, sep = "")) },
        interval_start_inclusive, interval_end_inclusive)
    interval_labels2 = c(interval_labels, "Unknown")

    if (nrow(vcf) == 0) 
    {
        result_list = lapply(interval_labels2, function(x) Rle())
    }
    else
    {
        depth = geno(vcf)$KCCG_PERF_DP_MIN
        cut_lengths = cut(depth, breaks = breaks, labels = interval_labels, right = FALSE)
        result_list = lapply(interval_labels2, function(l) { 
            if (l != "Unknown")
                return(Rle(!is.na(cut_lengths) & cut_lengths == l))
            else
                return(Rle(is.na(cut_lengths))) 
        } )
    }

    names(result_list) = interval_labels2
    result_list
}


tabulatePositiveCallTruth = cxxfunction(
    signature(scores_postruth = "numeric", scores_negtruth = "numeric", scores_uniq = "numeric"),
    body = '
        NumericVector _scores_postruth(scores_postruth);
        NumericVector _scores_negtruth(scores_negtruth);
        NumericVector _scores_uniq(scores_uniq);

        int n_uniq = _scores_uniq.length();
        int n_pos = _scores_postruth.length();
        int n_neg = _scores_negtruth.length();

        // Use floating point variables here, even though these tallies
        // are integers, to prevent overflow issues.
        NumericVector tally_postruth(n_uniq);
        NumericVector tally_negtruth(n_uniq);
        
        // The tally itself should still be an actual integer, to ensure that
        // increment operations always work.  Explicitly cast the completed
        // tallies to double, for assignation to the NumericVector objects.
        unsigned long tally;

        int i_uniq, i_pos, i_neg;
        i_pos = 0;
        i_neg = 0;

        for (i_uniq = 0; i_uniq < n_uniq; i_uniq++)
        {
            tally = 0;
            for (; i_pos < n_pos && _scores_postruth[i_pos] == _scores_uniq[i_uniq]; i_pos++)
                tally++;
            tally_postruth[i_uniq] = double(tally);

            tally = 0;
            for (; i_neg < n_neg && _scores_negtruth[i_neg] == _scores_uniq[i_uniq]; i_neg++)
                tally++;
            tally_negtruth[i_uniq] = double(tally);
        }

        List ret;
        ret["tally.postruth"] = tally_postruth;
        ret["tally.negtruth"] = tally_negtruth;

        return ret;
    ', plugin = "Rcpp"
)



vcfPerfWorker = function(scores.tp, scores.fp, n.tn.structural, n.fn.structural)
{
    scores.tp = sort(scores.tp)
    scores.fp = sort(scores.fp)

    # Currently does *not* handle non-finite scores.
    if (any(is.na(scores.tp)) || any(is.na(scores.fp)) || any(is.infinite(scores.tp)) || any(is.infinite(scores.fp)) || any(is.nan(scores.fp)) || any(is.nan(scores.fp)))
        stop("Infinite, NaN, or NA variant scores are not currently supported.")

    # Total counts of the tp and fp classes
    n.tp = length(scores.tp)
    n.fp = length(scores.fp)

    # Generate a tabulation of tp and fp counts as a function of 
    # a threshold.  Specifically, calculate:
    #   sum(scores.tp >= thresh)
    #   sum(scores.fp >= thresh)
    # for all possible thresholds.  By using increasing thresholds, 
    # a sweep of the count of numbers of TP and FP items above a swept
    # threshold results, and this is essentially the TP and FP counts
    # that are used in performance estimates.

    # First, generate the thresholds.  Select them so that every change
    # in score is captured.
    unique_vals = sort(unique(c(unique(scores.tp), unique(scores.fp))))
    unique_vals = unique_vals[is.finite(unique_vals)]
    if (length(unique_vals) == 1 && (unique_vals[1] == 1 || unique_vals[1] == 0))
        unique_vals = sort(unique(c(unique_vals, 0, 1)))
    thresholds = (unique_vals[-1] + unique_vals[-length(unique_vals)]) / 2

    # Test for an edge case, in which all scores are the same.
    # Return with a reasonable default, to enable plotting.  See below
    # for an explanation of the double Inf cutoffs.
    if (length(thresholds) == 0)
    {
        result = data.frame(
            #              "Everything is a variant."   "Everything non-NA in    "Everything's reference.  EVERYTHING."
            #                                           the VCF is a variant."        
            cutoff = c(    -Inf,                        -Inf,                     Inf,                     Inf),
            midpoint = c(  -Inf,                        -Inf,                     Inf,                     Inf),
            ntp = c(       n.tp + n.fn.structural,      n.tp,                     0,                       0),
            nfp = c(       n.fp + n.tn.structural,      n.fp,                     0,                       0),
            ntn = c(       0,                           n.tn.structural,          n.fp + n.tn.structural,  n.fp + n.tn.structural),
            nfn = c(       0,                           n.fn.structural,          n.tp + n.fn.structural,  n.tp + n.fn.structural))
        return(result)
    }

    # Generate the tallies.  This returns a list of vectors:
    #   tallies$tally.postruth == sapply(thresholds, function(thresh) sum(scores.tp == max(unique_vals[unique_vals < thresh])))
    #   tallies$tally.negtruth == sapply(thresholds, function(thresh) sum(scores.fp == max(unique_vals[unique_vals < thresh])))
    # Or more simply:
    #   tallies$tally.postruth == sapply(unique_vals, function(uniq) sum(scores.tp == uniq))
    #   tallies$tally.negtruth == sapply(unique_vals, function(uniq) sum(scores.fp == uniq))
    tallies = tabulatePositiveCallTruth(scores.tp, scores.fp, unique_vals)
    # The reason for the complexity of this approach is that the
    # tabulatePositiveCallTruth code is faster, and *far* more memory
    # efficient, than the equivalent sapply or table().

    # These tallies (effectively the counts of each unique_val in each
    # of scores.tp and scores.fp) are less useful than their cumulative
    # sums.  Calculate these sums, converting to float to avoid 
    # integer overflow.
    cumsum.tp = cumsum(as.numeric(tallies$tally.postruth))
    cumsum.fp = cumsum(as.numeric(tallies$tally.negtruth))

    # Use these cumulative sums to calculate the numbers of  
    # false negatives (fn), and true negatives (tn), then use these
    # to in turn derive the numbers of false positives (fp) and true
    # positives (tp).  Each quantity is calculated for a given 
    # threshold in thresholds.  The addition of two Inf values to each
    # end of thresholds is to accomodate two extreme cases:
    #  1) (inner Infs) a threshold so extreme, that every measured 
    #     variant falls beyond the threshold.  However, *structural 
    #     zeros* (as in n.fn.structural and n.tn.structural) are not included.
    #     For -Inf, this is equivalent to having no cutoffs at all: if
    #     a variant is present in the final VCF at all, it is 
    #     considered to be present.  However, variants *not* in the 
    #     VCF, are not considered -- such sites are supposed to be 
    #     reference.  
    #     For +Inf, this is equivalent to accepting absolutely no 
    #     variants, ignoring the genotyping output entirely (ie 
    #     calling all sites as reference, always).
    #  2) (outer Infs) thresholds so extreme that even structural 
    #     zeros are included.  For -Inf, this is equivalent to calling
    #     a variant at every single valid site in the genome, 
    #     regardless of whether a VCF entry is present for that
    #     variant or not (ie calling all sites as variant, always).  
    #     For +Inf, the result is that for the inner +Inf, but now 
    #     includes variants with NA scores.
    thresholds = c(-Inf, -Inf, thresholds, Inf, Inf)
    cutoffs = c(-Inf, unique_vals, Inf, Inf)

    total.neg.truth = n.fp + n.tn.structural
    total.pos.truth = n.tp + n.fn.structural

    path.n.tn = c(0, n.tn.structural, n.tn.structural + cumsum.fp, total.neg.truth)
    path.n.fn = c(0, n.fn.structural, n.fn.structural + cumsum.tp, total.pos.truth)

    path.n.fp = total.neg.truth - path.n.tn
    path.n.tp = total.pos.truth - path.n.fn

    result = data.frame(cutoff = cutoffs, midpoint = thresholds, ntp = path.n.tp, nfp = path.n.fp, ntn = path.n.tn, nfn = path.n.fn)

    result
}


vcfPerf = function(data, field_access_func)
{
    # Calculate the performance of the metric accessed by field_access_func,
    # at all observed levels, in the measurements data.  Essentially, extracts
    # a score from each item in data using field_access_func, then sorts these
    # scores, and uses the sorted order, and the true class of each item, to
    # derive performance metrics.
    #
    # The true classes are indicated by the fields of data:
    # Field         Truth       Description
    # data$vcf.tp   Variant     A Vcf object containing scores for true variants
    #                           that were present at all in the VCF.  Scores are
    #                           accessed using field_access_func(data$vcf.tp)
    # data$vcf.fp   Reference   A Vcf object containing scores for false 
    #                           variants that were present in any way in the VCF.
    #                           Scores are access using 
    #                           field_access_func(data$vcf.fp)
    # data$n.fn     Variant     The count of variants which are known to be 
    #                           present, but were not listed at all in the VCF.
    # data$n.tn     Reference   The count of variants which are known to NOT be
    #                           present, and were not listed in any way in the 
    #                           VCF.
    #
    # The fields, as well as the scores, define an ordering of the data.  From
    # lowest score to highest, this is:
    #   n.fn, n.tn, no internal ordering
    #   vcf.tp and vcf.fp -Inf scores, and vcf.tp and vcf.fp NA scores, no internal ordering
    #   vcf.fp and vcf.fp finite scores, internally ordered by increasing value
    #   vcf.tp and vcf.fp Inf scores
    # NOTE ABOVE: NA scores are at the lowest non-structural rank, equivalent
    # to a score of -Inf (in fact, in the code they are converted to this value).
    #
    # This function manipulates the data and counts to effectively simulate
    # that ordering, without needing to explicitly generate the very large
    # label vectors involved.
    #
    # Currently, the full ordering above is not implemented.  Particularly,
    # infinite, NaN, and NA scores are not supported at all.
    if (is.null(data))
    {
        return(data.frame(
            cutoff = c(    -Inf,  -Inf,   Inf,   Inf),
            midpoint = c(  -Inf,  -Inf,   Inf,   Inf),
            ntp = c(          0,     0,     0,     0),
            nfp = c(          0,     0,     0,     0),
            ntn = c(          0,     0,     0,     0),
            nfn = c(          0,     0,     0,     0)))
    }

    # First, extract variables
    # Raw scores
    scores.tp = field_access_func(data$vcf.tp)
    scores.fp = field_access_func(data$vcf.fp)

    if (is.null(scores.tp) || is.null(scores.fp))
    {
        # The required field is not present in the VCFs, or is all NA.
        stop("VCF does not support supplied field access function -- required field is probably missing")
    }

    # Despite the aspirational description above, currently does *not* handle non-finite scores.
    if (any(is.na(scores.tp)) || any(is.na(scores.fp)) || any(is.infinite(scores.tp)) || any(is.infinite(scores.fp)) || any(is.nan(scores.fp)) || any(is.nan(scores.fp)))
        stop("Infinite, NaN, or NA variant scores are not currently supported.")

    n.fn.structural = data$n.fn     # 'Structural' counts -- these are counts of variants that were
    n.tn.structural = data$n.tn     # not present in the VCF at all, and have the lowest rank.

    result = vcfPerfWorker(scores.tp, scores.fp, n.tn.structural, n.fn.structural)
    result
}


getPerfAtCutoff = function(perf, cutoff)
{
    perf = perf[order(perf$cutoff),]
    sel = min(which(perf$cutoff >= cutoff))
    unlist(perf[sel,])
}
