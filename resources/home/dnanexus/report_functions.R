library(ggplot2)
library(scales)
library(inline)
library(Rcpp)


texquote = function(str) gsub("_", "\\\\_", sub("\\s+$", "", str))


readpipe = function(command, n = 1e6) {
	conn = pipe(command)
	data = readChar(conn, n)
	close(conn)
	data
}


fileMD5 = function(path) {
	readpipe(sprintf("md5sum %s | sed 's/ .*//'", path))
}


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
	ambiguous = 		bed2GRanges(paste(prefix, "ambiguous.bed.gz", sep = ""), seqinfo)
	low_complexity = 	bed2GRanges(paste(prefix, "mdust.bed.gz", sep = ""), seqinfo)
	repetitive = 		bed2GRanges(paste(prefix, "repetitive.bed.gz", sep = ""), seqinfo)

	list(ambiguous = ambiguous, low_complexity = low_complexity, repetitive = repetitive)
}


readFunctionalRegions = function(prefix, seqinfo)
{
	coding = 		bed2GRanges(paste(prefix, "exonic_coding.bed.gz", sep = ""), seqinfo)
	intronic = 		bed2GRanges(paste(prefix, "intronic.bed.gz", sep = ""), seqinfo)
	utr = 			bed2GRanges(paste(prefix, "exonic_utr.bed.gz", sep = ""), seqinfo)
	intergenic = 	bed2GRanges(paste(prefix, "intergenic.bed.gz", sep = ""), seqinfo)
	genic = 		bed2GRanges(paste(prefix, "genes.bed.gz", sep = ""), seqinfo)
	splice = 		bed2GRanges(paste(prefix, "splice.bed.gz", sep = ""), seqinfo)

	list(coding = coding, intronic = intronic, utr = utr, intergenic = intergenic, genic = genic, splice = splice)
}


# Transform l, a list of lists, with outer list A and inner list
# B, into a list of lists, with outer list B and inner list A.
# For example, transforms:
#   list(A1 = list(B1 = "w", B2 = "x"), A2 = list(B1 = "y", B2 = "z", B3 = "p"))
# into
#   list(B1 = list(A1 = "w", A2 = "y"), B2 = list(A1 = "x", A2 = "z"), B3 = list(A1 = NULL, A2 = "p"))
# Returns the transformed list.  Does not preserve item ordering
# within levels.
#
# R> str(swapListOrder(list(A1 = list(B1 = "w", B2 = "x"), A2 = list(B1 = "y", B2 = "z", B3 = "p"))))
# List of 3
#  $ B1:List of 2
#   ..$ A1: chr "w"
#   ..$ A2: chr "y"
#  $ B2:List of 2
#   ..$ A1: chr "x"
#   ..$ A2: chr "z"
#  $ B3:List of 2
#   ..$ A1: NULL
#   ..$ A2: chr "p"
swapListOrder = function(l)
{
	outer_names = suppressWarnings(sort(names(l)))
	inner_names = sort(unique(unlist(lapply(l, names))))
	sapply(inner_names, function(inner_name) {
		sapply(outer_names, function(outer_name) {
			l[[outer_name]][[inner_name]]
		}, simplify = FALSE, USE.NAMES = TRUE)
	}, simplify = FALSE, USE.NAMES = TRUE)
}


classifyZygosity = function(vcfs)
{
	swapListOrder(lapply(vcfs, classifyZygosityVcf))
}


classifyZygosityVcf = function(vcf)
{
	genotypes = geno(vcf)$GT
	# Perform the following transformation:
	# Genotype 					Output
	# . 						NA
	# 0 						R (shouldn't occur in GIAB VCF)
	# 1 						A
	# 0/0 or 0|0 				R/R (shouldn't occur in GIAB VCF)
	# 0/[^0] or 0|[^0] 			R/A
	# [^0]/[^0] or [^0]|[^0]	A/A if both genotypes the same, A/B otherwise

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
			if (is.na(gt))		# as.numeric(".") = NA, so this condition catches the GT = "." case
				return(NA)
			
			if (gt == 0)
				return("R")
			return("A")
		}
		else if (length(gt) == 2)
		{
			if (gt[1] == gt[2])
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
		RRvsAA = Rle(result == "R" | result == "R/R" | result == "A" | result == "A/A"),
		RRvsRA = Rle(result == "R" | result == "R/R" | result == "R/A"),
		RRvsAB = Rle(result == "R" | result == "R/R" | result == "A/B"))
}


classifySomy = function(vcfs)
{
	swapListOrder(lapply(vcfs, classifySomyVcf))
}


classifySomyVcf = function(vcf)
{
	list(
		autosomal = seqnames(vcf) %in% as.character(1:22),
		X = seqnames(vcf) == "X",
		Y = seqnames(vcf) == "Y",
		MT = seqnames(vcf) %in% c("MT", "M"))
}


classifyMutationType = function(vcfs)
{
	swapListOrder(lapply(vcfs, classifyMutationTypeVcf))
}


classifyMutationTypeVcf = function(vcf)
{
	result = list(
		SNV = Rle(isSNV(vcf)),
		InsDelSubst = Rle(isInsertion(vcf) | isDeletion(vcf) | isSubstitution(vcf)))
	result$Other = !(result$SNV | result$InsDelSubst)
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
	swapListOrder(lapply(vcfs, classifyRegionOverlapVcf, regions = regions, mode = mode))
}



classifyMutationSize = function(vcfs)
{
	swapListOrder(lapply(vcfs, classifyMutationSizeVcf))
}


classifyMutationSizeVcf = function(vcf)
{
	size = getMutationSizeVcf(vcf)
	list(
		"[01,02)" = Rle(size == 1 & !is.na(size)),
		"[02,06)" = Rle(size >= 2 & size < 6 & !is.na(size)),
		"[06,11)" = Rle(size >= 6 & size < 11 & !is.na(size)),
		"[11,21)" = Rle(size >= 11 & size < 21 & !is.na(size)),
		"[21,inf)" = Rle(size >= 21 & !is.na(size)),
		"NA" = Rle(is.na(size)))
}


getMutationSizeVcf = function(vcf)
{
	# Return the 'size' of the mutation.
	# For SNVs, this is always 1
	# For insertions, the number of inserted bases (in the most likely
	# genotype, if multiple)
	# For deletions, the number of deleted bases
	# For other, NA.
	size = rep(NA, length(vcf))
	size[isSNV(vcf)] = 1
	size[isDeletion(vcf)] = width(vcf[isDeletion(vcf)])
	temp = alt(vcf[isInsertion(vcf)])
	size[isInsertion(vcf)] = nchar(unlist(temp)[start(PartitioningByEnd(temp))])
	size
}



# TODO: Update this doc to match the new class storage format.
# Subset the variant classes in class to just those variants for which
# subset is TRUE, and return as a list.  class must be a list 
# containing members "tp", "fp", and "fn"; each of these in turn is a 
# list of Rle objects, encoding whether or not each variant (in one of
# the classes TP, FP, or FN) is in that given list item's class.  
# subset is a list that also contains named members "tp", "fp", and
# "fn", each of which is also an Rle object encoding whether or not 
# the variant is in the single subset of interest.  The result is a 
# list of lists of logical Rle objects.  The outer list contains 
# elements named "tp", "fp", "tn", and "fn", and each in turn contains 
# corresponding elements from class, subset as per the elements in 
# subset.  The "tn" element of the result is constructed from the tn
# argument: if it is supplied, then result$tn is set to this value; if
# it is missing, then result$tn is set to zero for all classes.
subsetClass = function(class, subset, tn = NULL)
{
	class_group_names = names(class)

	result = sapply(class_group_names, function(class_group_name) {
		sapply(names(class[[class_group_name]]), function(call_type_name) {
			stopifnot(length(class[[class_group_name]][[call_type_name]]) == length(subset[[call_type_name]]))
			class[[class_group_name]][[call_type_name]][subset[[call_type_name]]]
		}, simplify = FALSE, USE.NAMES = TRUE)
	}, simplify = FALSE, USE.NAMES = TRUE)

	for (i in names(result))
	{
		if (!("tn" %in% names(result[[i]])))
		{
			if (is.null(tn))
				result[[i]]$tn = Rle()
			else
				result[[i]]$tn = tn[[i]]
		}
	}

	result
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
	', plugin = "Rcpp")


# An optimized function to calculate nTP, nFP, nTN, and nFN directly 
# from Vcf objects, without needing to convert them to score and truth
# vectors.  Also includes quite a bit of logic to deal with odd cases 
# mostly unique to genome comparisons.
vcfPerf = function(data, field_access_func)
{
	vcf.tp = data$vcf.tp
	vcf.fp = data$vcf.fp
	n.fn.always = data$n.fn
	n.tn.always = data$n.tn

	scores.tp = sort(field_access_func(vcf.tp))
	scores.fp = sort(field_access_func(vcf.fp))

	missing.tp = sum(as.numeric(is.na(scores.tp)))
	missing.fp = sum(as.numeric(is.na(scores.fp)))

	scores.tp = scores.tp[!is.na(scores.tp)]
	scores.fp = scores.fp[!is.na(scores.fp)]

	n.truth.pos = as.numeric(length(scores.tp)) + as.numeric(n.fn.always)
	n.truth.neg = as.numeric(length(scores.fp)) + as.numeric(n.tn.always)

	unique_vals = sort(unique(c(unique(scores.tp), unique(scores.fp))))
	thresholds = (unique_vals[-1] + unique_vals[-length(unique_vals)]) / 2

	# All scores are the same -- there is not much selection going on!
	# Return with a reasonable default, to enable plotting.  See below
	# for an explanation of the double Inf cutoffs.
	if (length(thresholds) == 0)
	{
		result = data.frame(
			#              "If it's in the              "If it has any      "Everything's reference.  EVERYTHING."
			#               VCF at all,                  valid score at
			#               call it as real."            all, call it."
			cutoff = c(    -Inf,                        -Inf,               Inf,                         Inf),
			tp = c(        n.truth.pos + missing.tp,    n.truth.pos,        0,                           0),
			fp = c(        n.truth.neg + missing.fp,    n.truth.neg,        0,                           0),
			tn = c(        0,                           0,                  n.truth.neg + missing.fp,    n.truth.neg + missing.fp),
			fn = c(        0,                           0,                  n.truth.pos + missing.tp,    n.truth.pos + missing.tp))
		return(result)
	}

	tallies = tabulatePositiveCallTruth(scores.tp, scores.fp, unique_vals)

	# Use the cumulative sums of the tallies to calculate numbers of 
	# false negatives (fn), and true negatives (tn), then use these
	# to in turn derive the numbers of false positives (fp) and true
	# positives (tp).  Each quantity is calculated for a given 
	# threshold in thresholds.  The addition of two Inf values to each
	# end of thresholds is to accomodate two extreme cases:
	#  1) (inner Infs) a threshold so extreme, that every measured 
	#     variant falls beyond the threshold.  However, *structural 
	#     zeros* (as in n.fn.always and n.tn.always) are not included.
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
	path.n.fn = c(0, c(0, cumsum(as.numeric(tallies$tally.pos))) + as.numeric(n.fn.always), n.truth.pos + missing.tp)
	path.n.tn = c(0, c(0, cumsum(as.numeric(tallies$tally.neg))) + as.numeric(n.tn.always), n.truth.neg + missing.fp)
	path.n.fp = n.truth.neg - path.n.tn
	path.n.tp = n.truth.pos - path.n.fn

	result = data.frame(cutoff = thresholds, tp = path.n.tp, fp = path.n.fp, tn = path.n.tn, fn = path.n.fn)

	result
}


# Run vcfPerf on calls, after subsetting it into groups defined
# by subgroups.  subgroups is a list of lists, with each member
# either being an Rle-logical vector of subgroup membership 
# indicators, or a single numeric value containing the total count
# of sites in the respective category.
vcfPerfGrouped = function(calls, field_access_func, subgroups)
{
	sapply(names(subgroups), function(subgroup_name) {
		this_group_data = list(
			vcf.tp = calls$vcf.tp[subgroups[[subgroup_name]]$tp],
			vcf.fp = calls$vcf.fp[subgroups[[subgroup_name]]$fp],
			n.fn = sum(subgroups[[subgroup_name]]$fn),
			n.tn = sum(subgroups[[subgroup_name]]$tn))
		vcfPerf(this_group_data, field_access_func)
	}, simplify = FALSE, USE.NAMES = TRUE)
}


# Helper function to calculate the ROC axes (either TPR, TPN, FPR, or 
# FPN), for plotROC.
calcROCaxes = function(perf, type.fp = c("rate", "count"), type.tp = c("rate", "count"))
{
	type.tp = match.arg(type.tp)
	type.fp = match.arg(type.fp)

	if (type.tp == "rate")
	{
		if (type.fp == "rate")
			perf2 = data.frame(TP = perf$tp / (perf$tp + perf$fn), FP = perf$fp / (perf$fp + perf$tn), cutoff = perf$cutoff)
		else
			perf2 = data.frame(TP = perf$tp / (perf$tp + perf$fn), FP = perf$fp, cutoff = perf$cutoff)
	}
	else
	{
		if (type.fp == "rate")
			perf2 = data.frame(TP = perf$tp, FP = perf$fp / (perf$fp + perf$tn), cutoff = perf$cutoff)
		else
			perf2 = data.frame(TP = perf$tp, FP = perf$fp, cutoff = perf$cutoff)
	}

	perf2
}


# Convert the performance data in perf to a data frame, for use with 
# ggplot2.  Parameters are as per plotROC.  Returns a list, with two 
# items:
#   data: the data frame for plotting
#   case: a character denoting the type of the data.  Either "A", "B", 
#         or "C"; definitions of these values are in the plotROC 
#         comments.
createROCData = function(perf, type.fp = c("rate", "count"), type.tp = c("rate", "count"))
{
	type.tp = match.arg(type.tp)
	type.fp = match.arg(type.fp)

	if (class(perf) == "data.frame")
	{
		# Case A
		case = "A"
		plotData = calcROCaxes(perf, type.fp, type.tp)
		plotData$group1 = NA
		plotData$group2 = NA
	}
	else if (class(perf) == "list" && length(perf) > 0 && all(sapply(perf, class) == "data.frame"))
	{
		# Case B
		case = "B"
		temp = lapply(perf, calcROCaxes, type.fp = type.fp, type.tp = type.tp)
		plotData = do.call(rbind, temp)
		plotData$group1 = factor(rep(names(temp), sapply(temp, nrow)))
		plotData$group2 = NA
	}
	else if (class(perf) == "list" && length(perf) > 0 && all(sapply(perf, class) == "list") && min(sapply(perf, length)) > 0 && all(sapply(perf, function(p1) all(sapply(p1, class) == "data.frame"))))
	{
		# Case C
		case = "C"
		temp = lapply(perf, function(perf_sub) { 
			temp = lapply(perf_sub, calcROCaxes, type.fp = type.fp, type.tp = type.tp)
			plotData_sub = do.call(rbind, temp)
			plotData_sub$group2 = rep(names(temp), sapply(temp, nrow))
			plotData_sub 
		})
		plotData = do.call(rbind, temp)
		plotData$group1 = factor(rep(names(temp), sapply(temp, nrow)))
		plotData$group2 = factor(plotData$group2)
	}
	else
	{
		stop("Error: perf is not of a supported type.  Must be either data.frame, list(data.frame), or list(list(data.frame))")
	}

	list(data = plotData, case = case)
}

# Create a ggplot2 ROC plot of the performance data in perf.
# perf may be one of three types:
#   A) a data frame
#   B) a list of data frames
#   C) a list of lists of data frames
# where each data frame contains the columns "cutoff", "tp", "fp", 
# "tn", and "fn".
#
# In case A, a plot is generated with just one curve; that described 
# by the single data frame.  In case B, one ROC curve is drawn per 
# list item, with each curve in a different colour.  In case C, each
# combination of upper-level and lower-level list is used to create
# a ROC curve, with the upper level determining curve colour, and the
# lower level curve line type.
#
# ROCs may use either total number, or rate, of false positives; this 
# is controlled independently type.fp parameter.
plotROC = function(perf, type.fp = c("rate", "count"), facet = c())
{
	type.fp = match.arg(type.fp)

	ROCData = createROCData(perf, type.fp)

	case = ROCData$case
	plotData = ROCData$data

	if (case == "A")
		plot = ggplot(plotData, aes(x = FP, y = TP))
	else if (case == "B")
	{
		if (1 %in% facet)
			plot = ggplot(plotData, aes(x = FP, y = TP)) + facet_wrap(~ group1, scales = ifelse(type.fp == "rate", "fixed", "free_x"))
		else
			plot = ggplot(plotData, aes(x = FP, y = TP, colour = group1))
	}
	else if (case == "C")
	{
		if (1 %in% facet)
		{
			if (2 %in% facet)
				plot = ggplot(plotData, aes(x = FP, y = TP)) + facet_grid(group1 ~ group2, scales = ifelse(type.fp == "rate", "fixed", "free_x"))
			else
				plot = ggplot(plotData, aes(x = FP, y = TP, linetype = group2)) + facet_wrap(~ group1, scales = ifelse(type.fp == "rate", "fixed", "free_x"))
		}
		else
		{
			if (2 %in% facet)
				plot = ggplot(plotData, aes(x = FP, y = TP, colour = group1)) + facet_wrap(~ group2, scales = ifelse(type.fp == "rate", "fixed", "free_x"))
			else
				plot = ggplot(plotData, aes(x = FP, y = TP, colour = group1, linetype = group2))
		}
	}
	else
		stop("Assertion Error: execution should not have reached this point.  Check plotROC case logic.")

 	plot = plot + geom_path()

	plot = plot + ylim(0, 1) + ylab("True positive rate")

	if (type.fp == "rate")
		plot = plot + xlim(0, 1) + xlab("False positive rate") + coord_fixed() + geom_abline(intercept = 0, slope = 1, linetype = "dotted", alpha = 0.5)
	else
		plot = plot + xlab("False positive count")

	plot
}


calcSensSpecAtCutoff = function(perf, cutoff)
{
	sel = which.max(perf$cutoff * (perf$cutoff <= cutoff))
	list(sens = perf$tp[sel] / (perf$tp[sel] + perf$fn[sel]), spec = perf$tn[sel] / (perf$tn[sel] + perf$fp[sel]), ntp = perf$tp[sel], nfp = perf$fp[sel], ntn = perf$tn[sel], nfn = perf$fn[sel])
}
