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


tabulatePositiveCallTruth = cxxfunction(
	signature(scores_postruth = "numeric", scores_negtruth = "numeric", scores_uniq = "numeric"),
	body = '
		NumericVector _scores_postruth(scores_postruth);
		NumericVector _scores_negtruth(scores_negtruth);
		NumericVector _scores_uniq(scores_uniq);

		int n_uniq = _scores_uniq.length();
		int n_pos = _scores_postruth.length();
		int n_neg = _scores_negtruth.length();

		IntegerVector tally_postruth(n_uniq);
		IntegerVector tally_negtruth(n_uniq);

		int i_uniq, i_pos, i_neg, tally;
		i_pos = 0;
		i_neg = 0;

		for (i_uniq = 0; i_uniq < n_uniq; i_uniq++)
		{
			tally = 0;
			for (; i_pos < n_pos && _scores_postruth[i_pos] == _scores_uniq[i_uniq]; i_pos++)
				tally++;
			tally_postruth[i_uniq] = tally;

			tally = 0;
			for (; i_neg < n_neg && _scores_negtruth[i_neg] == _scores_uniq[i_uniq]; i_neg++)
				tally++;
			tally_negtruth[i_uniq] = tally;
		}

		List ret;
		ret["tally.postruth"] = tally_postruth;
		ret["tally.negtruth"] = tally_negtruth;

		return ret;
	', plugin = "Rcpp")


# An optimized function to calculate TPR, FPR, TNR, and FNR directly from
# VCF objects, without needing to convert them to score and truth vectors.
# Also includes quite a bit of logic to deal with odd cases mostly unique 
# to genome comparisons.
vcfPerf = function(vcf.tp, vcf.fp, n.fn.always, n.tn.always, field_access_func)
{
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
		total = result$tp + result$fp + result$tn + result$fn
		result$tp = result$tp / total
		result$fp = result$fp / total
		result$tn = result$tn / total
		result$fn = result$fn / total
		return(result)
	}

	tallies = tabulatePositiveCallTruth(scores.tp, scores.fp, unique_vals)

	# Use the cumulative sums of the tallies to calculate numbers of 
	# false negatives (fn), and true negatives (tn), then use these
	# to in turn derive the numbers of false positives (fp) and true
	# positives (tp).  Each quantity is calculated for a given threshold
	# in thresholds.  The addition of two Inf values to each end of
	# thresholds is to accomodate two extreme cases:
	#  1) (inner Infs) a threshold so extreme, that every measured variant
	#     falls beyond the threshold.  However, *structural zeros* 
	#     (as in n.fn.always and n.tn.always) are not included.
	#     For -Inf, this is equivalent to having no cutoffs at all: if
	#     a variant is present in the final VCF at all, it is considered
	#     to be present.  However, variants *not* in the VCF, are not
	#     considered -- such sites are supposed to be reference.  
	#     For +Inf, this is equivalent to accepting absolutely no 
	#     variants, ignoring the genotyping output entirely (ie calling 
	#     all sites as reference, always).
	#  2) (outer Infs) thresholds so extreme that even structural zeros
	#     are included.  For -Inf, this is equivalent to calling a variant
	#     at every single valid site in the genome, regardless of whether
	#     a VCF entry is present for that variant or not (ie calling all
	#     sites as variant, always).  For +Inf, the result is that for the
	#     inner +Inf, but now includes variants with NA scores.
	thresholds = c(-Inf, -Inf, thresholds, Inf, Inf)
	path.n.fn = c(0, c(0, cumsum(as.numeric(tallies$tally.pos))) + as.numeric(n.fn.always), n.truth.pos + missing.tp)
	path.n.tn = c(0, c(0, cumsum(as.numeric(tallies$tally.neg))) + as.numeric(n.tn.always), n.truth.neg + missing.fp)
	path.n.fp = n.truth.neg - path.n.tn
	path.n.tp = n.truth.pos - path.n.fn

	result = data.frame(cutoff = thresholds, tp = path.n.tp, fp = path.n.fp, tn = path.n.tn, fn = path.n.fn)
	total = result$tp + result$fp + result$tn + result$fn
	result$tp = result$tp / total
	result$fp = result$fp / total
	result$tn = result$tn / total
	result$fn = result$fn / total

	result
}


plotROC = function(perf_list)
{
	perf2_list = lapply(perf_list, function(perf) data.frame(TPR = perf$tp / (perf$tp + perf$fn), FPR = perf$fp / (perf$fp + perf$tn), cutoff = perf$cutoff))
	data = data.frame(Name = rep(names(perf2_list), sapply(perf2_list, nrow)), TPR = unlist(sapply(perf2_list, function(x) x$TPR)), FPR = unlist(sapply(perf2_list, function(x) x$FPR)), Cutoff = unlist(sapply(perf2_list, function(x) x$cutoff)))

	ggplot(data, aes(x = FPR, y = TPR, colour = Name)) + geom_path() + 
		xlim(0, 1) + ylim(0, 1) + 
		coord_fixed() + geom_abline(intercept = 0, slope = 1, linetype = "dotted", alpha = 0.5) + 
		xlab("False positive rate") + ylab("True positive rate")
}


plotDET = function(perf_list)
{
	perf2_list = lapply(perf_list, function(perf) data.frame(FNR = perf$fn / (perf$tp + perf$fn), FPR = perf$fp / (perf$fp + perf$tn), cutoff = perf$cutoff))
	data = data.frame(Name = rep(names(perf2_list), sapply(perf2_list, nrow)), FPR = unlist(sapply(perf2_list, function(x) x$FPR)), FNR = unlist(sapply(perf2_list, function(x) x$FNR)), Cutoff = unlist(sapply(perf2_list, function(x) x$cutoff)))

	ggplot(data, aes(x = log10(FPR), y = log10(FNR), colour = Name)) + geom_path() + 
		xlab("log10(False positive rate)") + ylab("log10(False negative rate)")
}


plotLR = function(perf_list)
{
	perf2_list = lapply(perf_list, function(perf) data.frame(sens = perf$tp / (perf$tp + perf$fn), spec = perf$tn / (perf$fp + perf$tn), cutoff = perf$cutoff))
	data = data.frame(Name = rep(names(perf2_list), sapply(perf2_list, nrow)), sens = unlist(sapply(perf2_list, function(x) x$sens)), spec = unlist(sapply(perf2_list, function(x) x$spec)), Cutoff = unlist(sapply(perf2_list, function(x) x$cutoff)))
	data$LRP = data$sens / (1 - data$spec)
	data$LRN = (1 - data$sens) / data$spec

	ggplot(data, aes(x = LRN, y = LRP, colour = Name)) + geom_path() + 
		xlab("LR-") + ylab("LR+") + scale_x_log10() + scale_y_log10()
}
