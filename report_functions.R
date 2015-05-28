library(ggplot2)
library(scales)



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
	data = read.table(path, header = FALSE, stringsAsFactors = FALSE)

	if (ncol(data) > 6)
		data = data[,1:6]
	else if (ncol(data) < 3)
		stop(sprintf("#s is not a valid BED file: fewer than three columns", path))

	header = c("chr", "start", "end", "id", "score", "strand")
	colnames(data) = header[1:ncol(data)]

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


perflist2data = function(perf_list)
{
	xvals = do.call(c, lapply(perf_list, function(perf) perf@x.values[[1]]))
	yvals = do.call(c, lapply(perf_list, function(perf) perf@y.values[[1]]))
	cvals = do.call(c, lapply(perf_list, function(perf) perf@alpha.values[[1]]))
	names = rep(names(perf_list), sapply(perf_list, function(perf) length(perf@x.values[[1]])))
	data.frame(name = names, x = xvals, y = yvals, cutoff = cvals)
}


# Some necessary messing about as ROCR silently shits itself with -Inf and Inf scores.
# Basically, finds the finite range of all scores, and defines extreme scores 
# (smallest_score and largest_score) as just beyond this range.  Then, sets all
# non-finite values to the appropriate extreme score, and also assigns these to
# the TN samples.  Returns a concatenated score vector, combined as:
#   TP, FN, FP, TN
# This ordering agrees with that of makeTruthVector.
makeScoreVector = function(vcf.tp, vcf.fp, n.fn, n.tn, field_access_func)
{
	scores.tp = field_access_func(vcf.tp)
	scores.fp = field_access_func(vcf.fp)
	# scores.fn not included, as those variants were not given scores at all by the genotyper

	score_range = range(c(range(scores.tp, na.rm = TRUE), range(scores.fp, na.rm = TRUE)), na.rm = TRUE)
	smallest_score = score_range[1] - 1
	largest_score = score_range[1] + 1

	scores.tp[scores.tp <= -Inf] = smallest_score
	scores.tp[scores.tp >= Inf] = largest_score
	scores.fp[scores.fp <= -Inf] = smallest_score
	scores.fp[scores.fp >= Inf] = largest_score
	scores.fn = rep(smallest_score, n.fn)
	scores.tn = rep(smallest_score, n.tn)

	c(scores.tp, scores.fn, scores.fp, scores.tn)
}


# Creates an ordered factor vector (values wt < mut) of true mutation status.
# Returns concatenated values in the roder TP, FN, FP, TN, which matches the
# ordering from makeScoreVector.
makeTruthVector = function(n.tp, n.fp, n.fn, n.tn)
{
	ordered(c(rep("mut", n.tp), rep("mut", n.fn), rep("wt", n.fp), rep("wt", n.tn)), levels = c("wt", "mut"))
}


plotROC = function(pred_list)
{
	perf_list = lapply(pred_list, function(pred) performance(pred, measure = "tpr", x.measure = "fpr"))
	data = perflist2data(perf_list)
	colnames(data) = c("Name", "FPR", "TPR", "Cutoff")

	ggplot(data, aes(x = FPR, y = TPR, colour = Name)) + geom_line() + 
		xlim(0, 1) + ylim(0, 1) + 
		coord_fixed() + geom_abline(intercept = 0, slope = 1, linetype = "dotted", alpha = 0.5) + 
		xlab("False positive rate") + ylab("True positive rate")
}


plotDET = function(pred_list)
{
	perf_list = lapply(pred_list, function(pred) performance(pred, measure = "fnr", x.measure = "fpr"))
	data = perflist2data(perf_list)
	colnames(data) = c("Name", "FPR", "FNR", "Cutoff")

#	probit_trans = function() trans_new(name = "probit", transform = function(x) qnorm(x*0.99+0.005), inverse = function(y) (pnorm(y)-0.005)/0.99, domain = c(0, 1))
	ggplot(data, aes(x = log10(FPR), y = log10(FNR), colour = Name)) + geom_line() + 
		xlab("log10(False positive rate)") + ylab("log10(False negative rate)")
}


plotLR = function(pred_list)
{
	perf_list = lapply(pred_list, function(pred) performance(pred, measure = "sens", x.measure = "spec"))
	data = perflist2data(perf_list)
	colnames(data) = c("Name", "spec", "sens", "Cutoff")
	data$LRP = data$sens / (1 - data$spec)
	data$LRN = (1 - data$sens) / data$spec

	ggplot(data, aes(x = LRN, y = LRP, colour = Name)) + geom_line() + 
		xlab("LR-") + ylab("LR+") + scale_x_log10() + scale_y_log10()
}

