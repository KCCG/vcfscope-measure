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


classifyZygosity = function(vcf)
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
	genotypes = lapply(genotypes, as.numeric)

	# Perfom the zygosity classification
	result = sapply(genotypes, function(gt) {
		if (length(gt) == 1)
		{
			if (gt == 0)
			{
				warning(sprintf("Hemizygous reference genotype found (%s); this should not occur in a GIAB VCF.", gt), gt)
				return("R")
			}
			else if (is.na(gt))		# as.numeric(".") = NA, so this condition catches the GT = "." case
				return(NA)
			return("A")
		}
		else if (length(gt) == 2)
		{
			if (gt[1] == gt[2])
			{
				if (gt[1] == 0)
				{
					warning(sprintf("Homozygous reference genotype found (%s); this should not occur in a GIAB VCF.", gt), gt)
					return("R/R")
				}
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

	data.frame(
		R = result == "R",
		A = result == "A",
		RR = result == "R/R",
		RA = result == "R/A",
		AA = result == "A/A",
		AB = result == "A/B")
}


simplifyZygosityClass = function(class)
{
	data.frame(
		R = class$R | class$RR,
		A = class$A | class$AA,
		RA = class$RA,
		AB = class$AB)
}


classifySomy = function(vcf)
{
	result = as.data.frame(matrix(NA, nrow(vcf), 4))
	colnames(result) = c("autosomal", "X", "Y", "MT")
	result$autosomal = as.vector(seqnames(vcf) %in% as.character(1:22))
	result$X = as.vector(seqnames(vcf) == "X")
	result$Y = as.vector(seqnames(vcf) == "Y")
	result$MT = as.vector(seqnames(vcf) %in% c("MT", "M"))
	result
}


classifyMutationType = function(vcf)
{
	result = as.data.frame(matrix(NA, nrow(vcf), 4))
	colnames(result) = c("SNV", "Ins", "Del", "Other")
	result$SNV = isSNV(vcf)
	result$Ins = isInsertion(vcf)
	result$Del = isDeletion(vcf)
	result$Other = !(result$SNV | result$Ins | result$Del)
	result$Indel = result$Ins | result$Del
	result
}


# Classifies the variants in vcf by their overlap with
# the GRanges regions in the list regions.  Returns overlaps
# as a factor matrix, with rows equal to the number of
# elements in vcf, and columns equal to the number of
# elements in regions.  Each cell will contain one of the
# following values:
#   0L	The feature in this row does not overlap
#		any intervals in this column's GRanges, by
#		any amount.
#	1L	The feature in this row partially overlaps
#		an interval of this column's GRanges.
#	2L	The feature in this row is completely 
#		contained in an interval of this column's 
#		GRanges.
classifyRegion = function(vcf, regions)
{
	result = matrix(NA, nrow = length(rowData(vcf)), ncol = length(regions))
	colnames(result) = names(regions)
	rownames(result) = names(rowData(vcf))

	regions = lapply(regions, reduce)

	for (region_name in names(regions))
	{
		overlap_full = overlapsAny(rowData(vcf), regions[[region_name]], maxgap = 0L, minoverlap = 1L, type = "within", ignore.strand = TRUE)
		overlap_any = overlapsAny(rowData(vcf), regions[[region_name]], maxgap = 0L, minoverlap = 1L, type = "any", ignore.strand = TRUE)
		overlap_partial = overlap_any & !overlap_full
		overlap_none = !overlap_any
		result[overlap_none, region_name] = 0L
		result[overlap_partial, region_name] = 1L
		result[overlap_full, region_name] = 2L
	}

	result
}


simplifyRegionClass = function(region_classes, min_overlap_levels)
{
	t(t(region_classes) >= min_overlap_levels[colnames(region_classes)])
}


getMutationSize = function(vcf)
{
	# Return the 'size' of the mutation.
	# For SNVs, this is always 1
	# For insertions, the number of inserted bases (in the most likely genotype, if multiple)
	# For deletions, the number of deleted bases
	# For other, NA.
	size = rep(NA, length(vcf))
	size[isSNV(vcf)] = 1
	size[isDeletion(vcf)] = width(vcf[isDeletion(vcf)])
	temp = alt(vcf[isInsertion(vcf)])
	size[isInsertion(vcf)] = nchar(unlist(temp)[start(PartitioningByEnd(temp))])
	size
}


classifyMutationSize = function(size)
{
	result = as.data.frame(matrix(NA, nrow = length(size), ncol = 6))
	colnames(result) = c("[1,2)", "[2,6)", "[6,11)", "[11,21)", "[21,inf)", "NA")
	result[,"[1,2)"] = size == 1 & !is.na(size)
	result[,"[2,6)"] = size >= 2 & size < 6 & !is.na(size)
	result[,"[6,11)"] = size >= 6 & size < 11 & !is.na(size)
	result[,"[11,21)"] = size >= 11 & size < 21 & !is.na(size)
	result[,"[21,inf)"] = size >= 21 & !is.na(size)
	result[,"NA"] = is.na(size)
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


# An optimized function to calculate nTP, nFP, nTN, and nFN directly from
# VCF objects, without needing to convert them to score and truth vectors.
# Also includes quite a bit of logic to deal with odd cases mostly unique 
# to genome comparisons.
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

	result
}


vcfPerfGrouped = function(data, field_access_func, subgroups)
{
	group_names = unique(unlist(lapply(subgroups, colnames)))
	result = lapply(group_names, function(subgroup_name) {
		this_group_data = list(
			vcf.tp = data$vcf.tp[subgroups$tp[,subgroup_name]],
			vcf.fp = data$vcf.fp[subgroups$fp[,subgroup_name]],
			n.fn = sum(subgroups$fn[,subgroup_name]),
			n.tn = 0)
		vcfPerf(this_group_data, field_access_func)
	})
	names(result) = group_names
	result
}


plotROC = function(perf_list, type.fpr = c("rate", "count"), type.tpr = c("rate", "count"))
{
	type.tpr = match.arg(type.tpr)
	type.fpr = match.arg(type.fpr)

	if (type.tpr == "rate")
	{
		if (type.fpr == "rate")
			perf2_list = lapply(perf_list, function(perf) data.frame(TP = perf$tp / (perf$tp + perf$fn), FP = perf$fp / (perf$fp + perf$tn), cutoff = perf$cutoff))
		else
			perf2_list = lapply(perf_list, function(perf) data.frame(TP = perf$tp / (perf$tp + perf$fn), FP = perf$fp, cutoff = perf$cutoff))
	}
	else
	{
		if (type.fpr == "rate")
			perf2_list = lapply(perf_list, function(perf) data.frame(TP = perf$tp, FP = perf$fp / (perf$fp + perf$tn), cutoff = perf$cutoff))
		else
			perf2_list = lapply(perf_list, function(perf) data.frame(TP = perf$tp, FP = perf$fp, cutoff = perf$cutoff))
	}

	data = data.frame(Name = rep(names(perf2_list), sapply(perf2_list, nrow)), TP = unlist(sapply(perf2_list, function(x) x$TP)), FP = unlist(sapply(perf2_list, function(x) x$FP)), Cutoff = unlist(sapply(perf2_list, function(x) x$cutoff)))

	plot = ggplot(data, aes(x = FP, y = TP, colour = Name)) + geom_path()

	if (type.tpr == "rate")
		plot = plot + ylim(0, 1) + ylab("True positive rate")
	else
		plot = plot + ylab("True positive count")

	if (type.fpr == "rate")
		plot = plot + xlim(0, 1) + xlab("False positive rate")
	else
		plot = plot + xlab("False positive count")
	
	if (type.tpr == "rate" && type.fpr == "rate")
		plot = plot + coord_fixed() + geom_abline(intercept = 0, slope = 1, linetype = "dotted", alpha = 0.5)

	plot
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


plotTPRFNR = function(perf_list)
{
	perf2_list = lapply(perf_list, function(perf) data.frame(TPR = perf$tp / (perf$tp + perf$fn), FPR = perf$fp / (perf$fp + perf$tn), cutoff = perf$cutoff))
	data = data.frame(
		Name = rep(names(perf2_list), sapply(perf2_list, nrow)), 
		Value = c(unlist(sapply(perf2_list, function(x) x$TPR)), unlist(sapply(perf2_list, function(x) x$FPR))), 
		Type = rep(c("TPR", "FNR"), each = sum(sapply(perf2_list, nrow))),
		Cutoff = unlist(sapply(perf2_list, function(x) x$cutoff)))

	ggplot(data, aes(x = Cutoff, y = Value, colour = Name)) + geom_path() + 
		ylim(0, 1) + 
		xlab("Cutoff") + ylab("Performance value") + 
		facet_wrap( ~ Type, scales = "free_x")
}
