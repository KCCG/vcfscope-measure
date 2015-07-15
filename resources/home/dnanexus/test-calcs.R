options(echo = TRUE)

library(testthat)
context("Performance calculations")


# Intermediate variables
load("report_data.rda")
library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome)
source("report_functions.R")


# Test the JSON serialization
test_that("JSON serialization is accurate", {
	if (param$path.json.output == "")
		skip("Script was not asked to output JSON file")

	library(jsonlite)
	export.fromjson = fromJSON(readChar(param$path.json.output, file.info(param$path.json.output)$size))
	expect_equal(export, export.fromjson)
})


# A convenience report function for debugging
perfAroundCutoff = function(perf, cutoff) {
	nearest_cutoff = which.min(abs(perf$cutoff - cutoff))
	start = max(1, nearest_cutoff - 5)
	end = min(nrow(perf), nearest_cutoff + 5)
	cbind(perf[start:end,], perf$cutoff[start:end] >= cutoff)
}


# Run tests on the continuous QUAL metric

TEST_METRIC = "QUAL"
TEST_CUTOFF = 200
TEST_METRIC_FUNC = function(x) rowRanges(x)$QUAL

snv.het.perf = calcSensSpecAtCutoff(perf.snv$zyg[[TEST_METRIC]]$RRvsRA, TEST_CUTOFF)
indelsubst.het.perf = calcSensSpecAtCutoff(perf.indelsubst$zyg[[TEST_METRIC]]$RRvsRA, TEST_CUTOFF)
snv.hom.perf = calcSensSpecAtCutoff(perf.snv$zyg[[TEST_METRIC]]$RRvsAA, TEST_CUTOFF)
indelsubst.hom.perf = calcSensSpecAtCutoff(perf.indelsubst$zyg[[TEST_METRIC]]$RRvsAA, TEST_CUTOFF)

snv.coding10.het.perf = calcSensSpecAtCutoff(perf.snv.coding10$zyg[[TEST_METRIC]]$RRvsRA, TEST_CUTOFF)
indelsubst.coding10.het.perf = calcSensSpecAtCutoff(perf.indelsubst.coding10$zyg[[TEST_METRIC]]$RRvsRA, TEST_CUTOFF)
snv.coding10.hom.perf = calcSensSpecAtCutoff(perf.snv.coding10$zyg[[TEST_METRIC]]$RRvsAA, TEST_CUTOFF)
indelsubst.coding10.hom.perf = calcSensSpecAtCutoff(perf.indelsubst.coding10$zyg[[TEST_METRIC]]$RRvsAA, TEST_CUTOFF)

snv.all.perf = calcSensSpecAtCutoff(vcfPerf(perfdata.snv, TEST_METRIC_FUNC), TEST_CUTOFF)
indelsubst.all.perf = calcSensSpecAtCutoff(vcfPerf(perfdata.indelsubst, TEST_METRIC_FUNC), TEST_CUTOFF)
snv.coding10.all.perf = calcSensSpecAtCutoff(vcfPerf(perfdata.snv.coding10, TEST_METRIC_FUNC), TEST_CUTOFF)
indelsubst.coding10.all.perf = calcSensSpecAtCutoff(vcfPerf(perfdata.indelsubst.coding10, TEST_METRIC_FUNC), TEST_CUTOFF)

all.perf = calcSensSpecAtCutoff(vcfPerf(list(vcf.tp = calls$tp, vcf.fp = calls$fp, n.tn = 0, n.fn = nrow(calls$fn)), TEST_METRIC_FUNC), TEST_CUTOFF)

test_that("vcfPerf matches ROCR; QUAL 200", {
	if (!require(ROCR))
		skip("ROCR not available")

	vcfPerfWorker.Test <- function(scores.tp, scores.fp, n.tn.structural, n.fn.structural)
	{
		scores.tp = sort(scores.tp)
		scores.fp = sort(scores.fp)
		scores = c(scores.tp, scores.fp, rep(-Inf, n.tn.structural), rep(-Inf, n.fn.structural))
		classes = c(rep(TRUE, length(scores.tp)), rep(FALSE, length(scores.fp)), rep(FALSE, n.tn.structural), rep(TRUE, n.fn.structural))

		prediction(scores, classes)
	}

	test_rocr = vcfPerfWorker.Test(TEST_METRIC_FUNC(perfdata.snv.coding10$vcf.tp), TEST_METRIC_FUNC(perfdata.snv.coding10$vcf.fp), perfdata.snv.coding10$n.tn, perfdata.snv.coding10$n.fn)
	test_rocr = data.frame(cutoff = unlist(test_rocr@cutoffs), tp = unlist(test_rocr@tp), fp = unlist(test_rocr@fp), tn = unlist(test_rocr@tn), fn = unlist(test_rocr@fn))
	test_rocr = test_rocr[rev(1:nrow(test_rocr)),]
	test_vcfperf = vcfPerf(perfdata.snv.coding10, TEST_METRIC_FUNC)
	test_vcfperf = test_vcfperf[-c(1, nrow(test_vcfperf)),]
	expect_equal(test_vcfperf$tp, test_rocr$tp)
	expect_equal(test_vcfperf$fp, test_rocr$fp)
	expect_equal(test_vcfperf$tn - test_vcfperf$tn[1], test_rocr$tn)
	expect_equal(test_vcfperf$fn - test_vcfperf$fn[1], test_rocr$fn)
})


# vcfPerfGrouped testing
test_that("VcfPerfGrouped shard sum; QUAL 200", {
	perf_full <- vcfPerf(perfdata.all, TEST_METRIC_FUNC)
	perf_sharded <- vcfPerfGrouped(perfdata.all, TEST_METRIC_FUNC, class.all.zyg)
	perf_full_atcutoff <- unlist(calcSensSpecAtCutoff(perf_full, TEST_CUTOFF))
	perf_sharded_atcutoff <- sapply(perf_sharded, function(x) unlist(calcSensSpecAtCutoff(x, TEST_CUTOFF)))
	perf_sharded_atcutoff_sum <- rowSums(perf_sharded_atcutoff)
	expect_equal(perf_full_atcutoff[c(3,6)], perf_sharded_atcutoff_sum[c(3,6)])
})


# Overall summary table test
test_that("WG tables; QUAL 200", {
	if (any(is.na(all.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls$tp) > TEST_CUTOFF), all.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls$fp) > TEST_CUTOFF), all.perf$nfp)
	expect_equal(nrow(calls$fn) + sum(TEST_METRIC_FUNC(calls$tp) <= TEST_CUTOFF), all.perf$nfn)
})


# Overall summary table tests, performance split by type
test_that("SNV WG tables; QUAL 200", {
	if (any(is.na(snv.all.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv$tp) > TEST_CUTOFF), snv.all.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv$fp) > TEST_CUTOFF), snv.all.perf$nfp)
	expect_equal(sum(as.numeric(width(reduce(calls.snv$tn)))) + sum(TEST_METRIC_FUNC(calls.snv$fp) <= TEST_CUTOFF), snv.all.perf$ntn)
	expect_equal(nrow(calls.snv$fn) + sum(TEST_METRIC_FUNC(calls.snv$tp) <= TEST_CUTOFF), snv.all.perf$nfn)
})

test_that("IndelSubst WG tables; QUAL 200", {
	if (any(is.na(indelsubst.all.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst$tp) > TEST_CUTOFF), indelsubst.all.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst$fp) > TEST_CUTOFF), indelsubst.all.perf$nfp)
	expect_equal(nrow(calls.indelsubst$fn) + sum(TEST_METRIC_FUNC(calls.indelsubst$tp) <= TEST_CUTOFF), indelsubst.all.perf$nfn)
})


# Overall summary table tests, performance split by type, zygosity
test_that("SNV WG Zyg tables; QUAL 200", {
	if (any(is.na(snv.het.perf)) || any(is.na(snv.hom.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv$tp[class.snv.zyg$RRvsRA$tp]) > TEST_CUTOFF), snv.het.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv$tp[class.snv.zyg$RRvsAA$tp]) > TEST_CUTOFF), snv.hom.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv$fp[class.snv.zyg$RRvsRA$fp]) > TEST_CUTOFF), snv.het.perf$nfp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv$fp[class.snv.zyg$RRvsAA$fp]) > TEST_CUTOFF), snv.hom.perf$nfp)
	expect_equal(sum(as.numeric(width(reduce(calls.snv$tn)))) + sum(TEST_METRIC_FUNC(calls.snv$fp[class.snv.zyg$RRvsRA$fp]) <= TEST_CUTOFF), snv.het.perf$ntn)
	expect_equal(sum(as.numeric(width(reduce(calls.snv$tn)))) + sum(TEST_METRIC_FUNC(calls.snv$fp[class.snv.zyg$RRvsRA$fp]) <= TEST_CUTOFF), snv.hom.perf$ntn)
	expect_equal(nrow(calls.snv$fn[class.snv.zyg$RRvsRA$fn]) + sum(TEST_METRIC_FUNC(calls.snv$tp[class.snv.zyg$RRvsRA$tp]) <= TEST_CUTOFF), snv.het.perf$nfn)
	expect_equal(nrow(calls.snv$fn[class.snv.zyg$RRvsAA$fn]) + sum(TEST_METRIC_FUNC(calls.snv$tp[class.snv.zyg$RRvsAA$tp]) <= TEST_CUTOFF), snv.hom.perf$nfn)
})


test_that("IndelSubst WG Zyg tables; QUAL 200", {
	if (any(is.na(indelsubst.het.perf)) || any(is.na(indelsubst.hom.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst$tp[class.indelsubst.zyg$RRvsRA$tp]) > TEST_CUTOFF), indelsubst.het.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst$tp[class.indelsubst.zyg$RRvsAA$tp]) > TEST_CUTOFF), indelsubst.hom.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst$fp[class.indelsubst.zyg$RRvsRA$fp]) > TEST_CUTOFF), indelsubst.het.perf$nfp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst$fp[class.indelsubst.zyg$RRvsAA$fp]) > TEST_CUTOFF), indelsubst.hom.perf$nfp)
	expect_equal(nrow(calls.indelsubst$fn[class.indelsubst.zyg$RRvsRA$fn]) + sum(TEST_METRIC_FUNC(calls.indelsubst$tp[class.indelsubst.zyg$RRvsRA$tp]) <= TEST_CUTOFF), indelsubst.het.perf$nfn)
	expect_equal(nrow(calls.indelsubst$fn[class.indelsubst.zyg$RRvsAA$fn]) + sum(TEST_METRIC_FUNC(calls.indelsubst$tp[class.indelsubst.zyg$RRvsAA$tp]) <= TEST_CUTOFF), indelsubst.hom.perf$nfn)
})


test_that("SNV C10 Zyg tables; QUAL 200", {
	if (any(is.na(snv.coding10.het.perf)) || any(is.na(snv.coding10.hom.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv.coding10$tp[class.snv.coding10.zyg$RRvsRA$tp]) > TEST_CUTOFF), snv.coding10.het.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv.coding10$tp[class.snv.coding10.zyg$RRvsAA$tp]) > TEST_CUTOFF), snv.coding10.hom.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv.coding10$fp[class.snv.coding10.zyg$RRvsRA$fp]) > TEST_CUTOFF), snv.coding10.het.perf$nfp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv.coding10$fp[class.snv.coding10.zyg$RRvsAA$fp]) > TEST_CUTOFF), snv.coding10.hom.perf$nfp)
	expect_equal(sum(as.numeric(width(reduce(calls.snv.coding10$tn)))) + sum(TEST_METRIC_FUNC(calls.snv.coding10$fp[class.snv.coding10.zyg$RRvsRA$fp]) <= TEST_CUTOFF), snv.coding10.het.perf$ntn)
	expect_equal(sum(as.numeric(width(reduce(calls.snv.coding10$tn)))) + sum(TEST_METRIC_FUNC(calls.snv.coding10$fp[class.snv.coding10.zyg$RRvsRA$fp]) <= TEST_CUTOFF), snv.coding10.hom.perf$ntn)
	expect_equal(nrow(calls.snv.coding10$fn[class.snv.coding10.zyg$RRvsRA$fn]) + sum(TEST_METRIC_FUNC(calls.snv.coding10$tp[class.snv.coding10.zyg$RRvsRA$tp]) <= TEST_CUTOFF), snv.coding10.het.perf$nfn)
	expect_equal(nrow(calls.snv.coding10$fn[class.snv.coding10.zyg$RRvsAA$fn]) + sum(TEST_METRIC_FUNC(calls.snv.coding10$tp[class.snv.coding10.zyg$RRvsAA$tp]) <= TEST_CUTOFF), snv.coding10.hom.perf$nfn)
})


test_that("IndelSubst C10 Zyg tables; QUAL 200", {
	if (any(is.na(indelsubst.coding10.het.perf)) || any(is.na(indelsubst.coding10.hom.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$tp[class.indelsubst.coding10.zyg$RRvsRA$tp]) > TEST_CUTOFF), indelsubst.coding10.het.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$tp[class.indelsubst.coding10.zyg$RRvsAA$tp]) > TEST_CUTOFF), indelsubst.coding10.hom.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$fp[class.indelsubst.coding10.zyg$RRvsRA$fp]) > TEST_CUTOFF), indelsubst.coding10.het.perf$nfp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$fp[class.indelsubst.coding10.zyg$RRvsAA$fp]) > TEST_CUTOFF), indelsubst.coding10.hom.perf$nfp)
	expect_equal(nrow(calls.indelsubst.coding10$fn[class.indelsubst.coding10.zyg$RRvsRA$fn]) + sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$tp[class.indelsubst.coding10.zyg$RRvsRA$tp]) <= TEST_CUTOFF), indelsubst.coding10.het.perf$nfn)
	expect_equal(nrow(calls.indelsubst.coding10$fn[class.indelsubst.coding10.zyg$RRvsAA$fn]) + sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$tp[class.indelsubst.coding10.zyg$RRvsAA$tp]) <= TEST_CUTOFF), indelsubst.coding10.hom.perf$nfn)
})



# Redo all tests, but now on the edge-case FILTER score

TEST_METRIC = "FILTER"
TEST_CUTOFF = 0.5
TEST_METRIC_FUNC = function(x) { (rowRanges(x)$FILTER == "PASS")*1 }


snv.het.perf = calcSensSpecAtCutoff(perf.snv$zyg[[TEST_METRIC]]$RRvsRA, TEST_CUTOFF)
indelsubst.het.perf = calcSensSpecAtCutoff(perf.indelsubst$zyg[[TEST_METRIC]]$RRvsRA, TEST_CUTOFF)
snv.hom.perf = calcSensSpecAtCutoff(perf.snv$zyg[[TEST_METRIC]]$RRvsAA, TEST_CUTOFF)
indelsubst.hom.perf = calcSensSpecAtCutoff(perf.indelsubst$zyg[[TEST_METRIC]]$RRvsAA, TEST_CUTOFF)

snv.coding10.het.perf = calcSensSpecAtCutoff(perf.snv.coding10$zyg[[TEST_METRIC]]$RRvsRA, TEST_CUTOFF)
indelsubst.coding10.het.perf = calcSensSpecAtCutoff(perf.indelsubst.coding10$zyg[[TEST_METRIC]]$RRvsRA, TEST_CUTOFF)
snv.coding10.hom.perf = calcSensSpecAtCutoff(perf.snv.coding10$zyg[[TEST_METRIC]]$RRvsAA, TEST_CUTOFF)
indelsubst.coding10.hom.perf = calcSensSpecAtCutoff(perf.indelsubst.coding10$zyg[[TEST_METRIC]]$RRvsAA, TEST_CUTOFF)

snv.all.perf = calcSensSpecAtCutoff(vcfPerf(perfdata.snv, TEST_METRIC_FUNC), TEST_CUTOFF)
indelsubst.all.perf = calcSensSpecAtCutoff(vcfPerf(perfdata.indelsubst, TEST_METRIC_FUNC), TEST_CUTOFF)
snv.coding10.all.perf = calcSensSpecAtCutoff(vcfPerf(perfdata.snv.coding10, TEST_METRIC_FUNC), TEST_CUTOFF)
indelsubst.coding10.all.perf = calcSensSpecAtCutoff(vcfPerf(perfdata.indelsubst.coding10, TEST_METRIC_FUNC), TEST_CUTOFF)

all.perf = calcSensSpecAtCutoff(vcfPerf(list(vcf.tp = calls$tp, vcf.fp = calls$fp, n.tn = 0, n.fn = nrow(calls$fn)), TEST_METRIC_FUNC), TEST_CUTOFF)


test_that("vcfPerf matches ROCR; FILTER PASS", {
	if (!require(ROCR))
		skip("ROCR not available")

	vcfPerfWorker.Test <- function(scores.tp, scores.fp, n.tn.structural, n.fn.structural)
	{
		scores.tp = sort(scores.tp)
		scores.fp = sort(scores.fp)
		scores = c(scores.tp, scores.fp, rep(-Inf, n.tn.structural), rep(-Inf, n.fn.structural))
		classes = c(rep(TRUE, length(scores.tp)), rep(FALSE, length(scores.fp)), rep(FALSE, n.tn.structural), rep(TRUE, n.fn.structural))

		prediction(scores, classes)
	}

	test_rocr = vcfPerfWorker.Test(TEST_METRIC_FUNC(perfdata.snv.coding10$vcf.tp), TEST_METRIC_FUNC(perfdata.snv.coding10$vcf.fp), perfdata.snv.coding10$n.tn, perfdata.snv.coding10$n.fn)
	test_rocr = data.frame(cutoff = unlist(test_rocr@cutoffs), tp = unlist(test_rocr@tp), fp = unlist(test_rocr@fp), tn = unlist(test_rocr@tn), fn = unlist(test_rocr@fn))
	test_rocr = test_rocr[rev(1:nrow(test_rocr)),]
	test_vcfperf = vcfPerf(perfdata.snv.coding10, TEST_METRIC_FUNC)
	test_vcfperf = test_vcfperf[-c(1, nrow(test_vcfperf)),]
	expect_equal(test_vcfperf$tp, test_rocr$tp)
	expect_equal(test_vcfperf$fp, test_rocr$fp)
	expect_equal(test_vcfperf$tn - test_vcfperf$tn[1], test_rocr$tn)
	expect_equal(test_vcfperf$fn - test_vcfperf$fn[1], test_rocr$fn)
})


test_that("VcfPerfGrouped shard sum; FILTER PASS", {
	perf_full <- vcfPerf(perfdata.all, TEST_METRIC_FUNC)
	perf_sharded <- vcfPerfGrouped(perfdata.all, TEST_METRIC_FUNC, class.all.zyg)
	perf_full_atcutoff <- unlist(calcSensSpecAtCutoff(perf_full, TEST_CUTOFF))
	perf_sharded_atcutoff <- sapply(perf_sharded, function(x) unlist(calcSensSpecAtCutoff(x, TEST_CUTOFF)))
	perf_sharded_atcutoff_sum <- rowSums(perf_sharded_atcutoff)
	expect_equal(perf_full_atcutoff[c(3,6)], perf_sharded_atcutoff_sum[c(3,6)])
})


# Overall summary table test
test_that("WG tables; FILTER PASS", {
	if (any(is.na(all.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls$tp) > TEST_CUTOFF), all.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls$fp) > TEST_CUTOFF), all.perf$nfp)
	expect_equal(nrow(calls$fn) + sum(TEST_METRIC_FUNC(calls$tp) <= TEST_CUTOFF), all.perf$nfn)
})


# Overall summary table tests, performance split by type
test_that("SNV WG tables; FILTER PASS", {
	if (any(is.na(snv.all.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv$tp) > TEST_CUTOFF), snv.all.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv$fp) > TEST_CUTOFF), snv.all.perf$nfp)
	expect_equal(sum(as.numeric(width(reduce(calls.snv$tn)))) + sum(TEST_METRIC_FUNC(calls.snv$fp) <= TEST_CUTOFF), snv.all.perf$ntn)
	expect_equal(nrow(calls.snv$fn) + sum(TEST_METRIC_FUNC(calls.snv$tp) <= TEST_CUTOFF), snv.all.perf$nfn)
})

test_that("IndelSubst WG tables; FILTER PASS", {
	if (any(is.na(indelsubst.all.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst$tp) > TEST_CUTOFF), indelsubst.all.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst$fp) > TEST_CUTOFF), indelsubst.all.perf$nfp)
	expect_equal(nrow(calls.indelsubst$fn) + sum(TEST_METRIC_FUNC(calls.indelsubst$tp) <= TEST_CUTOFF), indelsubst.all.perf$nfn)
})


# Overall summary table tests, performance split by type, zygosity
test_that("SNV WG Zyg tables; FILTER PASS", {
	if (any(is.na(snv.het.perf)) || any(is.na(snv.hom.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv$tp[class.snv.zyg$RRvsRA$tp]) > TEST_CUTOFF), snv.het.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv$tp[class.snv.zyg$RRvsAA$tp]) > TEST_CUTOFF), snv.hom.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv$fp[class.snv.zyg$RRvsRA$fp]) > TEST_CUTOFF), snv.het.perf$nfp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv$fp[class.snv.zyg$RRvsAA$fp]) > TEST_CUTOFF), snv.hom.perf$nfp)
	expect_equal(sum(as.numeric(width(reduce(calls.snv$tn)))) + sum(TEST_METRIC_FUNC(calls.snv$fp[class.snv.zyg$RRvsRA$fp]) <= TEST_CUTOFF), snv.het.perf$ntn)
	expect_equal(sum(as.numeric(width(reduce(calls.snv$tn)))) + sum(TEST_METRIC_FUNC(calls.snv$fp[class.snv.zyg$RRvsRA$fp]) <= TEST_CUTOFF), snv.hom.perf$ntn)
	expect_equal(nrow(calls.snv$fn[class.snv.zyg$RRvsRA$fn]) + sum(TEST_METRIC_FUNC(calls.snv$tp[class.snv.zyg$RRvsRA$tp]) <= TEST_CUTOFF), snv.het.perf$nfn)
	expect_equal(nrow(calls.snv$fn[class.snv.zyg$RRvsAA$fn]) + sum(TEST_METRIC_FUNC(calls.snv$tp[class.snv.zyg$RRvsAA$tp]) <= TEST_CUTOFF), snv.hom.perf$nfn)
})


test_that("IndelSubst WG Zyg tables; FILTER PASS", {
	if (any(is.na(indelsubst.het.perf)) || any(is.na(indelsubst.hom.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst$tp[class.indelsubst.zyg$RRvsRA$tp]) > TEST_CUTOFF), indelsubst.het.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst$tp[class.indelsubst.zyg$RRvsAA$tp]) > TEST_CUTOFF), indelsubst.hom.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst$fp[class.indelsubst.zyg$RRvsRA$fp]) > TEST_CUTOFF), indelsubst.het.perf$nfp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst$fp[class.indelsubst.zyg$RRvsAA$fp]) > TEST_CUTOFF), indelsubst.hom.perf$nfp)
	expect_equal(nrow(calls.indelsubst$fn[class.indelsubst.zyg$RRvsRA$fn]) + sum(TEST_METRIC_FUNC(calls.indelsubst$tp[class.indelsubst.zyg$RRvsRA$tp]) <= TEST_CUTOFF), indelsubst.het.perf$nfn)
	expect_equal(nrow(calls.indelsubst$fn[class.indelsubst.zyg$RRvsAA$fn]) + sum(TEST_METRIC_FUNC(calls.indelsubst$tp[class.indelsubst.zyg$RRvsAA$tp]) <= TEST_CUTOFF), indelsubst.hom.perf$nfn)
})


test_that("SNV C10 Zyg tables; FILTER PASS", {
	if (any(is.na(snv.coding10.het.perf)) || any(is.na(snv.coding10.hom.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv.coding10$tp[class.snv.coding10.zyg$RRvsRA$tp]) > TEST_CUTOFF), snv.coding10.het.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv.coding10$tp[class.snv.coding10.zyg$RRvsAA$tp]) > TEST_CUTOFF), snv.coding10.hom.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv.coding10$fp[class.snv.coding10.zyg$RRvsRA$fp]) > TEST_CUTOFF), snv.coding10.het.perf$nfp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.snv.coding10$fp[class.snv.coding10.zyg$RRvsAA$fp]) > TEST_CUTOFF), snv.coding10.hom.perf$nfp)
	expect_equal(sum(as.numeric(width(reduce(calls.snv.coding10$tn)))) + sum(TEST_METRIC_FUNC(calls.snv.coding10$fp[class.snv.coding10.zyg$RRvsRA$fp]) <= TEST_CUTOFF), snv.coding10.het.perf$ntn)
	expect_equal(sum(as.numeric(width(reduce(calls.snv.coding10$tn)))) + sum(TEST_METRIC_FUNC(calls.snv.coding10$fp[class.snv.coding10.zyg$RRvsRA$fp]) <= TEST_CUTOFF), snv.coding10.hom.perf$ntn)
	expect_equal(nrow(calls.snv.coding10$fn[class.snv.coding10.zyg$RRvsRA$fn]) + sum(TEST_METRIC_FUNC(calls.snv.coding10$tp[class.snv.coding10.zyg$RRvsRA$tp]) <= TEST_CUTOFF), snv.coding10.het.perf$nfn)
	expect_equal(nrow(calls.snv.coding10$fn[class.snv.coding10.zyg$RRvsAA$fn]) + sum(TEST_METRIC_FUNC(calls.snv.coding10$tp[class.snv.coding10.zyg$RRvsAA$tp]) <= TEST_CUTOFF), snv.coding10.hom.perf$nfn)
})


test_that("IndelSubst C10 Zyg tables; FILTER PASS", {
	if (any(is.na(indelsubst.coding10.het.perf)) || any(is.na(indelsubst.coding10.hom.perf)))
		skip("Insufficient samples for test")
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$tp[class.indelsubst.coding10.zyg$RRvsRA$tp]) > TEST_CUTOFF), indelsubst.coding10.het.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$tp[class.indelsubst.coding10.zyg$RRvsAA$tp]) > TEST_CUTOFF), indelsubst.coding10.hom.perf$ntp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$fp[class.indelsubst.coding10.zyg$RRvsRA$fp]) > TEST_CUTOFF), indelsubst.coding10.het.perf$nfp)
	expect_equal(sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$fp[class.indelsubst.coding10.zyg$RRvsAA$fp]) > TEST_CUTOFF), indelsubst.coding10.hom.perf$nfp)
	expect_equal(sum(as.numeric(width(reduce(calls.indelsubst.coding10$tn)))) + sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$fp[class.indelsubst.coding10.zyg$RRvsRA$fp]) <= TEST_CUTOFF), indelsubst.coding10.het.perf$ntn)
	expect_equal(sum(as.numeric(width(reduce(calls.indelsubst.coding10$tn)))) + sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$fp[class.indelsubst.coding10.zyg$RRvsRA$fp]) <= TEST_CUTOFF), indelsubst.coding10.hom.perf$ntn)
	expect_equal(nrow(calls.indelsubst.coding10$fn[class.indelsubst.coding10.zyg$RRvsRA$fn]) + sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$tp[class.indelsubst.coding10.zyg$RRvsRA$tp]) <= TEST_CUTOFF), indelsubst.coding10.het.perf$nfn)
	expect_equal(nrow(calls.indelsubst.coding10$fn[class.indelsubst.coding10.zyg$RRvsAA$fn]) + sum(TEST_METRIC_FUNC(calls.indelsubst.coding10$tp[class.indelsubst.coding10.zyg$RRvsAA$tp]) <= TEST_CUTOFF), indelsubst.coding10.hom.perf$nfn)
})
