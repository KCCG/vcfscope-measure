options(echo = TRUE)

library(testthat)
context("Performance calculations")

source("report_functions.R")

temp = readRDS("report_data.rds")
params = temp$params
class_subsets.performance_thresholded = temp$class_subsets.performance_thresholded
sampleid = temp$sampleid


