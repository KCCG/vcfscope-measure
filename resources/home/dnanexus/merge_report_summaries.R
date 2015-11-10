options(
    stringsAsFactors = FALSE, 
    warn = 2,                   # Treat warnings as errors -- this will be production code
    echo = TRUE)                # Emit lines to aid debugging


env = as.list(Sys.getenv())
param = list()
param$knitr.scratch = env$PARAM_KNITR_SCRATCH
param$sample.count = as.integer(env$LOOP_NUM_SAMPLES)
param$path.rds.output = env$PARAM_OUTPUT_RDS_PATH
param$path.json.output = env$PARAM_OUTPUT_JSON_PATH


merged_data = list()

for (sample.index in 0:(param$sample.count-1))
{
    sample.rds_path = paste(param$knitr.scratch, "_", sample.index, "/report_data.rds", sep = "")
    message(paste("  ", sample.rds_path, "...", sep = ""))
    sample.exported_summary = readRDS(sample.rds_path)
    sample.id = sample.exported_summary$params$sample.id
    if (sample.id %in% merged_data)
    {
        stop(sprintf("Error: Duplicate VCF sample IDs found (offending ID: %s).  Sample IDs must be unique.", sample.id))
    }
    merged_data[[sample.id]] = sample.exported_summary
}


if (param$path.rds.output != "")
    saveRDS(merged_data, param$path.rds.output, version = 2, compress = "xz")


if (param$path.json.output != "")
{
    library(jsonlite)
    merged_data_summary = lapply(merged_data, function(x) x$report_summary)
    cat(toJSON(merged_data_summary), file = param$path.json.output)
}
