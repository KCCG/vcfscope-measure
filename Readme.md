# KCCG Performance Reporter

Create a WGS performance report, by comparing a test VCF of NA12878 calls to a GiaB consortium gold-standard reference.

## Inputs
* [Required] Variants: A `.vcf.gz` of variants called in an NA12878 sample.
* [Optional] Region BED: A `.bed` file of genomic intervals; only variants overlapping at least one of these intervals will be considered in calculations.  If missing, the entire genome and reportable range are considered.
* [Optional] VCF Sample IDs: A string containing a comma-separated list of the sample IDs in the VCF for which to generate the report.  If '*', the report will be based on all samples in the VCF.
* [Optional] Extended report: A boolean indicating whether an extended report, with additional diagnostic plots, should be generated.  Default: false.

## Outputs
* Performance report PDF.
* JSON of summary performance statistics.
* RDS of performance metrics.


# Running

### Multi-sample VCFs
A separate analysis is performed for each sample matching the supplied VCF Sample IDs, and all sub-reports are concatenated to form the final performance report.  All samples should be NA12878, although this is not required or checked.

## Instance type
`mem2_ssd1_x4` is recommended.
