# KCCG Validation Reporter

Create a WGS validation report, by comparing a test VCF of NA12878 calls to a GiaB consortium gold-standard reference.

## Inputs
* [Required] Variants: A .vcf.gz of variants called in an NA12878 sample.  At this time, this should be a post-VQSR VCF; see below.
* [Optional] Region BED: A .bed file of genomic intervals; only variants overlapping at least one of these intervals will be considered in calculations.  If missing, the entire genome is considered.
* [Optional] Extended report: A boolean indicating whether an extended report, with additional diagnostic plots, should be generated.  Default: false.

## Outputs
* A validation report PDF.

# Running

## VCF requirements
The following fields are *required* to be present in the supplied VCF:

* QD
* MQ
* FS
* MQRankSum
* ReadPosRankSum
* VQSLOD
* QUAL
* FILTER (the pass-fail column)

The script does *not* currently verify the presence of the above fields.  If one is missing, the behaviour is undefined: the script may fail, or it may produce incorrect results for metrics involving the missing value.

### Multi-sample VCFs
These are currently handled, in which case a separate analysis is performed for each sample, and all sub-reports are concatenated to form the final validation report.  All samples should be NA12878, although this is not required or checked.

## Instance type
mem2_ssd1_x4 is recommended.
