# KCCG Performance Measure

Calculate variant caller performance, by comparing a test VCF of NA12878 calls to a GiaB consortium gold-standard reference.  Produces a raw data file that is then processed by KCCG Performance Report to make a human-readable report.

## Inputs
* [Required] Variants: A `.vcf.gz` of variants called in an NA12878 sample.
* [Required] Reads: A `.bam` of mapped reads used to generate the variants above.
* [Required] Read index: A `.bai` index for the BAM file above.
* [Optional] Region BED: A `.bed` file of genomic intervals; only variants overlapping at least one of these intervals will be considered in calculations.  If missing, the entire genome and reportable range are considered.

## Output
* RDS of performance metrics.
