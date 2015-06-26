# KCCG Validation Reporter

Create a WGS validation report, by comparing a test VCF of NA12878 calls, to a gold-standard reference from the GiaB consortium.

## Inputs
* [Required] Variants: A .vcf.gz of variants called in an NA12878 sample
* [Optional] Region BED: A .bed file of genomic intervals; only variants overlapping at least one of these intervals will be considered in calculations.  If missing, the entire genome is considered.
* [Optional] Extended report: A boolean indicating whether an extended report, with additional diagnostic plots, should be generated.  Default: false.

## Outputs
* A validation report PDF.