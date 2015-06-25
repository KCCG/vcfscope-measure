# KCCG Validation Reporter

Create a WGS validation report from a VCF of NA12878 calls.



## Usage

```
Usage: validation_report.sh [-o OUTFILE] [-r BEDFILE] [-x] <INFILE>

Create a WGS validation report from a VCF of NA12878 calls.

    INFILE       Input NA12878 genotype calls, in vcf.gz format.
    -o OUTFILE   Write the report to OUTFILE (default: report.pdf)
    -r BEDFILE   Restrict analysis to the regions in BEDFILE only.
                 Default: the full genome is considered.
    -x           Generate an extended report, with threshold and
                 score diagnostics appended to the standard report.
                 Default: generate the standard report only.
    -h           Display this help and exit.

v1.1.0

Mark Pinese
20150625
```



### Examples

```
./validation_report.sh /directflow/ClinicalGenomicsPipeline/projects/validation-reporter/test_data/HiSeqX_v2_TKCC/R_150203_DAVMIL1_FGS_M001.hc.vqsr.vep.vcf.gz
```



## Requirements

* R, with the following packages:
 * VariantAnnotation (Bioconductor)
 * GenomicRanges (Bioconductor)
 * BSgenome (Bioconductor)
 * Rcpp (CRAN)
 * inline (CRAN)
 * ggplot2 (CRAN)
 * knitr (CRAN)
 * xtable (CRAN)
 * BSgenome.HSapiens.1000g.37d5 (this is a custom package, for installation instructions see below)
* A functional LaTeX installation
* Real Time Genomics Core (from [here](http://realtimegenomics.com/products/rtg-core-downloads/), currently included in the resources bundle)
* The following resources (see below for how to access a ready-made resource bundle)
 * NIST GIAB NA12878 VCF and valid region BED
 * KCCG hardmask callable region BED (adapted from the mask described [here](https://ccg.garvan.org.au/confluence/display/TxGen/Depth+and+quality+requirements+for+clinical+sequencing))



### BSgenome.HSapiens.1000g.37d5 package

This package is available at `/share/ClusterShare/biodata/contrib/marpin/reference/hs37d5/build/BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz`

To install, run the following:

```
cd /share/ClusterShare/biodata/contrib/marpin/reference/hs37d5/build/
R CMD INSTALL BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz
```

If you encounter a permissions error, it is likely because you are trying to install to the system-wide R libraries location, and not your user-specific one.  Google what to do about that and try again.



### Resources bundle

A bundle of the required resources is currently available at `/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources`.  This location is specified by the value of `$RESOURCES_HEAD` in `validation_report.sh`, which is already set by default to that bundle location.  So, provided the script will have access to directflow, the existing resource bundle will be accessed by default.  However, if the bundle 
location is changed, the value of `$RESOURCES_HEAD` variable will need to be changed to match.



## Code points and modifying

The main script is `validation_report.sh`.  This first runs the RTG variant evaluator, then calls `report_calculations.R` to perform the main calculations, and finally calls R with `report.Rnw` to generate the report.

To modify the default locations of resources and software, edit `validation_report.sh`.  Adding sections to the report will need changes to `report.Rnw`, and possibly also `report_calculations.R` or `report_functions.R` if further computation is required.


# Blame

Currently very much hacked together by Mark Pinese.
