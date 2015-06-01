# KCCG Validation Reporter

A script to build a validation report from a VCF of NA12878 calls.



## Usage

```
./validation_report.sh [-d CHROM] [-o OUTFILE] <INFILE>

    INFILE       Input NA12878 genotype calls, in vcf.gz format
    -o OUTFILE   Write the report to OUTFILE (default: report.pdf)
    -h           Display this help and exit
    -d CHROM     Debug mode.  Currently, adds additional debug 
                 information to the report, and examines chromosome
                 CHROM only.
```



### Example

```bash
./validation_report.sh /directflow/ClinicalGenomicsPipeline/projects/validation-reporter/test_data/HiSeqX_v2_TKCC/R_150203_DAVMIL1_FGS_M001.hc.vqsr.vep.vcf.gz
```

Running with the -d flag will dramatically speed things up (by only considering chromosome 22 in calculations), but will also pollute the output pdf with debug messages.



## Requirements

* R, with the following packages:
  - VariantAnnotation (Bioconductor)
  - GenomicRanges (Bioconductor)
  - BSgenome (Bioconductor)
  - Rcpp (CRAN)
  - inline (CRAN)
  - ggplot2 (CRAN)
  - knitr (CRAN)
  - BSgenome.HSapiens.1000g.37d5 (this is a custom package, for installation instructions see below)
* A functional LaTeX installation
* Real Time Genomics Core (from [here](http://realtimegenomics.com/products/rtg-core-downloads/), currently included in the resources bundle)
* The following resources (see below for how to access a ready-made resource bundle)
  - NIST GIAB NA12878 VCF and valid region BED
  - KCCG hardmask callable region BED (adapted from the mask described [here](https://ccg.garvan.org.au/confluence/display/TxGen/Depth+and+quality+requirements+for+clinical+sequencing))



### BSgenome.HSapiens.1000g.37d5 package

This package is available at `/share/ClusterShare/biodata/contrib/marpin/reference/hs37d5/build/BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz`

To install, run the following:

```bash
cd /share/ClusterShare/biodata/contrib/marpin/reference/hs37d5/build/
R CMD INSTALL BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz
```

If you encounter a permissions error, it is likely because you are trying to install to the system-wide R libraries location, and not your user-specific one.  Google what to do about that and try again.



### Resources bundle

A bundle of the required resources is currently available at `/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources`.  This location is specified by the value of `$RESOURCES_HEAD` in `validation_report.sh`, which is already set by default to that bundle location.  So, provided the script will have access to directflow, the existing resource bundle will be accessed by default.  However, if the bundle 
location is changed, the value of `$RESOURCES_HEAD` variable will need to be changed to match.


## Limitations

Currently many.  The report is far from complete.  Both the report and the R calculations consider only SNVs at this time.  There is no stratification of performance by variant type or location.  There is no pass/fail logic on performance metrics.  There is much yet to do.


## Code points and modifying

The main script is `validation_report.sh`.  This first runs the RTG variant evaluator, then calls `report_calculations.R` to perform the main calculations, and finally calls R with `report.Rnw` to generate the report.

To modify the default locations of resources and software, edit `validation_report.sh`.  Adding sections to the report will require changes to `report.Rnw`, and possibly also `report_calculations.R` or `report_functions.R` if further computation is required.


# Blame

Currently very much hacked together by Mark Pinese
