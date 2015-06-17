# KCCG Validation Reporter

A script to build a validation report from a VCF of NA12878 calls.



## Usage

```
./validation_report.sh [-d CHROM] [-t TMPDIR] [-r] [-f] [-x] [-o OUTFILE] <INFILE>

Create a WGS validation report.

    INFILE       Input NA12878 genotype calls, in vcf.gz format.
    -o OUTFILE   Write the report to OUTFILE (default: report.pdf)
    -d CHROM     Debug mode.  Currently, adds additional debug 
                 information to the report, and examines chromosome
                 CHROM only.
    -t TMPDIR    Path to a temporary scratch space directory.  If 
                 missing, mktemp will be used to create one.
    -f           Force noninteractive mode; all Y/N prompts will be
                 automatically answered Y.  Default is to run in 
                 interactive mode.
    -r           Resume flag.  If set, and TMPDIR points to a 
                 partially-completed run, will attempt to resume
                 this run.  Default is to not resume; if a partial
                 run is found in TMPDIR and -f is not set, query
                 the user as to what to do.  With -f set, any partial
                 run data in TMPDIR is automatically cleared, and
                 the validation report generation started from the 
                 beginning.  Regardless of the resume setting, the
                 final report generation is always repeated with
                 every invocation.
    -x           Generate an extended report, with threshold and
                 score diagnostics.
    -h           Display this help and exit.

v20150615-1

Mark Pinese
```



### Examples

Interactive, debug restricted to chromosome 22:
```bash
./validation_report.sh -d 22 /directflow/ClinicalGenomicsPipeline/projects/validation-reporter/test_data/HiSeqX_v2_TKCC/R_150203_DAVMIL1_FGS_M001.hc.vqsr.vep.vcf.gz
```

Non-interactive, full run:
```bash
qsub -pe smp 4 -cwd -b y -j y bash ./validation_report.sh -f /directflow/ClinicalGenomicsPipeline/projects/validation-reporter/test_data/HiSeqX_v2_TKCC/R_150203_DAVMIL1_FGS_M001.hc.vqsr.vep.vcf.gz
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
 * BSgenome.HSapiens.1000g.37d5 (this is a custom package, for installation instructions see below)
* A functional LaTeX installation
* Real Time Genomics Core (from [here](http://realtimegenomics.com/products/rtg-core-downloads/), currently included in the resources bundle)
* The following resources (see below for how to access a ready-made resource bundle)
 * NIST GIAB NA12878 VCF and valid region BED
 * KCCG hardmask callable region BED (adapted from the mask described [here](https://ccg.garvan.org.au/confluence/display/TxGen/Depth+and+quality+requirements+for+clinical+sequencing))



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



## Code points and modifying

The main script is `validation_report.sh`.  This first runs the RTG variant evaluator, then calls `report_calculations.R` to perform the main calculations, and finally calls R with `report.Rnw` to generate the report.

To modify the default locations of resources and software, edit `validation_report.sh`.  Adding sections to the report will need changes to `report.Rnw`, and possibly also `report_calculations.R` or `report_functions.R` if further computation is required.


# Blame

Currently very much hacked together by Mark Pinese.
