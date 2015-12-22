#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail


main() {
  #
  # Fetch inputs
  #
  dx-download-all-inputs --parallel

  # Move the BAI to the location of the BAM
  mv ${bai_path} ${bam_path}.bai
  touch ${bam_path}.bai

  #
  # Locate the assets bundle, the location of which varies, depending on whether
  # this is an app or an applet.
  #
  if [[ "$DX_RESOURCES_ID" != "" ]]; then
    # This is an app; fetch assets from the app's private asset container
    DX_ASSETS_ID="$DX_RESOURCES_ID"
  else
    # This is an applet; fetch assets from the parent project
    DX_ASSETS_ID="$DX_PROJECT_CONTEXT_ID"
  fi

  #
  # Stream and unpack assets bundle
  #
  mkdir ~/resources
  cd ~/resources
  dx cat "${DX_ASSETS_ID}:/assets/kccg_performance_reporter_resources_bundle-2.0.tar" | tar xf -

  #
  # setup R
  #
  cd ~
  # This R is Aaron's R-3.2.0 with pre-installed packages, with the following 
  # additional packages pre-installed:
  # CRAN: inline, RSQLite, png, gsalib
  # BioC: VariantAnnotation, GenomicRanges, BSgenome
  dx cat "${DX_ASSETS_ID}:/assets/R-3.2.0.compiled.packages_v2.tar.gz" | tar -zxf -
  export PATH="$PWD/bin:$PATH"
  export RHOME=${HOME} # This is needed to make RScript work, since it was compiled in a different dir.

  dx get "${DX_ASSETS_ID}:/assets/BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz"
  R CMD INSTALL BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz

  dx get "${DX_ASSETS_ID}:/assets/testthat_0.10.0.tar.gz"
  dx get "${DX_ASSETS_ID}:/assets/memoise_0.2.1.tar.gz"
  dx get "${DX_ASSETS_ID}:/assets/crayon_1.3.0.tar.gz"
  dx get "${DX_ASSETS_ID}:/assets/digest_0.6.8.tar.gz"
  R CMD INSTALL digest_0.6.8.tar.gz
  R CMD INSTALL memoise_0.2.1.tar.gz
  R CMD INSTALL crayon_1.3.0.tar.gz
  R CMD INSTALL testthat_0.10.0.tar.gz

  # An updated VariantAnnotation package is required to handle HAS2.0 VCFs
  dx get "${DX_ASSETS_ID}:/assets/VariantAnnotation_1.14.6.tar.gz"
  R CMD INSTALL VariantAnnotation_1.14.6.tar.gz

  #
  # process options
  #
  opts=()
  if [ -n "${region}" ]; then
    opts+=("-r" "${region_path}")
  fi

  #
  # run report
  #
  mkdir -p ~/out/report/ ~/out/rds
  sample_basename=$(basename ${vcfgz_path} .vcf.gz)
  ./performance_measure.sh "${opts[@]}" "${vcfgz_path}" "${bam_path}" "/home/dnanexus/out/rds/${sample_basename}.perfmeas.rds"

  #
  # upload results
  #
  dx-upload-all-outputs
  propagate-user-meta vcfgz rds
}
