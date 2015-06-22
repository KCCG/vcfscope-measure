#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail


main() {
  #
  # Fetch inputs (~/in/vcfgz/*)
  #
  dx-download-all-inputs --parallel

  #
  # Stream and unpack assets bundle
  #
  [[ -n "$DX_RESOURCES_ID" ]] && DX_ASSETS_ID="$DX_RESOURCES_ID" || DX_ASSETS_ID="$DX_PROJECT_CONTEXT_ID"
  mkdir ./resources
  cd ./resources
  dx cat "${DX_ASSETS_ID}:/assets/kccg_validation_reporter_resources_bundle-1.0.tar.gz" | tar zxf - # => functional_regions/* gold_standard/* kccg/* mask_regions/* orig/*
  cd -

  #
  # setup R
  #
  dx cat "${DX_ASSETS_ID}:/assets/R-3.2.0.compiled.packages_v2.tar.gz" | tar -zxf -
  export PATH="$PWD/bin:$PATH"
  export RHOME=${HOME} # This is needed to make RScript work, since it was compiled in a different dir.

  dx get "${DX_ASSETS_ID}:/assets/BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz"
  R CMD INSTALL BSgenome.HSapiens.1000g.37d5_1.0.0.tar.gz

  #
  # process options
  #
  if [ "extended" == "true" ]; then
    opts+=("-x")
  fi
  if [ -n "${region}" ]; then
    opts+=("-r" "${region}")
  fi

  #
  # run report
  #
  export PATH=".:$PATH"
  validation_report.sh -o report.pdf "${opts[@]}" ${vcfgz_path}

  #
  # upload results
  #
  mkdir -p ~/out/report/
  mv report.pdf ~/out/report/${outfile}
  dx-upload-all-outputs
  propagate-user-meta vcfgz report
}

#
# do this only once, during development.
# It facilitates both the Validation Report, and R plots
#
dev_update_aaron_R () {
  mkdir R-tmp; cd R-tmp
  dx cat "$DX_PROJECT_CONTEXT_ID:/assets/R-3.2.0.compiled.packages_v1.tar.gz" | tar -zxf -
  sed -i "s%/home/dnanexus%$PWD%" bin/R # Set R_HOME_DIR
  bin/R --quiet -e 'install.packages(c("inline", "RSQLite", "png", "gsalib"), repos="http://cran.rstudio.com")'
  # sudo apt-get install libxml2-dev
  bin/R --quiet -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("VariantAnnotation", "GenomicRanges", "BSgenome"))'
  bin/R --quiet -e 'paste(sort(rownames(installed.packages())), collapse=", ")'
  tar -pzcf ../R-3.2.0.compiled.packages_v2.tar.gz ./
  cd ..
  dx upload --path /assets/ R-3.2.0.compiled.packages_v2.tar.gz
  rm -rf R-tmp
}

build_R311 () {
sudo apt-get install libx11-dev libxml-dev openjdk-7-jre-headless
./configure --prefix=/home/dnanexus/R-tmp

}
