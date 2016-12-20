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
  dx cat "${DX_ASSETS_ID}:/assets/vcfscope_reporter_resources_bundle-2.0.tar" | tar xf -

  #
  # process options
  #
  opts=()
  if [ -n "${region}" ]; then
    opts+=("-r" "${region_path}")
  fi

  if [ "$store_variants" == "true" ]; then
    opts+=("-a")
  fi

  #
  # run vcfscope_measure.sh to create an rds (ie R data file)
  #
  mkdir -p ~/out/rds
  sample_basename=$(basename ${vcfgz_path} .vcf.gz)
  ./vcfscope_measure.sh "${opts[@]}" "${vcfgz_path}" "${bam_path}" "/home/dnanexus/out/rds/${sample_basename}.vcfscope.rds"

  #
  # upload results
  #
  dx-upload-all-outputs
  propagate-user-meta vcfgz rds
}
