#!/bin/bash -e

dx select "VCFscope resources"

dx run vcfscope-measure/9.9.7 \
    -ivcfgz=/test/input/GKX_30.realigned.recalibrated.hc.vqsr.chr21.vcf.gz \
    -ibam=/test/input/GKX_05.realigned.chr21.aroundvariants.bam \
    -ibai=/test/input/GKX_05.realigned.chr21.aroundvariants.bam.bai \
    -iregion=/test/input/reportable_range.chr21.bed \
    --folder /test/output --yes --watch
# SUCCESS job-F1JpY9Q0pVj0BgpYBp14f31Q
