#!/bin/bash -e

dx select kccg
dx mkdir -p /test/vcfscope-measure/input /test/vcfscope-measure/output

dx run vcfscope-measure/9.9.7 \
    -ivcfgz=/test/kccg-performance-measure/input/GKX_30.realigned.recalibrated.hc.vqsr.chr21.vcf.gz \
    -ibam=/test/kccg-performance-measure/input/GKX_05.realigned.chr21.aroundvariants.bam \
    -ibai=/test/kccg-performance-measure/input/GKX_05.realigned.chr21.aroundvariants.bam.bai \
    -iregion=/test/kccg-performance-measure/input/reportable_range.chr21.bed \
    --folder /test/vcfscope-measure/output --yes --watch
# SUCCESS job-Bkkjv3j098GfBQ63GxY099vP
