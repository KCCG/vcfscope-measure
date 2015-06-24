# DNANexus Tests
## HiSeqX V2 (R_150203_DAVMIL1_FGS_M001), chromosome 21 only, standard report
dx run /kccg-validation-reporter-dx -ivcfgz=project-Bb9KVk8029vp1qzXz4yx4xB3:file-BbBg1g00z1f8vYJq606Yjv7K -iregion=project-BZ4JvjQ0K74XK3bP71gykXKQ:file-Bf4JFV00K74xkqBP2Qzbq1p4 --yes --ssh --debug-on AppInternalError

## HiSeqX V2 (R_150203_DAVMIL1_FGS_M001), whole genome, extended report
dx run /kccg-validation-reporter-dx -ivcfgz=project-Bb9KVk8029vp1qzXz4yx4xB3:file-BbBg1g00z1f8vYJq606Yjv7K -iextended=true --yes --ssh --debug-on AppInternalError
