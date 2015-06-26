# DNANexus Tests
## HiSeqX V2 (R_150203_DAVMIL1_FGS_M001), chromosome 21 only, standard report
dx run /kccg-validation-reporter -ivcfgz=file-Bf24XjQ0jxF5fJ7Gb83Fvp2f -iregion=file-Bf4JFV00K74xkqBP2Qzbq1p4 --yes --ssh --debug-on AppInternalError

## HiSeqX V2 (R_150203_DAVMIL1_FGS_M001), whole genome, extended report
dx run /kccg-validation-reporter -ivcfgz=file-Bf24XjQ0jxF5fJ7Gb83Fvp2f -iregion=file-Bf4JFV00K74xkqBP2Qzbq1p4 -iextended=true --yes --ssh --debug-on AppInternalError
