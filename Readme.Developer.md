# tests
## chr 22
dx run /kccg-validation-reporter-dx -ivcfgz=file-BbBg1g00z1f8vYJq606Yjv7K -itbi=file-BbBg1g00z1f27y9Y0Q6vfb02 -idebug=22 --yes --ssh --debug-on AppInternalError
dx ssh job-Bf15zx00p67J7g4F31ky8zQV
dx terminate job-Bf15zx00p67J7g4F31ky8zQV

<!-- qsub -pe smp 4 -N VR_1_v10 -cwd -j y -b y bash validation_report.sh -f -x -o v10.pdf /directflow/ClinicalGenomicsPipeline/projects/validation-reporter/test_data/HiSeqX_v1_TKCC/NA12878_v1.hc.vqsr.vep.vcf.gz
qsub -pe smp 4 -N VR_1_v25 -cwd -j y -b y bash validation_report.sh -f -x -o v25.pdf /directflow/ClinicalGenomicsPipeline/projects/validation-reporter/test_data/HiSeqX_v2_TKCC/R_150203_DAVMIL1_FGS_M001.hc.vqsr.vep.vcf.gz
qsub -pe smp 4 -N VR_1_vEx -cwd -j y -b y bash validation_report.sh -f -o vEx.pdf -x /directflow/ClinicalGenomicsPipeline/projects/validation-reporter/test_data/HiSeq2500_NexteraRapid_Illumina/Basespace_NA12878_HiSeq_2500_Nextera_Rapid_Capture_Exome_CEPH_TRIO.hc.vqsr.vep.vcf.gz

genome = "BSgenome.HSapiens.1000g.37d5"
extendedflag = 1
DEBUG=1
DEBUG.chrom="21"
path.input = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/test_data/HiSeq2500_NexteraRapid_Illumina/Basespace_NA12878_HiSeq_2500_Nextera_Rapid_Capture_Exome_CEPH_TRIO.hc.vqsr.vep.vcf.gz"
path.tp = "../overlap/tp.vcf.gz"
path.fp = "../overlap/fp.vcf.gz"
path.fn = "../overlap/fn.vcf.gz"
path.gold = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/gold_standard/calls-2.19.vcf.gz"
path.gold_regions = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/gold_standard/valid_regions-2.19.bed.gz"
path.function_regions_prefix = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/functional_regions/"
path.mask_regions_prefix = "/directflow/ClinicalGenomicsPipeline/projects/validation-reporter/resources/mask_regions/"
 

Result 1:
ID              file-BbBg1g00z1f8vYJq606Yjv7K
Class           file
Project         project-Bb9KVk8029vp1qzXz4yx4xB3
Folder          /variants
Name            R_150203_DAVMIL1_FGS_M001.hc.vqsr.vcf.gz
State           closed
Visibility      visible
Types           -
Properties      -
Tags            -
Outgoing links  -
Created         Fri May 15 14:06:52 2015
Created by      joecop
 via the job    job-Bb9yfpQ029vVkjkQ0Zz9B5z5
Last modified   Fri May 15 14:07:04 2015
archivalState   null
Media type      application/x-gzip
Size            253.19 MB


Result 1:
ID              file-BbBg1g00z1f27y9Y0Q6vfb02
Class           file
Project         project-Bb9KVk8029vp1qzXz4yx4xB3
Folder          /variants
Name            R_150203_DAVMIL1_FGS_M001.hc.vqsr.vcf.gz.tbi
State           closed
Visibility      visible
Types           -
Properties      -
Tags            -
Outgoing links  -
Created         Fri May 15 14:06:52 2015
Created by      joecop
 via the job    job-Bb9yfpQ029vVkjkQ0Zz9B5z5
Last modified   Fri May 15 14:07:04 2015
archivalState   null
Media type      application/x-gzip
Size            1.58 MB
-->
