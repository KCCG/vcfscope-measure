import os
import unittest
import re
import dxpy
import sys
import subprocess

app_name = "kccg-performance-measure"
project_name = os.getenv('PROJ_NAME')
kccg_project_id = "project-BVJz7k0098GX43GV9ZPFXVY2"
RUN_JOB_ON_DX = os.getenv('RUN_JOB_ON_DX', "True") != "False"

#
# To test this locally, without running the job on DX:
#
# test -d venv || python virtualenv-1.11.6/virtualenv.py -p python2.7 venv
# source venv/bin/activate
# pip install -r requirements.txt
#
# export RUN_JOB_ON_DX=False
# nosetests
#
class TestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if RUN_JOB_ON_DX:
            if not project_name:
                print "'PROJ_NAME' environment variable must be defined!"
                sys.exit(1)
            project_id = dxpy.find_one_project(more_ok=False, name=project_name)["id"]
            run_args = {}
            run_args["project"] = project_id
            run_args["name"] = "kccg-performance-measure on chr21"
            run_args["folder"] = "/purge/" + app_name
            input_hash = {}
            input_hash["vcfgz"] = dxpy.dxlink("file-BkkjFkj098Gb2jZ1Yx533JFv", kccg_project_id)
            input_hash["bam"] = dxpy.dxlink("file-Bkkjj5Q098Gkvkb3Xx5Pxj1J", kccg_project_id)
            input_hash["bai"] = dxpy.dxlink("file-Bkkjj5Q098GzYx2bG5YJ3z34", kccg_project_id)
            input_hash["region"] = dxpy.dxlink("file-Bkkj22Q098Gz5yK1Q955G5gX", kccg_project_id)

            app = dxpy.DXApp(name=app_name, alias="9.9.7")
            cls.job = app.run(input_hash, **run_args)

        else:
            job_id = "job-Bkkjv3j098GfBQ63GxY099vP"
            cls.job = dxpy.DXJob(job_id)

        cls.job.wait_on_done()

    def test_AAA_DownloadResultResults(self):
        job_hash = self.job.describe()
        output_hash = job_hash["output"]["rds"]
        for output in output_hash:
            f = dxpy.DXFile(output, project=job_hash["project"])
            print "TestCase: Downloading %s" % f.name
            dxpy.download_dxfile(f.id, f.name, project=job_hash["project"])
            self.assertTrue(os.path.isfile(f.name))

    @classmethod
    def tearDownClass(self):
        os.remove("GKX_30.realigned.recalibrated.hc.vqsr.chr21.perfmeas.rds")


    def test_check_perf_table(self):
        subprocess.call(['R', '--vanilla', '-e', 'write.csv(readRDS("GKX_30.realigned.recalibrated.hc.vqsr.chr21.perfmeas.rds")$class_subsets.performance_thresholded, file = "GKX_30.realigned.recalibrated.hc.vqsr.chr21.perfmeas.csv")'])
        md5sum = os.popen("md5sum GKX_30.realigned.recalibrated.hc.vqsr.chr21.perfmeas.csv | sed 's/ .*//'").read().strip()
        self.assertEqual(md5sum, "8114e3688e6b3d0e7e3f2365af5c540c")
