import os
import unittest
import re

import dxpy


class SmokeTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project_id = str(dxpy.PROJECT_CONTEXT_ID)
        cls.project_name = os.environ.get("PROJ_NAME", None)
        cls.project = dxpy.DXProject(dxpy.find_one_project(zero_ok=False, more_ok=False, name=cls.project_name)["id"].encode('ascii', 'ignore'))
        for output_file in cls.project.list_folder(describe=True)["objects"]:
            description = output_file["describe"]
            if description["name"].endswith("pdf"):
                dxpy.download_dxfile(description["id"], "report.pdf")

    # Output pdf must exist and be larger than 0 bytes.
    def test_check_output_files_made(self):
        for output_file in self.project.list_folder(describe=True)["objects"]:
            description = output_file["describe"]
            if description["name"].endswith("pdf"):
                self.assertGreater(description["size"], 0)
                break
        else:
            self.fail('No pdf file produced.')
