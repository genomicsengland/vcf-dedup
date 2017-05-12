import unittest
from vcf_dedup.runner import VcfDedupRunner
import logging


class VcfDedupTests(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.starling_vcf = "../resources/starling_duplicated_variants.vcf"
        self.strelka_vcf = "../resources/strelka_duplicated_variants.vcf"
        self.unsorted_vcf = "../resources/duplicated_unsorted.vcf"
        self.platypus_vcf = "../resources/platypus_duplicated_variants.vcf"

    def test1_0(self):
        input_vcf = self.starling_vcf
        output_vcf = "../resources/testVcfDedup1_0.vcf"
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": "strelka",
            "selection_method": "af",
            "equality_mode": "1",
            "sample_idx": "0",
            "sample_name": "sample",
            "sort_vcf": True,
            "sort_threads": 1,
            "temp_folder": ""
        }
        vcf_dedup_runner = VcfDedupRunner(config)
        vcf_dedup_runner.run()
        del vcf_dedup_runner

    def test1_1(self):
        input_vcf = self.starling_vcf
        output_vcf = "../resources/testVcfDedup1_1.vcf"
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": "strelka",
            "selection_method": "af",
            "equality_mode": "1",
            "sample_idx": "0",
            "sample_name": "sample",
            "sort_vcf": False,
            "sort_threads": 1,
            "temp_folder": ""
        }
        vcf_dedup_runner = VcfDedupRunner(config)
        vcf_dedup_runner.run()
        del vcf_dedup_runner

    def test1_2(self):
        input_vcf = self.starling_vcf
        output_vcf = "../resources/testVcfDedup1_2.vcf"
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": "strelka",
            "selection_method": "af",
            "equality_mode": "1",
            "sample_idx": "0",
            "sample_name": "sample",
            "sort_vcf": False,
            "sort_threads": -2,
            "temp_folder": ""
        }
        try:
            vcf_dedup_runner = VcfDedupRunner(config)
            vcf_dedup_runner.run()
            self.assertTrue(False)
        except:
            self.assertTrue(True)

    def test1_3(self):
        input_vcf = self.platypus_vcf
        output_vcf = "../resources/testVcfDedup1_3.vcf"
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": "platypus",
            "selection_method": "allele_calls",
            "equality_mode": "1",
            "sample_idx": "0",
            "sample_name": "sample",
            "sort_vcf": True,
            "sort_threads": 1,
            "temp_folder": ""
        }
        vcf_dedup_runner = VcfDedupRunner(config)
        vcf_dedup_runner.run()

    def test1_4(self):
        input_vcf = self.starling_vcf
        output_vcf = "../resources/testVcfDedup1_4.vcf"
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": "starling",
            "selection_method": "allele_calls",
            "equality_mode": "1",
            "sample_idx": "0",
            "sample_name": "sample",
            "sort_vcf": True,
            "sort_threads": 1,
            "temp_folder": ""
        }
        vcf_dedup_runner = VcfDedupRunner(config)
        vcf_dedup_runner.run()

    def test1_5(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../resources/testVcfDedup1_5.vcf"
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": "strelka",
            "selection_method": "allele_calls",
            "equality_mode": "1",
            "sample_idx": "0",
            "sample_name": "tumor",
            "sort_vcf": True,
            "sort_threads": 1,
            "temp_folder": ""
        }
        try:
            vcf_dedup_runner = VcfDedupRunner(config)
            vcf_dedup_runner.run()
            self.assertTrue(False)
        except:
            self.assertTrue(True)
