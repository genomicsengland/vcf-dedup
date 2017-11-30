import os
import sys
import unittest
from vcf_dedup.runner import VcfDedupRunner
from vcf_dedup.constants import *
import logging


class VcfDedupTests(unittest.TestCase):

    INPUT_FOLDER = "resources/input"
    OUTPUT_FOLDER = "resources/output"

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.starling_vcf = os.path.join(self.INPUT_FOLDER, "starling_duplicated_variants.vcf")
        self.strelka_vcf = os.path.join(self.INPUT_FOLDER, "strelka_duplicated_variants.vcf")
        self.unsorted_vcf = os.path.join(self.INPUT_FOLDER, "duplicated_unsorted.vcf")
        self.platypus_vcf = os.path.join(self.INPUT_FOLDER, "platypus_duplicated_variants.vcf")
        self.generic1 = os.path.join(self.INPUT_FOLDER, "duplicates.chr11.no_head.vcf")
        self.generic2 = os.path.join(self.INPUT_FOLDER, "uk10k.dup.no_head.vcf")
        self.no_duplicated = os.path.join(self.INPUT_FOLDER, "vcf_without_duplicated_variants.vcf")

    def test1_0(self):
        input_vcf = self.starling_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": STRELKA_VARIANT_CALLER,
            "selection_method": SELECTION_METHOD_AF,
            "equality_mode": INCLUDE_ALTERNATE,
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
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": STRELKA_VARIANT_CALLER,
            "selection_method": SELECTION_METHOD_AF,
            "equality_mode": INCLUDE_ALTERNATE,
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
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": STRELKA_VARIANT_CALLER,
            "selection_method": SELECTION_METHOD_AF,
            "equality_mode": INCLUDE_ALTERNATE,
            "sample_idx": "0",
            "sample_name": "sample",
            "sort_vcf": False,
            "sort_threads": -2,
            "temp_folder": ""
        }
        try:
            print 1
            vcf_dedup_runner = VcfDedupRunner(config)
            print 2
            vcf_dedup_runner.run()
            print 3
            self.assertTrue(False)
            print 4
        except:
            self.assertTrue(True)

    def test1_3(self):
        input_vcf = self.platypus_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": PLATYPUS_VARIANT_CALLER,
            "selection_method": SELECTION_METHOD_ALLELE_CALLS,
            "equality_mode": INCLUDE_ALTERNATE,
            "sort_vcf": True,
            "sort_threads": 1,
            "temp_folder": ""
        }
        vcf_dedup_runner = VcfDedupRunner(config)
        vcf_dedup_runner.run()
        del vcf_dedup_runner

    def test1_4(self):
        input_vcf = self.starling_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": STARLING_VARIANT_CALLER,
            "selection_method": SELECTION_METHOD_ALLELE_CALLS,
            "equality_mode": INCLUDE_ALTERNATE,
            "sample_idx": "0",
            "sample_name": "sample",
            "sort_vcf": True,
            "sort_threads": 1,
            "temp_folder": ""
        }
        vcf_dedup_runner = VcfDedupRunner(config)
        vcf_dedup_runner.run()
        del vcf_dedup_runner

    def test1_5(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": STRELKA_VARIANT_CALLER,
            "selection_method": SELECTION_METHOD_ALLELE_CALLS,
            "equality_mode": INCLUDE_ALTERNATE,
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

    def test1_6(self):
        input_vcf = self.generic1
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": GENERIC_VARIANT_CALLER,
            "selection_method": SELECTION_METHOD_AF,
            "equality_mode": INCLUDE_ALTERNATE,
            "sort_vcf": True,
            "sort_threads": 1,
            "temp_folder": ""
        }
        vcf_dedup_runner = VcfDedupRunner(config)
        vcf_dedup_runner.run()
        del vcf_dedup_runner

    def test1_7(self):
        input_vcf = self.generic2
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": GENERIC_VARIANT_CALLER,
            "selection_method": SELECTION_METHOD_AF,
            "equality_mode": INCLUDE_ALTERNATE,
            "sort_vcf": True,
            "sort_threads": 1,
            "temp_folder": ""
        }
        vcf_dedup_runner = VcfDedupRunner(config)
        vcf_dedup_runner.run()
        del vcf_dedup_runner

    def test1_8(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": STARLING_VARIANT_CALLER,
            "selection_method": SELECTION_METHOD_AF,
            "equality_mode": INCLUDE_ALTERNATE,
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

    def test1_9(self):
        input_vcf = self.no_duplicated
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        config = {
            "input_vcf": input_vcf,
            "output_vcf": output_vcf,
            "variant_caller": GENERIC_VARIANT_CALLER,
            "selection_method": SELECTION_METHOD_AF,
            "equality_mode": INCLUDE_ALTERNATE,
            "sort_vcf": True,
            "sort_threads": 1,
            "temp_folder": ""
        }
        try:
            vcf_dedup_runner = VcfDedupRunner(config)
            vcf_dedup_runner.run()
            self.assertTrue(True)
        except ValueError, ex:
            print ex
            self.assertTrue(False)