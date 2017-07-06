import unittest
from vcf_dedup.tools.vcf_transformer import StrelkaVcfDedupper, StarlingVcfDedupper, \
    DuplicationFinder, VcfFormatError, GenericVcfDedupper
from vcf_dedup.tools.variant_comparer import VariantComparerNoAlternate, VariantComparerWithAlternate
from vcf_dedup.tools.vcf_sorter import VcfSorter
from vcf_dedup.constants import *
import logging
import os
import sys

class VcfDedupTests2(unittest.TestCase):

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

    #####
    ## Strelka VCF dedupper
    #####
    def test1_0(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(),
                                             SELECTION_METHOD_AF, 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_0_0(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(),
                                             SELECTION_METHOD_AF, 1, "tumor")
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_0_1(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(),
                                             SELECTION_METHOD_AF, 1, "none")
        try:
            vcf_transformer.process_vcf()
            self.assertTrue(False)
        except VcfFormatError:
            self.assertTrue(True)
        del vcf_transformer

    def test1_0_2(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StrelkaVcfDedupper(input_vcf, None, VariantComparerNoAlternate(),
                                             SELECTION_METHOD_AF, 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_1(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(),
                                             SELECTION_METHOD_AF, 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_2(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(),
                                             SELECTION_METHOD_QUALITY, 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_3(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(),
                                             SELECTION_METHOD_QUALITY, 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_4(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(),
                                             SELECTION_METHOD_ARBITRARY, 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_5(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(),
                                             SELECTION_METHOD_ARBITRARY, 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    #####
    ## Starling VCF dedupper
    #####
    def test2_0(self):
        input_vcf = self.starling_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(),
                                              SELECTION_METHOD_AF,
                                              sample_idx=0)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_1(self):
        input_vcf = self.starling_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(),
                                              SELECTION_METHOD_AF,
                                              sample_idx=0)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_2(self):
        input_vcf = self.starling_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(),
                                              SELECTION_METHOD_QUALITY,
                                              sample_idx=0)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_3(self):
        input_vcf = self.starling_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(),
                                              SELECTION_METHOD_QUALITY,
                                              sample_idx=0)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_4(self):
        input_vcf = self.starling_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(),
                                              SELECTION_METHOD_ARBITRARY,
                                              sample_idx=0)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_5(self):
        input_vcf = self.starling_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(),
                                              SELECTION_METHOD_ARBITRARY,
                                              sample_idx=0)
        vcf_transformer.process_vcf()
        del vcf_transformer

    #####
    ## Generic VCF dedupper
    #####
    def test5_0(self):
        input_vcf = self.generic1
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = GenericVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), SELECTION_METHOD_AF)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test5_1(self):
        input_vcf = self.generic1
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = GenericVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), SELECTION_METHOD_AF)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test5_2(self):
        input_vcf = self.generic1
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = GenericVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(),
                                             SELECTION_METHOD_QUALITY)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test5_3(self):
        input_vcf = self.generic1
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = GenericVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(),
                                             SELECTION_METHOD_QUALITY)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test5_4(self):
        input_vcf = self.generic1
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = GenericVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(),
                                             SELECTION_METHOD_ARBITRARY)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test5_5(self):
        input_vcf = self.generic1
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = GenericVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(),
                                             SELECTION_METHOD_ARBITRARY)
        vcf_transformer.process_vcf()
        del vcf_transformer

    #####
    ## Duplicate finder
    #####
    def test3_0(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_transformer = DuplicationFinder(input_vcf, output_vcf, VariantComparerWithAlternate())
        vcf_transformer.process_vcf()
        del vcf_transformer


    #####
    ##
    #####
    def test4_0(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_sorter = VcfSorter(input_vcf, output_vcf)
        vcf_sorter.sort()
        del vcf_sorter

    def test4_1(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_sorter = VcfSorter(input_vcf, output_vcf, threads=8)
        vcf_sorter.sort()
        del vcf_sorter

    def test4_2(self):
        input_vcf = self.strelka_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_sorter = VcfSorter(input_vcf, output_vcf, temp_folder="../resources/tmp")
        vcf_sorter.sort()
        del vcf_sorter

    def test4_3(self):
        input_vcf = self.starling_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_sorter = VcfSorter(input_vcf, output_vcf)
        vcf_sorter.sort()
        del vcf_sorter

    def test4_4(self):
        input_vcf = self.starling_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_sorter = VcfSorter(input_vcf, output_vcf, threads=8)
        vcf_sorter.sort()
        del vcf_sorter

    def test4_5(self):
        input_vcf = self.starling_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_sorter = VcfSorter(input_vcf, output_vcf, temp_folder="../resources/tmp")
        vcf_sorter.sort()
        del vcf_sorter

    def test4_6(self):
        input_vcf = self.unsorted_vcf
        class_name = type(self).__name__
        test_name = sys._getframe().f_code.co_name
        output_vcf = os.path.join(self.OUTPUT_FOLDER, "%s.%s.vcf" % (class_name, test_name))
        vcf_sorter = VcfSorter(input_vcf, output_vcf)
        vcf_sorter.sort()
        del vcf_sorter