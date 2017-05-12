import unittest
from vcf_dedup.tools.vcf_transformer import StrelkaVcfDedupper, StarlingVcfDedupper, DuplicationFinder, VcfFormatError
from vcf_dedup.tools.variant_comparer import VariantComparerNoAlternate, VariantComparerWithAlternate
from vcf_dedup.tools.vcf_sorter import VcfSorter
import logging

class VcfDedupTests(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.starling_vcf = "../../resources/starling_duplicated_variants.vcf"
        self.strelka_vcf = "../../resources/strelka_duplicated_variants.vcf"
        self.unsorted_vcf = "../../resources/duplicated_unsorted.vcf"

    #####
    ## Strelka VCF dedupper
    #####
    def test1_0(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test1_0.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), "af", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_0_0(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test1_0_0.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), "af", 1, "tumor")
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_0_1(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test1_0_1.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), "af", 1, "none")
        try:
            vcf_transformer.process_vcf()
            self.assertTrue(False)
        except VcfFormatError:
            self.assertTrue(True)
        del vcf_transformer

    def test1_0_2(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test1_0_2.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, None, VariantComparerNoAlternate(), "af", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_1(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test1_1.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), "af", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_2(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test1_2.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), "quality", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_3(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test1_3.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), "quality", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_4(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test1_4.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), "arbitrary", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_5(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test1_5.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), "arbitrary", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    #####
    ## Starling VCF dedupper
    #####
    def test2_0(self):
        input_vcf = self.starling_vcf
        output_vcf = "../../resources/test2_0.vcf"
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), "af")
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_1(self):
        input_vcf = self.starling_vcf
        output_vcf = "../../resources/test2_1.vcf"
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), "af")
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_2(self):
        input_vcf = self.starling_vcf
        output_vcf = "../../resources/test2_2.vcf"
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), "quality")
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_3(self):
        input_vcf = self.starling_vcf
        output_vcf = "../../resources/test2_3.vcf"
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), "quality")
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_4(self):
        input_vcf = self.starling_vcf
        output_vcf = "../../resources/test2_4.vcf"
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), "arbitrary")
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_5(self):
        input_vcf = self.starling_vcf
        output_vcf = "../../resources/test2_5.vcf"
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(),
                                                   "arbitrary")
        vcf_transformer.process_vcf()
        del vcf_transformer

    #####
    ## Duplicate finder
    #####
    def test3_0(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test3_0.vcf"
        vcf_transformer = DuplicationFinder(input_vcf, output_vcf, VariantComparerWithAlternate())
        vcf_transformer.process_vcf()
        del vcf_transformer


    #####
    ##
    #####
    def test4_0(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test4_0.vcf"
        vcf_sorter = VcfSorter(input_vcf, output_vcf)
        vcf_sorter.sort()
        del vcf_sorter

    def test4_1(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test4_1.vcf"
        vcf_sorter = VcfSorter(input_vcf, output_vcf, threads=8)
        vcf_sorter.sort()
        del vcf_sorter

    def test4_2(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../../resources/test4_2.vcf"
        vcf_sorter = VcfSorter(input_vcf, output_vcf, temp_folder="../resources/tmp")
        vcf_sorter.sort()
        del vcf_sorter

    def test4_3(self):
        input_vcf = self.starling_vcf
        output_vcf = "../../resources/test4_3.vcf"
        vcf_sorter = VcfSorter(input_vcf, output_vcf)
        vcf_sorter.sort()
        del vcf_sorter

    def test4_4(self):
        input_vcf = self.starling_vcf
        output_vcf = "../../resources/test4_4.vcf"
        vcf_sorter = VcfSorter(input_vcf, output_vcf, threads=8)
        vcf_sorter.sort()
        del vcf_sorter

    def test4_5(self):
        input_vcf = self.starling_vcf
        output_vcf = "../../resources/test4_5.vcf"
        vcf_sorter = VcfSorter(input_vcf, output_vcf, temp_folder="../resources/tmp")
        vcf_sorter.sort()
        del vcf_sorter

    def test4_6(self):
        input_vcf = self.unsorted_vcf
        output_vcf = "../../resources/test4_6.vcf"
        vcf_sorter = VcfSorter(input_vcf, output_vcf)
        vcf_sorter.sort()
        del vcf_sorter