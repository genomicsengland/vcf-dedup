import unittest
from vcf_dedup.tools.vcf_transformer import StrelkaVcfDedupper, StarlingVcfDedupper, DuplicationFinder
from vcf_dedup.tools.variant_comparer import VariantComparerNoAlternate, VariantComparerWithAlternate
import logging

class VcfDedupTests(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.starling_vcf = "../resources/starling_duplicated_variants.vcf"
        self.strelka_vcf = "../resources/strelka_duplicated_variants.vcf"

    #####
    ## Strelka VCF dedupper
    #####
    def test1_0(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../resources/test1_0.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), "af", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_1(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../resources/test1_1.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), "af", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_2(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../resources/test1_2.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), "quality", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_3(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../resources/test1_3.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), "quality", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_4(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../resources/test1_4.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), "arbitrary", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test1_5(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../resources/test1_5.vcf"
        vcf_transformer = StrelkaVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), "arbitrary", 1)
        vcf_transformer.process_vcf()
        del vcf_transformer

    #####
    ## Starling VCF dedupper
    #####
    def test2_0(self):
        input_vcf = self.starling_vcf
        output_vcf = "../resources/test2_0.vcf"
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), "af")
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_1(self):
        input_vcf = self.starling_vcf
        output_vcf = "../resources/test2_1.vcf"
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), "af")
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_2(self):
        input_vcf = self.starling_vcf
        output_vcf = "../resources/test2_2.vcf"
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(), "quality")
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_3(self):
        input_vcf = self.starling_vcf
        output_vcf = "../resources/test2_3.vcf"
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), "quality")
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_4(self):
        input_vcf = self.starling_vcf
        output_vcf = "../resources/test2_4.vcf"
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerWithAlternate(), "arbitrary")
        vcf_transformer.process_vcf()
        del vcf_transformer

    def test2_5(self):
        input_vcf = self.starling_vcf
        output_vcf = "../resources/test2_5.vcf"
        vcf_transformer = StarlingVcfDedupper(input_vcf, output_vcf, VariantComparerNoAlternate(),
                                                   "arbitrary")
        vcf_transformer.process_vcf()
        del vcf_transformer

    #####
    ## Duplicate finder
    #####
    def test3_0(self):
        input_vcf = self.strelka_vcf
        output_vcf = "../resources/test3_0.vcf"
        vcf_transformer = DuplicationFinder(input_vcf, output_vcf, VariantComparerWithAlternate())
        vcf_transformer.process_vcf()
        del vcf_transformer