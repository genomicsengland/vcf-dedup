import unittest
from vcf_dedup.tools.vcf_transformer import VcfDedupper
from vcf_dedup.tools.variant_comparer import VariantComparerNoAlternate, VariantComparerWithAlternate
import logging

class VcfDedupTests(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)


    def test1_0(self):
        self.input_vcf = "../resources/LP2000906-DNA_F01.somatic.atomic.left.split.vcf.gz"
        self.output_vcf = "../resources/test1_0.vcf"
        self.vcf_transformer = VcfDedupper(self.input_vcf, self.output_vcf, VariantComparerNoAlternate())
        self.vcf_transformer.process_vcf()

    def test1_1(self):
        self.input_vcf = "../resources/LP2000906-DNA_F01.somatic.atomic.left.split.vcf"
        self.output_vcf = "../resources/test1_1.vcf"
        self.vcf_transformer = VcfDedupper(self.input_vcf, self.output_vcf, VariantComparerNoAlternate())
        self.vcf_transformer.process_vcf()

    def test2(self):
        self.input_vcf = "../resources/LP2000906-DNA_F01.somatic.atomic.left.split.vcf.gz"
        self.output_vcf = "../resources/test2_0.vcf"
        self.vcf_transformer = VcfDedupper(self.input_vcf, self.output_vcf, VariantComparerWithAlternate())
        self.vcf_transformer.process_vcf()

    def test2_1(self):
        self.input_vcf = "../resources/LP2000906-DNA_F01.somatic.atomic.left.split.vcf"
        self.output_vcf = "../resources/test2_1.vcf"
        self.vcf_transformer = VcfDedupper(self.input_vcf, self.output_vcf, VariantComparerWithAlternate())
        self.vcf_transformer.process_vcf()