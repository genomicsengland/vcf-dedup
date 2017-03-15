import unittest
from vcf_dedup.tools.vcf_transformer import VcfDedupper
import logging

class VcfDedupTests(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.input_vcf = "../resources/LP2000906-DNA_F01.somatic.atomic.left.split.vcf.gz"
        self.output_vcf = "../resources/LP2000906-DNA_F01.somatic.transformed.atomic.left.split.vcf"



    def test1(self):
        self.vcf_transformer = VcfDedupper(self.input_vcf, self.output_vcf)
        self.vcf_transformer.process_vcf()