import logging
from vcf_dedup.tools.vcf_transformer import StrelkaVcfDedupper, StarlingVcfDedupper, DuplicationFinder
from vcf_dedup.tools.variant_comparer import VariantComparerNoAlternate, VariantComparerWithAlternate


class VcfDedupInputError(Exception):
    pass


class VcfDedupRunner(object):

    # TODO: add the rare diseases caller here, platypus?
    # TODO: add others
    SUPPORTED_VARIANT_CALLERS = ["strelka", "starling", "duplication_finder"]
    SUPPORTED_EQUALITY_MODE = ["1", "2"]

    def __init__(self, config):
        self.config = config
        self.input_vcf = config["input_vcf"]
        self.output_vcf = config["output_vcf"]
        self.variant_caller = config["variant_caller"]
        self.equality_mode = config["equality_mode"]

        # checks that configuration received is correct
        self.sanity_checks()

    def sanity_checks(self):
        """
        Checks the input configuration and throws an error if necessary
        :return:
        """
        if self.variant_caller not in self.SUPPORTED_VARIANT_CALLERS:
            raise VcfDedupInputError("Non supported variant caller [%s]. The list of supported variant callers is %s" %
                                     (self.variant_caller, ", ".join(self.SUPPORTED_VARIANT_CALLERS)))
        if self.equality_mode not in self.SUPPORTED_EQUALITY_MODE:
            raise VcfDedupInputError("Non supported variant equality method [%s]. The list of supported equality methods is %s" %
                                     (self.variant_caller, ", ".join(self.SUPPORTED_EQUALITY_MODE)))

    def run(self):
        logging.info("Starting the VCF duplication removal...")
        # selects the appropriate comparer
        if self.equality_mode == "1":
            comparer = VariantComparerWithAlternate()
        elif self.equality_mode == "2":
            comparer = VariantComparerNoAlternate()
        # selects the appropriate transformer
        if self.variant_caller == "strelka":
            transformer = StrelkaVcfDedupper(self.input_vcf, self.output_vcf, comparer)
        elif self.variant_caller == "starling":
            transformer = StarlingVcfDedupper(self.input_vcf, self.output_vcf, comparer)
        elif self.variant_caller == "duplication_finder":
            transformer = DuplicationFinder(self.input_vcf, self.output_vcf, comparer)
        # run the transformation
        transformer.process_vcf()

