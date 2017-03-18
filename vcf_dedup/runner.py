import logging
from vcf_dedup.tools.vcf_transformer import StrelkaVcfDedupper, StarlingVcfDedupper, DuplicationFinder
from vcf_dedup.tools.variant_comparer import VariantComparerNoAlternate, VariantComparerWithAlternate


class VcfDedupInputError(Exception):
    pass


class VcfDedupRunner(object):

    # TODO: add the rare diseases caller here, platypus?
    # TODO: add others
    SUPPORTED_VARIANT_CALLERS = ["strelka", "starling", "duplication_finder"]
    SUPPORTED_EQUALITY_MODES = ["1", "2"]
    SUPPORTED_SELECTION_METHODS = ["af", "quality", "arbitrary"]

    def __init__(self, config):
        self.config = config
        # checks that the configuration received is correct
        self.sanity_checks()
        # loads configuration data
        self.input_vcf = config["input_vcf"]
        self.output_vcf = config["output_vcf"]
        self.variant_caller = config["variant_caller"]
        self.selection_method = config["selection_method"]
        self.equality_mode = config["equality_mode"]
        self.sample_idx = config["sample_idx"]

    def sanity_checks(self):
        """
        Checks the input configuration and throws an error if necessary
        :return:
        """
        if "input_vcf" not in self.config:
            raise VcfDedupInputError("Missing parameter 'input_vcf'")
        if "output_vcf" not in self.config:
            raise VcfDedupInputError("Missing parameter 'output_vcf'")
        if "variant_caller" not in self.config:
            raise VcfDedupInputError("Missing parameter 'variant_caller'")
        if "selection_method" not in self.config:
            raise VcfDedupInputError("Missing parameter 'selection_method'")
        if "equality_mode" not in self.config:
            raise VcfDedupInputError("Missing parameter 'equality_mode'")
        if "sample_idx" not in self.config:
            raise VcfDedupInputError("Missing parameter 'sample_idx'")
        if self.config["variant_caller"] not in self.SUPPORTED_VARIANT_CALLERS:
            raise VcfDedupInputError("Non supported variant caller [%s]. The list of supported variant callers is %s" %
                                     (self.config["variant_caller"], ", ".join(self.SUPPORTED_VARIANT_CALLERS)))
        if self.config["selection_method"] not in self.SUPPORTED_SELECTION_METHODS:
            raise VcfDedupInputError("Non supported selection method [%s]. The list of supported selection methods is %s" %
                                     (self.config["selection_method"], ", ".join(self.SUPPORTED_SELECTION_METHODS)))
        if self.config["equality_mode"] not in self.SUPPORTED_EQUALITY_MODES:
            raise VcfDedupInputError(
                "Non supported variant equality method [%s]. The list of supported equality methods is %s" %
                (self.config["variant_caller"], ", ".join(self.SUPPORTED_EQUALITY_MODES)))

    def run(self):
        logging.info("Starting the VCF duplication removal...")
        # selects the appropriate comparer
        if self.equality_mode == "1":
            comparer = VariantComparerWithAlternate()
        elif self.equality_mode == "2":
            comparer = VariantComparerNoAlternate()
        # selects the appropriate transformer
        if self.variant_caller == "strelka":
            transformer = StrelkaVcfDedupper(
                self.input_vcf,
                self.output_vcf,
                comparer,
                self.selection_method,
                self.tumor_sample_idx)
        elif self.variant_caller == "starling":
            transformer = StarlingVcfDedupper(
                self.input_vcf,
                self.output_vcf,
                comparer,
                self.selection_method)
        elif self.variant_caller == "duplication_finder":
            transformer = DuplicationFinder(
                self.input_vcf,
                self.output_vcf,
                comparer)
        # run the transformation
        transformer.process_vcf()

