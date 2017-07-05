import logging
from vcf_dedup.tools.vcf_transformer import StrelkaVcfDedupper, StarlingVcfDedupper, \
    DuplicationFinder, PlatypusVcfDedupper, GenericVcfDedupper
from vcf_dedup.tools.variant_comparer import VariantComparerNoAlternate, VariantComparerWithAlternate
from vcf_dedup.tools.vcf_sorter import VcfSorter
from vcf_dedup.constants import *
import os


class VcfDedupInputError(Exception):
    pass


class VcfDedupRunner(object):

    # TODO: add the rare diseases caller here, platypus?
    # TODO: add others
    SUPPORTED_VARIANT_CALLERS = [
        STRELKA_VARIANT_CALLER,
        STARLING_VARIANT_CALLER,
        DUPLICATION_FINDER,
        PLATYPUS_VARIANT_CALLER,
        GENERIC_VARIANT_CALLER
    ]
    SUPPORTED_EQUALITY_MODES = [
        INCLUDE_ALTERNATE,
        EXCLUDE_ALTERNATE
    ]
    SUPPORTED_SELECTION_METHODS = [
        SELECTION_METHOD_AF,
        SELECTION_METHOD_ALLELE_CALLS,
        SELECTION_METHOD_ARBITRARY,
        SELECTION_METHOD_QUALITY
    ]

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
        try:
            self.sample_idx = int(config["sample_idx"])
        except ValueError:
            raise VcfDedupInputError("'sample_idx' must be numeric")
        self.sample_name = config["sample_name"]
        try:
            self.sort_vcf = bool(config["sort_vcf"])
        except ValueError:
            raise VcfDedupInputError("'sort_vcf' must be boolean")
        try:
            self.sort_threads = int(config["sort_threads"])
            if self.sort_threads < 1:
                raise VcfDedupInputError("'sort_threads' must be positive and greater than zero [%s]" %
                                         str(self.sort_threads))
        except ValueError:
            raise VcfDedupInputError("'sort_threads' must be numeric")
        self.temp_folder = config["temp_folder"] if config["temp_folder"] != "" else \
            os.path.dirname(os.path.realpath(self.output_vcf if self.output_vcf is not None and self.output_vcf != ""
                                             else self.input_vcf))

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
        if "sample_name" not in self.config:
            raise VcfDedupInputError("Missing parameter 'sample_name'")
        if "sort_vcf" not in self.config:
            raise VcfDedupInputError("Missing parameter 'sort_vcf'")
        if "sort_threads" not in self.config:
            raise VcfDedupInputError("Missing parameter 'sort_threads'")
        if "temp_folder" not in self.config:
            raise VcfDedupInputError("Missing parameter 'temp_folder'")
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
        if self.config["variant_caller"] == "strelka" and self.config["selection_method"] == "allele_calls":
            raise VcfDedupInputError("Non supported selection method for Strelka which has no genotypes")

    def __sort(self):
        logging.info("Sorts the VCF before processing...")
        input_basename = os.path.splitext(os.path.basename(self.input_vcf))[0]
        self.sorted_vcf = os.path.join(self.temp_folder, "." + input_basename + ".sorted.vcf")
        vcf_sorter = VcfSorter(self.input_vcf, self.sorted_vcf, temp_folder=self.temp_folder, threads=self.sort_threads)
        vcf_sorter.sort()

    def run(self):
        logging.info("Starting the VCF duplication removal...")
        if self.sort_vcf:
            # sorts the input VCF
            self.__sort()
        else:
            # assumes input is sorted
            self.sorted_vcf = self.input_vcf
        # selects the appropriate comparer
        if self.equality_mode == INCLUDE_ALTERNATE:
            comparer = VariantComparerWithAlternate()
        elif self.equality_mode == EXCLUDE_ALTERNATE:
            comparer = VariantComparerNoAlternate()
        # selects the appropriate transformer
        if self.variant_caller == STRELKA_VARIANT_CALLER:
            transformer = StrelkaVcfDedupper(
                self.sorted_vcf,
                self.output_vcf,
                comparer,
                self.selection_method,
                self.sample_idx,
                self.sample_name
            )
        elif self.variant_caller == STARLING_VARIANT_CALLER:
            transformer = StarlingVcfDedupper(
                self.sorted_vcf,
                self.output_vcf,
                comparer,
                self.selection_method,
                self.sample_idx,
                self.sample_name
            )
        elif self.variant_caller == PLATYPUS_VARIANT_CALLER:
            transformer = PlatypusVcfDedupper(
                self.sorted_vcf,
                self.output_vcf,
                comparer,
                self.selection_method,
                self.sample_idx,
                self.sample_name
            )
        elif self.variant_caller == GENERIC_VARIANT_CALLER:
            transformer = GenericVcfDedupper(
                self.sorted_vcf,
                self.output_vcf,
                comparer,
                self.selection_method,
                self.sample_idx,
                self.sample_name
            )
        elif self.variant_caller == DUPLICATION_FINDER:
            transformer = DuplicationFinder(
                self.sorted_vcf,
                self.output_vcf,
                comparer)
        # run the transformation
        transformer.process_vcf()
        # forces closing the resources
        transformer.__del__()
        # delete temporary files
        if self.sort_vcf:
            os.remove(self.sorted_vcf)
        logging.info("Duplication removal finished!")

