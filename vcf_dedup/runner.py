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

    SORT_MEM_PERCENTAGE_DEFAULT = 80
    LOGGING_FORMAT = '%(asctime)s %(levelname)s %(message)s'

    def __init__(self, config):
        self.config = config

        # Sets logging level
        log_level = 20 if self.config.get('verbose') is True else 40
        log_file = self.config.get('log_file')
        if log_file != "" and log_file is not None:
            logging.basicConfig(level=log_level, filename=log_file, format=self.LOGGING_FORMAT)
        else:
            logging.basicConfig(level=log_level, format=self.LOGGING_FORMAT)
        logging.debug("Configuration:")
        logging.debug(str(self.config))
        # checks that the configuration received is correct
        self.sanity_checks()
        # loads configuration data
        self.input_vcf = config["input_vcf"]
        self.output_vcf = config["output_vcf"]
        self.variant_caller = config["variant_caller"]
        self.selection_method = config["selection_method"]
        self.equality_mode = config["equality_mode"]
        self.sample_idx = int(config["sample_idx"]) if "sample_idx" in config else None
        self.sample_name = config["sample_name"] if "sample_name" in config else None
        self.sort_vcf = bool(config["sort_vcf"])
        self.sort_threads = int(config["sort_threads"])
        if config["temp_folder"] is not None:
            self.temp_folder = config["temp_folder"]
        else:
            self.temp_folder = os.path.dirname(os.path.realpath(
                self.output_vcf if self.output_vcf is not None and self.output_vcf != "" else self.input_vcf
            ))
        if config.get("sort_mem_percentage") is not None:
            self.sort_mem_percentage = int(config["sort_mem_percentage"])
        else:
            self.sort_mem_percentage = self.SORT_MEM_PERCENTAGE_DEFAULT

    def sanity_checks(self):
        """
        Checks the input configuration and throws an error if necessary
        :return:
        """
        if "input_vcf" not in self.config:
            message = "Missing parameter 'input_vcf'"
            logging.error(message)
            raise VcfDedupInputError(message)
        if "output_vcf" not in self.config:
            message = "Missing parameter 'output_vcf'"
            logging.error(message)
            raise VcfDedupInputError(message)
        if "variant_caller" not in self.config:
            message = "Missing parameter 'variant_caller'"
            logging.error(message)
            raise VcfDedupInputError(message)
        if "selection_method" not in self.config:
            message = "Missing parameter 'selection_method'"
            logging.error(message)
            raise VcfDedupInputError(message)
        if "equality_mode" not in self.config:
            message = "Missing parameter 'equality_mode'"
            logging.error(message)
            raise VcfDedupInputError(message)
        if "sort_vcf" not in self.config:
            message = "Missing parameter 'sort_vcf'"
            logging.error(message)
            raise VcfDedupInputError(message)
        try:
            sort_vcf = bool(self.config["sort_vcf"])
        except ValueError:
            message = "'sort_vcf' must be boolean and found [%s]" % str(sort_vcf)
            logging.error(message)
            raise VcfDedupInputError(message)
        if "sort_threads" not in self.config:
            message = "Missing parameter 'sort_threads'"
            logging.error(message)
            raise VcfDedupInputError(message)
        try:
            sort_threads = int(self.config["sort_threads"])
            if sort_threads < 1:
                message = "'sort_threads' must be positive and greater than zero [%s]" % str(sort_threads)
                logging.error(message)
                raise VcfDedupInputError(message)
        except ValueError:
            message = "'sort_threads' must be numeric"
            logging.error(message)
            raise VcfDedupInputError(message)
        if "temp_folder" not in self.config:
            message = "Missing parameter 'temp_folder'"
            logging.error(message)
            raise VcfDedupInputError(message)
        if self.config["variant_caller"] not in self.SUPPORTED_VARIANT_CALLERS:
            message = "Non supported variant caller [%s]. The list of supported variant callers is %s" % \
                      (self.config["variant_caller"], ", ".join(self.SUPPORTED_VARIANT_CALLERS))
            logging.error(message)
            raise VcfDedupInputError(message)
        if self.config["selection_method"] not in self.SUPPORTED_SELECTION_METHODS:
            message = "Non supported selection method [%s]. The list of supported selection methods is %s" %\
                      (self.config["selection_method"], ", ".join(self.SUPPORTED_SELECTION_METHODS))
            logging.error(message)
            raise VcfDedupInputError(message)
        if self.config["equality_mode"] not in self.SUPPORTED_EQUALITY_MODES:
            message = "Non supported variant equality method [%s]. The list of supported equality methods is %s" %\
                      (self.config["variant_caller"], ", ".join(self.SUPPORTED_EQUALITY_MODES))
            logging.error(message)
            raise VcfDedupInputError(message)
        if self.config["variant_caller"] == STRELKA_VARIANT_CALLER \
                and self.config["selection_method"] == SELECTION_METHOD_ALLELE_CALLS:
            message = "Non supported selection method [%s] for variant caller [%s] which has no genotypes" \
                      % (SELECTION_METHOD_ALLELE_CALLS, STRELKA_VARIANT_CALLER)
            logging.error(message)
            raise VcfDedupInputError(message)
        has_valid_sample_idx = "sample_idx" in self.config and self.config["sample_idx"] is not None
        if has_valid_sample_idx:
            try:
                sample_idx = int(self.config["sample_idx"])
            except ValueError:
                message = "'sample_idx' must be numeric and found [%s]" % str(sample_idx)
                logging.error(message)
                raise VcfDedupInputError(message)
        has_valid_sample_name = "sample_name" in self.config and self.config["sample_name"] is not None
        if self.config["variant_caller"] in [STRELKA_VARIANT_CALLER, STARLING_VARIANT_CALLER] \
                and not has_valid_sample_idx and not has_valid_sample_name:
            message = "Variant caller [%s] requires that either [sample_name] or [sample_idx] are provided" \
                      % self.config["variant_caller"]
            logging.error(message)
            raise VcfDedupInputError(message)
        if "sort_mem_percentage" in self.config:
            try:
                sort_mem_percentage = int(self.config["sort_mem_percentage"])
                if sort_mem_percentage <= 0 or sort_mem_percentage > 100:
                    message = "'sort_mem_percentage' must be in the range [1, 100] and it was found [%s]" % \
                              str(sort_mem_percentage)
                    logging.error(message)
                    raise VcfDedupInputError(message)
            except ValueError:
                message = "'sort_mem_percentage' must be numeric"
                logging.error(message)
                raise VcfDedupInputError(message)

    def __sort(self):
        logging.info("Sorts the VCF before processing...")
        input_basename = os.path.splitext(os.path.basename(self.input_vcf))[0]
        self.sorted_vcf = os.path.join(self.temp_folder, "." + input_basename + ".sorted.vcf.gz")
        vcf_sorter = VcfSorter(
            self.input_vcf,
            self.sorted_vcf,
            temp_folder=self.temp_folder,
            threads=self.sort_threads,
            mem_percentage=self.sort_mem_percentage
        )
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
                self.selection_method
            )
        elif self.variant_caller == GENERIC_VARIANT_CALLER:
            transformer = GenericVcfDedupper(
                self.sorted_vcf,
                self.output_vcf,
                comparer,
                self.selection_method
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

