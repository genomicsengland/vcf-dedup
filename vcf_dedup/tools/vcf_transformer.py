import vcf
import os.path
import logging
from abc import ABCMeta, abstractmethod


class AbstractVcfTransformer(object):

    __metaclass__ = ABCMeta

    def __init__(self, input_vcf_file, output_vcf_file):
        """
        Loads VCF reader and writer.
        :param input_vcf_file: the input file
        :param output_vcf_file: the output file
        """

        # loads reader
        logging.info("Initialising VCF transformer...")
        logging.info("Input VCF file: %s" % input_vcf_file)
        logging.info("Output VCF file: %s" % output_vcf_file)
        if os.path.exists(input_vcf_file):
            try:
                self.input_vcf_file = input_vcf_file
                self.reader = vcf.VCFReader(filename = self.input_vcf_file)
            except Exception, e:
                logging.error("Error opening input VCF file: " + str(e))
                raise ValueError("Error opening input VCF file: " + str(e))
        else:
            logging.error("Input VCF file does not exist!")
            raise ValueError("Input VCF file does not exist!")
        # loads writer
        try:
            self.output_vcf_file = output_vcf_file
            self.writer = vcf.VCFWriter(open(self.output_vcf_file, 'w'), self.reader)
        except Exception, e:
            logging.error("Error opening output VCF file: " + str(e))
            raise ValueError("Error opening output VCF file: " + str(e))
        logging.info("Initialised!")

    def __del__(self):
        """
        Closes reader and writer.
        :return:
        """

        logging.info("Closing the VCF reader and writer...")
        self.reader.close()
        self.writer.close()
        logging.info("Resources closed!")

    def process_vcf(self):
        """
        PRE: the VCF is sorted so equal variants must be adjacent
        :return:
        """

        logging.info("Starts processing the VCF...")
        # iterates all variants
        for variant in self.reader:
            self.transform_variant(variant)
        logging.info("Finished processing the VCF!")


    @abstractmethod
    def transform_variant(self, variant):
        pass


class VcfDedupper(AbstractVcfTransformer):

    variants = []

    def __init__(self, input_vcf_file, output_vcf_file, variant_comparer):

        # calls the parent constructor
        AbstractVcfTransformer.__init__(self, input_vcf_file, output_vcf_file)
        # stores the variant comparer
        self.variant_comparer = variant_comparer

    def transform_variant(self, variant):
        """
        Variant transformation. Stores all variants with equal coordinates and writes a merged result.
        :param variant: the current variant
        :return: None
        """

        # TODO: do something to the list variants
        # Reads previous variant
        prev_variant = None
        if len(self.variants) > 0:
            prev_variant = self.variants[len(self.variants) - 1]
        if prev_variant is not None and \
                not self.variant_comparer.equals(prev_variant, variant):
            # transforms variants
            merged_variant = None
            is_passed = False
            # merges all variants in one, takes the first one not filtered or just the first one
            for variant in self.variants:
                if merged_variant is None:
                    merged_variant = variant
                    is_passed = len(merged_variant.FILTER) == 0
                elif not is_passed:
                    if len(variant.FILTER) == 0:
                        merged_variant = variant
            # writes transformed variant
            self.writer.write_record(merged_variant)
            # resets variants memory
            self.variants = [variant]
        else:
            self.variants.append(variant)








