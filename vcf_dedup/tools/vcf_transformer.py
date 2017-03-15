import vcf
import os.path
import logging
from abc import ABCMeta, abstractmethod
import time
import multiprocessing
import itertools


def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

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

        import copy_reg
        import types
        copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

        logging.info("Initialised!")

    def __del__(self):
        """
        Closes reader and writer.
        :return:
        """
        logging.info("Closing the VCF writer...")
        self.writer.flush()
        self.writer.close()
        logging.info("Resources closed!")

    def process_vcf(self):
        """
        PRE: the VCF is sorted so equal variants must be adjacent
        :return:
        """
        logging.info("Starts processing the VCF...")
        before = time.time()
        # TODO: determine from the configuration the number of threads to use
        #pool = multiprocessing.Pool()
        #pool.map(process_contig, itertools.izip(itertools.repeat(self), contigs_records))

        #FIXME: error reading contig chr11_KI270721v1_random
        #for contig in self.reader.contigs:
        #    self.process_contig(contig)
        count = 0
        for variant in self.reader:
            self.transform_variant(variant)
            count += 1
            if count % 10000 == 0:
                logging.info("Processed %s variants..." % str(count))

        after = time.time()
        logging.info("Finished processing %s variants in %s seconds" % (str(count), str(after - before)))


    def process_contig(self, contig):
        """

        :param contig: the contig to process
        :return:        None
        """
        # iterates all variants
        logging.info("Processing contig %s" % contig)
        before = time.time()
        for variant in self.reader.fetch(contig):
            self.transform_variant(variant)
        after = time.time()
        logging.info("Finished contig %s in %s seconds" % (contig, str(after - before)))

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
        # Reads previous variant if any
        prev_variant = None
        if len(self.variants) > 0:
            prev_variant = self.variants[len(self.variants) - 1]
        # if previous and current variants are not equal it merges the block and stores the current variant for the
        # next block
        if prev_variant is not None and \
                not self.variant_comparer.equals(prev_variant, variant):
            # writes merged variant
            self.writer.write_record(self.merge_variants(self.variants))
            # resets variants memory
            self.variants = [variant]
        # current variant forms part of the same block, stores and continues
        else:
            self.variants.append(variant)

    def merge_variants(self, variants):
        """
        Get all variants with PASS
        # If only 1, return that one
        # if more than 1, calculate AF over all PASS and get highest
        # if 0, calculate AF over all and get highest
        :param variants:
        :return: the merged variant
        """
        merged_variant = None
        # gets all passed variants
        passed_variants = [variant for variant in variants if len(variant.FILTER) == 0]
        if len(passed_variants) == 1:
            merged_variant = passed_variants[0]
        elif len(passed_variants) > 1:
            allele_frequencies = {passed_variant : self.calculate_somatic_AF(passed_variant)
                                  for passed_variant in passed_variants}
            merged_variant = max(allele_frequencies, key=allele_frequencies.get)
        elif len(passed_variants) == 0:
            allele_frequencies = {variant: self.calculate_somatic_AF(variant) for variant in variants}
            merged_variant = max(allele_frequencies, key=allele_frequencies.get)

        return merged_variant


    def calculate_somatic_AF(self, variant):
        """
        Return the ratio of supporting reads for somatic variants called by Strelka
        For SNVs: alternate AC / (alternate AC + reference AC)
        For indels: indel AC / (alternate AC + indel AC)
        :param variant:
        :return:
        """
        logging.debug("Calculating somatic ratio of supporting reads for %s:%s" % (variant.CHROM, str(variant.POS)))
        af = 0
        if variant.is_snp:
            reference = str(variant.REF)
            # TODO: capture error when multiallelic gets here
            alternate = str(variant.ALT[0])
            # TODO: capture error when no samples available
            # TODO: consider case when multiple samples are available
            # TODO: capture error when no AU, CU, etc in sample data
            reference_ac = variant.samples[0].data._asdict()[reference + "U"][0]
            alternate_ac = variant.samples[0].data._asdict()[alternate + "U"][0]
            af = alternate_ac / (alternate_ac + reference_ac)
        elif variant.is_indel:
            # TODO: capture error when no TIR or TAR in sample data
            indel_ac = variant.samples[0].data._asdict()["TIR"][0]
            alternate_ac = variant.samples[0].data._asdict()["TAR"][0]
            af = indel_ac / (alternate_ac + indel_ac)
        elif variant.is_sv:
            # TODO: what should we do with these ones?
            pass
        else:
            # TODO: will this ever happen?
            pass

        return af








