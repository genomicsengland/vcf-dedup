import vcf
import vcf.model
import os.path
import logging
from abc import ABCMeta, abstractmethod
import time
import sys
from collections import OrderedDict
from vcf_dedup.constants import *


class VcfFormatError(Exception):
    """
    A exception to raise when VCF format errors are found
    """
    pass

class VcfDedupperError(Exception):
    """
    An unexpected exception in VCF dedupper
    """
    pass


class AbstractVcfTransformer(object):

    __metaclass__ = ABCMeta
    variants_buffer = {}
    BUFFER_SIZE = 10000

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
            if output_vcf_file is None or output_vcf_file == "":
                output_vcf = sys.stdout
            else:
                output_vcf = open(output_vcf_file, 'w')
            self.writer = vcf.VCFWriter(output_vcf, self.reader)
        except Exception, e:
            logging.error("Error opening output VCF file: " + str(e))
            raise ValueError("Error opening output VCF file: " + str(e))
        # loads writer for duplicated variants
        try:
            if output_vcf_file is None or output_vcf_file == "":
                output_duplicated_vcf = sys.stderr
            else:
                duplicated_vcf_file = os.path.join(
                    os.path.dirname(os.path.realpath(output_vcf_file)),
                    os.path.splitext(os.path.basename(output_vcf_file))[0] + ".duplicated.vcf"
                )
                output_duplicated_vcf = open(duplicated_vcf_file, 'w')
            self.writer_duplicated = vcf.VCFWriter(output_duplicated_vcf, self.reader)
        except Exception, e:
            logging.error("Error opening output VCF file: " + str(e))
            raise ValueError("Error opening output VCF file: " + str(e))

        logging.info("Initialised!")

    def __del__(self):
        """
        Closes reader and writer.
        :return:
        """
        logging.info("Closing the VCF writer...")
        self.writer.flush()
        self.writer.close()
        self.writer_duplicated.flush()
        self.writer_duplicated.close()
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

        # FIXME: error reading contig chr11_KI270721v1_random
        #for contig in self.reader.contigs:
        #    self.process_contig(contig)
        count = 0
        for variant in self.reader:
            if len(variant.ALT) > 1:
                logging.error("The VCF contains multi-allelic variants %s:%s%s>%s. Please split them before processing!" % (variant.CHROM, str(variant.POS), str(variant.REF), str(variant.ALT)))
                raise VcfFormatError("The VCF contains multi-allelic variants. Please split them before processing!")
            self._transform_variant(variant)
            count += 1
            if count % 10000 == 0:
                logging.info("Processed %s variants..." % str(count))
        # this methods should take care of writing the last group of variants and flushing the buffers
        self._write_last_block()
        after = time.time()
        logging.info("Finished processing %s variants in %s seconds" % (str(count), str(after - before)))

    def _write_variant(self, writer, variant):
        """
        Writes variants in batches of a considerable size
        :param writer:
        :param variant:
        :return:
        """
        if writer in self.variants_buffer:
            self.variants_buffer[writer].append(self._variant2dict(writer, variant))
            # writes a batch of variants when reached the buffer size
            if len(self.variants_buffer[writer]) >= self.BUFFER_SIZE:
                self._flush_buffer(writer)
                self.variants_buffer[writer] = []
        else:
            self.variants_buffer[writer] = [self._variant2dict(writer, variant)]

    def _flush_buffer(self, writer):
        """
        Flushes the buffer for a given writer
        :param writer:
        :return:
        """
        if writer in self.variants_buffer:
            logging.info("Writing %s variants" % str(len(self.variants_buffer[writer])))
            writer.writer.writerows(self.variants_buffer[writer])
            self.variants_buffer[writer] = []
        else:
            error_message = "Trying to flush a writer with no buffer"
            logging.warning(error_message)

    def _variant2dict(self, writer, variant):
        """
        Transforms a variant object into a dictionary ready to be written
        :param variant:
        :return:
        """
        ffs = writer._map(str, [variant.CHROM, variant.POS, variant.ID, variant.REF])
        ffs += [
            writer._format_alt(variant.ALT),
            variant.QUAL or '.',
            writer._format_filter(variant.FILTER),
            writer._format_info(variant.INFO)
        ]
        if variant.FORMAT:
            ffs.append(variant.FORMAT)
        samples = [writer._format_sample(variant.FORMAT, sample)
                   for sample in variant.samples]
        return ffs + samples

    def _process_contig(self, contig):
        """

        :param contig: the contig to process
        :return:        None
        """
        # iterates all variants
        logging.info("Processing contig %s" % contig)
        before = time.time()
        for variant in self.reader.fetch(contig):
            self._transform_variant(variant)
        after = time.time()
        logging.info("Finished contig %s in %s seconds" % (contig, str(after - before)))

    @abstractmethod
    def _transform_variant(self, variant):
        pass

    def _write_last_block(self):
        pass




class DuplicationFinder(AbstractVcfTransformer):

    variants = []

    def __init__(self, input_vcf_file, output_vcf_file, variant_comparer):
        """

        :param input_vcf_file:
        :param output_vcf_file:
        :param variant_comparer:
        """
        # calls the parent constructor
        AbstractVcfTransformer.__init__(self, input_vcf_file, output_vcf_file)
        # stores the variant comparer
        self.variant_comparer = variant_comparer

    def _transform_variant(self, variant):
        """
        Variant transformation. Stores all variants with equal coordinates and writes duplicated variants.
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
            # writes variants only duplicated variants
            if len(self.variants) > 1:
                for duplicated_variant in self.variants:
                    self._write_variant(self.writer, duplicated_variant)
            # resets variants memory
            self.variants = [variant]
        # current variant forms part of the same block, stores and continues
        else:
            self.variants.append(variant)

    def _write_last_block(self):
        """
        Writes the last block stored in self.variants
        :return:
        """
        if len(self.variants) > 1:
            # writes merged variant
            for variant in self.variants:
                self._write_variant(self.writer, variant)
        # clears it
        self.variants = []
        # flush the buffer
        self._flush_buffer(self.writer)


class AbstractVcfDedupper(AbstractVcfTransformer):

    variants = []

    def __init__(self,
                 input_vcf_file,
                 output_vcf_file,
                 variant_comparer,
                 selection_method,
                 sample_idx = None,
                 sample_name = None
                 ):

        # calls the parent constructor
        AbstractVcfTransformer.__init__(self, input_vcf_file, output_vcf_file)
        # stores the variant comparer
        self.variant_comparer = variant_comparer
        # sets the selection method
        if selection_method == SELECTION_METHOD_AF:
            self._select_variants = self._select_mode_af
        elif selection_method == SELECTION_METHOD_QUALITY:
            self._select_variants = self._select_mode_quality
        elif selection_method == SELECTION_METHOD_ALLELE_CALLS:
            self._select_variants = self._select_mode_allele_calls
        elif selection_method == SELECTION_METHOD_ARBITRARY:
            self._select_variants = self._select_mode_arbitrary
        self.sample_idx = sample_idx
        self.sample_name = sample_name
        self.duplications_count = 0

    def _transform_variant(self, variant):
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
            if len(self.variants) > 1:
                # writes merged variant
                self._write_variant(self.writer, self._merge_variants(self.variants))
                # writes to VCF of duplicated variants
                for duplicated_variant in self.variants:
                    self._write_variant(self.writer_duplicated, duplicated_variant)
            else:
                # writes the variant when there is no duplication
                self._write_variant(self.writer, self.variants[0])
            # resets variants memory
            self.variants = [variant]
        # current variant forms part of the same block, stores and continues
        else:
            self.variants.append(variant)

    @abstractmethod
    def _calculate_AF(self, variant):
        pass

    @abstractmethod
    def _get_variant_calling_quality(self, variant):
        pass

    def _select_mode_allele_calls(self, variants):
        """
        Returns the variants with the greater number of allele calls. If there is a collision returns the maximum AF.
        Otherwise returns an arbitrary variant.
        :param variants:
        :return:
        """
        allele_calls = OrderedDict()
        for variant in variants:
            allele_calls[variant] = self._calculate_allele_calls(variant)
        max_ac = allele_calls[max(allele_calls, key=allele_calls.get)]
        max_acs = [i for i, x in allele_calls.iteritems() if x == max_ac]
        if (len(max_acs) > 1):
            # collision with maximum allele call
            afs = OrderedDict()
            for variant in variants:
                afs[variant] = self._calculate_AF(variant)
            max_af = afs[max(afs, key=afs.get)]
            max_afs = [i for i, x in afs.iteritems() if x == max_af]
            if (len(max_afs) > 1):
                # collision with maximum AF
                qualities = OrderedDict()
                for variant in variants:
                    qualities[variant] = self._get_variant_calling_quality(variant)
                max_qual = qualities[max(qualities, key=qualities.get)]
                max_quals = [i for i, x in qualities.iteritems() if x == max_qual]
                # returns highest quality
                merged_variant = max_quals[0]
            else:
                # returns highest AFs
                merged_variant = max_afs[0]
        else:
            # returns highest allele calls
            merged_variant = max_acs[0]

        return merged_variant

    def _select_mode_af(self, variants):
        """
        Returns the variants with the highest allele frequency. If there is a collision returns the variants with
        the maximum quality. Otherwise returns an arbitrary variant.
        :param variants:
        :return:
        """
        allele_frequencies = OrderedDict()
        for variant in variants:
            allele_frequencies[variant] = self._calculate_AF(variant)
        max_af = allele_frequencies[max(allele_frequencies, key=allele_frequencies.get)]
        max_afs = [i for i, x in allele_frequencies.iteritems() if x == max_af]
        if (len(max_afs) > 1):
            # collision with maximum AF
            qualities = OrderedDict()
            for variant in variants:
                qualities[variant] = self._get_variant_calling_quality(variant)
            max_qual = qualities[max(qualities, key=qualities.get)]
            max_quals = [i for i, x in qualities.iteritems() if x == max_qual]
            # gets the first
            merged_variant = max_quals[0]
        else:
            merged_variant = max_afs[0]

        return merged_variant

    def _select_mode_quality(self, variants):
        """
        Returns the variant with the highest variant calling quality. If there is a collision returns the maximum AF.
        Otherwise returns an arbitrary variant.
        :param variants:
        :return:
        """
        qualities = OrderedDict()
        for variant in variants:
            qualities[variant] = self._get_variant_calling_quality(variant)
        max_qual = qualities[max(qualities, key=qualities.get)]
        max_quals = [i for i, x in qualities.iteritems() if x == max_qual]
        if (len(max_quals) > 1):
            # collision with maximum AF
            allele_frequencies = OrderedDict()
            for variant in variants:
                allele_frequencies[variant] = self._calculate_AF(variant)
            max_af = allele_frequencies[max(allele_frequencies, key=allele_frequencies.get)]
            max_afs = [i for i, x in allele_frequencies.iteritems() if x == max_af]
            # gets the first
            merged_variant = max_afs[0]
        else:
            merged_variant = max_quals[0]
        return merged_variant

    def _select_mode_arbitrary(self, variants):
        """
        Returns the first variant from the list
        :param variants:
        :return:
        """
        merged_variant = variants[0]
        return merged_variant

    def _merge_variants(self, variants):
        """
        Get all variants with PASS
        # If only 1, return that one
        # if more than 1, run selection mode on passed variants
        # if 0, calculate AF run selection mode on all
        :param variants:
        :return: the merged variant
        """
        logging.info("Merging duplicate variants at %s:%s:%s>%s" % (
            variants[0].CHROM,
            int(variants[0].POS),
            variants[0].REF,
            variants[0].ALT
        ))
        merged_variant = None
        # gets all passed variants
        passed_variants = [variant for variant in variants if len(variant.FILTER) == 0]
        if len(passed_variants) == 1:
            merged_variant = passed_variants[0]
        elif len(passed_variants) > 1:
            merged_variant = self._select_variants(passed_variants)
        elif len(passed_variants) == 0:
            merged_variant = self._select_variants(variants)
        self.duplications_count += 1

        return merged_variant

    def _get_sample(self, variant):
        """
        Returns the sample identified by index or sample name.
        :param variant:         the variant
        :param sample_idx:      the index
        :param sample_name:     the sample name
        :return:                the selected sample
        """
        if self.sample_name is not None and self.sample_name != "":
            samples = [sample for sample in variant.samples if sample.sample == self.sample_name]
            if len(samples) == 0:
                raise VcfFormatError("Sample '%s' is not available in the VCF" % self.sample_name)
            sample = samples[0]
        else:
            sample = variant.samples[self.sample_idx]
        return sample

    def _write_last_block(self):
        """
        Writes the last block stored in self.variants
        :return:
        """
        if len(self.variants) > 1:
            # writes merged variant
            self._write_variant(self.writer, self._merge_variants(self.variants))
            # writes duplictade variants
            for duplicated_variant in self.variants:
                self._write_variant(self.writer_duplicated, duplicated_variant)
        else:
            # writes the variant when there is no duplication
            self._write_variant(self.writer, self.variants[0])
        # clears it
        self.variants = []
        logging.info("Found %s positions having duplicated variants" % (str(self.duplications_count)))
        self._flush_buffer(self.writer)
        self._flush_buffer(self.writer_duplicated)

    def _calculate_allele_calls(self, variant):
        """
        Return the ratio of allele calls for the alternate allele.
        :param variant:
        :return:
        """
        logging.debug("Calculating ratio of alternate genotypes for %s:%s" % (variant.CHROM, str(variant.POS)))
        alternate_alleles = 0
        total_alleles = 0
        for sample in variant.samples:
            alternate_alleles += sample.gt_alleles.count('1')
            total_alleles += len(sample.gt_alleles)
        if total_alleles == 0:
            # when no samples are detected in the VCF we don't want to fail with a division by zero
            return 0
        else:
            return float(alternate_alleles) / total_alleles


class StrelkaVcfDedupper(AbstractVcfDedupper):
    """
    Strelka is the variant caller for somatic variants in cancer program.
    This class performs the merging of duplicated variants from Strelka.
    PRE: VCFs are single sample
    """

    def _calculate_AF(self, variant):
        """
        Return the ratio of supporting reads for somatic variants called by Strelka
        For SNVs: alternate AC / (alternate AC + reference AC)
        For indels: indel AC / (alternate AC + indel AC)
        For SNVs if the fields reference + "U" or alternate + "U" do not exist in format returns 0
        For indels if the fields "TIR" or "TAR" do not exist in format returns 0
        :param variant:
        :return:
        """
        logging.debug("Calculating somatic ratio of supporting reads for %s:%s" % (variant.CHROM, str(variant.POS)))
        af = 0
        format = self._get_sample(variant).data._asdict()
        if variant.is_snp:
            reference = str(variant.REF)
            # TODO: capture error when multiallelic gets here
            alternate = str(variant.ALT[0])
            if reference + "U" in format and alternate + "U" in format:
            # TODO: capture error when no samples available
            # TODO: consider case when multiple samples are available
                reference_ac = int(format[reference + "U"][0])
                alternate_ac = int(format[alternate + "U"][0])
                af = alternate_ac / (alternate_ac + reference_ac) if alternate_ac + reference_ac > 0 else 0
        elif variant.is_indel:
            if ("TIR" in format and "TAR" in format):
                indel_ac = int(format["TIR"][0])
                alternate_ac = int(format["TAR"][0])
                af = float(indel_ac) / (alternate_ac + indel_ac) if alternate_ac + indel_ac > 0 else 0
        elif variant.is_sv:
            # TODO: what should we do with these ones?
            pass
        else:
            # TODO: will this ever happen?
            pass
        return af

    def _get_variant_calling_quality(self, variant):
        """
        Retrieves the variant calling quality from the VQSR info field
        'Recalibrated quality score expressing the phred scaled probability of the somatic call being a FP observation.'
        :param variant:     the variant
        :return:            the variant calling quality
        """
        variant_calling_quality = 0.0
        if "VQSR" in variant.INFO:
            variant_calling_quality = variant.INFO["VQSR"]
        try:
            variant_calling_quality = float(variant_calling_quality)
        except Exception:
            variant_calling_quality = 0.0
        return variant_calling_quality



class StarlingVcfDedupper(AbstractVcfDedupper):
    """
    Starling is the variant caller for cancer germline variants.
    This class performs the merging of duplicated variants from Starling.
    PRE: VCFs are single sample
    """

    def _calculate_AF(self, variant):
        """
        Return the ratio of supporting reads for variants called by Starling
        For SNVs and indels: alternate AC / (alternate AC + reference AC)
        If the field AC is not present it will set AF to 0
        :param variant:
        :return:
        """
        logging.debug("Calculating ratio of supporting reads for %s:%s" % (variant.CHROM, str(variant.POS)))
        af = 0
        format = self._get_sample(variant).data._asdict()
        if ("AD" in format and len(format["AD"]) == 2):
            reference_ac = int(format["AD"][0])
            alternate_ac = int(format["AD"][1])
            af = float(alternate_ac) / (alternate_ac + reference_ac) if alternate_ac + reference_ac > 0 else 0
        return af

    def _get_variant_calling_quality(self, variant):
        """
        Retrieves the variant calling quality from the GQX format field
        'Empirically calibrated variant quality score for variant sites, otherwise Minimum of {Genotype quality
        assuming variant position,Genotype quality assuming non-variant position}'
        :param variant:     the variant
        :return:            the variant calling quality
        """
        variant_calling_quality = 0
        format = self._get_sample(variant).data._asdict()
        if "GQX" in format:
            variant_calling_quality = format["GQX"]
        try:
            variant_calling_quality = int(variant_calling_quality)
        except Exception:
            variant_calling_quality = 0
        return variant_calling_quality


class PlatypusVcfDedupper(AbstractVcfDedupper):
    """
    Platypus is the variant caller for rare disease germline variants.
    This class performs the merging of duplicated variants from Platypus.
    PRE: VCFs are multi-sample
    """

    def _calculate_AF(self, variant):
        """
        Return the ratio of supporting reads for variants called by Platypus
        For SNVs and indels: TR / TC
        :param variant:
        :return:
        """
        logging.debug("Calculating ratio of supporting reads for %s:%s" % (variant.CHROM, str(variant.POS)))
        af = 0
        if ("TC" in variant.INFO and "TR" in variant.INFO):
            total_ac = int(variant.INFO["TC"])
            # NOTE: we assume that multi-allelic variants do not reach this stage in normalisation pipeline
            alternate_ac = int(variant.INFO["TR"][0])
            af = float(alternate_ac) / (total_ac)
        return af

    def _get_variant_calling_quality(self, variant):
        """
        Retrieves the variant calling quality from the QUAL field
        :param variant:     the variant
        :return:            the variant calling quality
        """
        variant_calling_quality = variant.QUAL if variant.QUAL is not None else 0
        return variant_calling_quality


class GenericVcfDedupper(AbstractVcfDedupper):
    """
    This generic VCF dedupper is aimed to support any VCF based on standard annotations in the VCF specification.
    """

    def _calculate_AF(self, variant):
        """
        Return the allele frequency as annotated in AF field
        :param variant:
        :return:
        """
        af = 0
        if "AF" in variant.INFO:
            af = variant.INFO["AF"]
        return af

    def _get_variant_calling_quality(self, variant):
        """
        Retrieves the variant calling quality from the QUAL field
        :param variant:     the variant
        :return:            the variant calling quality
        """
        variant_calling_quality = variant.QUAL if variant.QUAL is not None else 0
        return variant_calling_quality





