import logging
import os
import sys
from subprocess import call


class VcfSorter(object):

    GREP = "grep"
    ZGREP = "zgrep"

    def __init__(self, input_vcf_file, output_vcf_file, temp_folder = None, threads = 1):
        """
        Constructtor for the VCF sorter
        :param input_vcf_file:
        :param output_vcf_file:
        """
        logging.info("Initialising VCF sorter...")
        logging.info("Input VCF file: %s" % input_vcf_file)
        logging.info("Output VCF file: %s" % output_vcf_file)
        if not os.path.exists(input_vcf_file):
            logging.error("Input VCF file does not exist!")
            raise ValueError("Input VCF file does not exist!")
        if output_vcf_file is None or output_vcf_file == "":
            logging.error("Output VCF file is not provided!")
            raise ValueError("Output VCF file is not provided!")
        self.input_vcf = input_vcf_file
        if self.input_vcf.endswith(".gz"):
            self.is_compressed = True
        else:
            self.is_compressed = False
        self.output_vcf = output_vcf_file
        self.threads = threads
        logging.info("Running sort with %s threads" % str(threads))
        if temp_folder is None or temp_folder == "":
            self.temp_folder = os.path.dirname(os.path.realpath(self.output_vcf))
        else:
            self.temp_folder = temp_folder
        logging.info("Temporary folder [%s]" % temp_folder)
        logging.info("Initialised!")

    def __sort(self):
        sort_command = "{0} -v '#' {1} | sort {2} {3} -k1,1V -k2,2n -k4,4V -k5,5V >> {4}".format(
            self.ZGREP if self.is_compressed else self.GREP,
            self.input_vcf,
            "--parallel=%s" % str(self.threads) if self.threads > 1 else "",
            "--temporary-directory=%s" % self.temp_folder,
            self.output_vcf)
        logging.info("Running [%s]" % sort_command)
        os.system(sort_command)

    def __write_header(self):
        get_header_command = "{0} '#' {1} > {2}".format(
            self.ZGREP if self.is_compressed else self.GREP,
            self.input_vcf,
            self.output_vcf)
        logging.info("Running [%s]" % get_header_command)
        os.system(get_header_command)

    def sort(self):
        """
        Runs the VCF sorting
        :return:
        """
        logging.info("Sorting...")
        self.__write_header()
        self.__sort()
        logging.info("Sorted!")