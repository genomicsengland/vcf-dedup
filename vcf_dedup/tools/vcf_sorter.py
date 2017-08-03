import logging
import os
import subprocess


class VcfSortError(Exception):
    """
    A exception to raise when an error sorting happens
    """
    pass


class VcfSorter(object):

    GREP = "grep"
    ZGREP = "zgrep"

    def __init__(self, input_vcf_file, output_vcf_file, temp_folder=None, threads=1, mem_percentage=80):
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
        if threads is not None and threads > 0:
            self.threads = threads
        else:
            self.threads = 1
        logging.info("Running sort with %s threads" % str(self.threads))
        if temp_folder is None or temp_folder == "":
            self.temp_folder = os.path.dirname(os.path.realpath(self.output_vcf))
        else:
            self.temp_folder = temp_folder
        logging.info("Temporary folder [%s]" % self.temp_folder)
        if mem_percentage is not None and mem_percentage > 0 and mem_percentage <= 100:
            self.mem_percentage = mem_percentage
        else:
            self.mem_percentage = 80
        logging.info("Memory percentage to be used by sort [%s]" % str(self.mem_percentage))
        logging.info("Initialised!")

    def __sort(self):
        '''
        Pastes the header, sorts the variants by chromosome, position, reference and alternate and compress them
        :return:
        '''
        sort_command = "({0} '#' {1} && {0} -v '#' {1} | sort {2} {3} {4} -k1,1V -k2,2n -k4,4V -k5,5V) | bgzip > {5}"\
            .format(
            self.ZGREP if self.is_compressed else self.GREP,
            self.input_vcf,
            "--parallel=%s" % str(self.threads) if self.threads > 1 else "",
            "--temporary-directory=%s" % self.temp_folder,
            "--buffer-size=%s%%" % str(self.mem_percentage) if self.mem_percentage is not None else "",
            self.output_vcf
        )
        logging.info("Running [%s]" % sort_command)
        sort = subprocess.Popen(sort_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sort.communicate()
        # raise an error if sort return code is other than 0
        if sort.returncode:
            error_message = 'Sort returned error code {0}'.format(sort.returncode)
            logging.error(error_message)
            raise VcfSortError(error_message)

    def sort(self):
        """
        Runs the VCF sorting
        :return:
        """
        logging.info("Sorting...")
        self.__sort()
        logging.info("Sorted!")