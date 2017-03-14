import vcf


class VcfTransformer(object):


    def __init__(self, input_vcf_file, output_vcf_file):

        # loads reader
        self.input_vcf_file = input_vcf_file
        self.reader = vcf.VCFReader(filename = self.input_vcf_file)
        # loads writer
        self.output_vcf_file = output_vcf_file
        self.writer = vcf.VCFWriter(open(self.output_vcf_file, 'w'), self.reader)

    def iterate_vcf(self):
        """
        PRE: the VCF is sorted so equal variants must be adjacent
        :return:
        """

        variants = []
        prev_variant = None
        # iterates all variants
        for variant in self.reader:
            if prev_variant is not None and not self.equal_variants(prev_variant, variant):
                # transforms variants
                self.transform_variants(variants)
                # initialize variants memory
                variants = [variant]
            else:
                variants.append(variant)
            prev_variant = variant

    def transform_variants(self, variants):

        # TODO: do something to the list variants
        transformed_variant = None
        # merges all variants in one, takes the first one not filtered or just the first one
        for variant in variants:
            if transformed_variant is None:
                transformed_variant = variant
            else:
                if len(variant.FILTER) == 0:
                    transformed_variant.FILTER.append("PASS")
                else:
                    transformed_variant.FILTER = transformed_variant.FILTER + variant.FILTER

            #elif transformed_variant.is_filtered and not variant.is_filtered:
            #    transformed_variant = variant
        # writes transformed variant
        self.writer.write_record(transformed_variant)



    def equal_variants(self, variant1, variant2):
        """

        :param variant1:
        :param variant2:
        :return:  a boolean indicating if the variants are equal
        """
        #TODO: implement something
        return (variant1.CHROM == variant2.CHROM and variant1.POS == variant2.POS and
                variant1.REF == variant2.REF)
        # TODO: use also ALT, beware that this is a list
        ## and variant1.ALT == variant2.ALT)

