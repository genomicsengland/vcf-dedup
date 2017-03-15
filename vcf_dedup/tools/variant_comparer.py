from abc import ABCMeta, abstractmethod


class AbstractVariantComparer(object):

    __metaclass__ = ABCMeta

    @abstractmethod
    def equals(self, variant1, variant2):
        pass


class VariantComparerNoAlternate(AbstractVariantComparer):

    def equals(self, variant1, variant2):

        return (variant1.CHROM == variant2.CHROM and
                variant1.POS == variant2.POS and
                variant1.REF == variant2.REF)


class VariantComparerWithAlternate(AbstractVariantComparer):

    def equals(self, variant1, variant2):

        return (variant1.CHROM == variant2.CHROM and
                variant1.POS == variant2.POS and
                variant1.REF == variant2.REF and
                set(variant1.ALT) == set(variant2.ALT))