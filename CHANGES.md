version 0.4.1 (11 Jul 2017)
----------------------------

* Minor changes
    - The output of sort operation is now bgzipped to better process big VCFs

* Known issues
    - The sort temporary folder parameter is broken

version 0.4.0 (6 Jul 2017)
----------------------------

* Major changes
    - Adding support for a generic mode that relies on standard VCF annotations. QUAL column for variant calling quality and AF INFO annotation for the allele frequency.
    - Platypus support now uses QUAL instead of QD as variant calling quality value.
* Minor changes
    - Options not used for Platypus and generic modes `--sample-name` and `--sample-idx` are not required now. They are still required for Strelka and Starling.

version 0.3.1 (5 Jul 2017)
----------------------------

* Minor changes
    - Making sure the writers are flushed and closed.


version 0.3.0 (12 May 2017)
----------------------------

* Major Changes:
    - Option to sort the input VCF by chromosome, position, reference and alternate bases
    - Platypus is now supported
    - Additional selection method `allele_calls`. This selection method is not available as it does not report genotypes.

* Minor changes
    - Duplicated variants are now written into file `${output_vcf}.duplicated.vcf`

version 0.2.2 (22 March 2017)
----------------------------

* Minor Changes:
    - `duplication_finder` failed to write properly last set of duplicates.
    - Fixed bug in calculation of allele frequencies for Starling. We were using AC, which did not exist, instead of AD.

version 0.2.1 (22 March 2017)
----------------------------

* Minor Changes:
    - `duplication_finder` failed to write last set of duplicates.

version 0.2.0 (20 March 2017)
----------------------------

* Major Changes:
    - `sample-name` parameter added in addition to `sample-idx`
    - When output is not provdied it writes to standard output
    - Fixed bug when calculating allele frequencies

version 0.1.5 (20 March 2017)
----------------------------

* Minor Changes:
    - `sample-idx` must be forced to numeric

version 0.1.4 (20 March 2017)
----------------------------

* Minor Changes:
    - Bugfix in command line parameter `sample-idx`

version 0.1.3 (20 March 2017)
----------------------------

* Minor Changes:
    - Fixed bug on script `vcf-dedupper`

version 0.1.1 (18 March 2017)
----------------------------

* Minor Changes:
    - Fixed bug on `sample-idx` parameter


version 0.1.0 (17 March 2017)
----------------------------

* Major Changes:
    - First implementation
    - Removal of duplicates for VCFs from Strelka and Starling
    - Duplicates are removed based on the filtering status, passed variant is always selected
    - Ties in filtering status are resolved by using alternatively the maximum allele frequency, the maximum variant calling quality or either arbitrarily
    - No consideration for multisample VCFs. The parameter `sample-idx` indicates the sample to analyse.
