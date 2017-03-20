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
