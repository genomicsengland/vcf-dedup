# VCF dedupper

VCF files may contain duplicated after normalisation or merging of VCFs. This module intends to remove those duplicates in a sensible way. This process is very specific to the variant caller, we only support Strelka, Starling and Platypus variant callers. The variant caller might be selected with the parameter `--variant-caller`.

Variant duplicates are identified by having the same genomic coordinates and optionally the same alternate allele. This is controlled with the parameter `--equality-mode`.

The currrent implementation takes into account:
1. The filtering information. Whenever there is one and only one variant non filtered (i.e.: PASS) This will be the resulting variant.
2. Other criteria are used when there is a collision in filtering status (i.e.: two or more non filtered variants or all filtered variants). This is set with `--selection-method`
    * Highest allele frequency (i.e.: ratio of supporting reads for the variant call)
      - if there is a tie, it selects the variant with the highest variant calling quality
      - if there is a tie again, it selects an arbitrary variant
    * Highest variant calling quality
      - if there is a tie, it selects the variant with the highest allele frequency
      - if there is a tie again, it selects an arbitrary variant
    * Greater number of allele calls
      - if there is a tie, it selects the variant with the highest allele frequency
      - if there is a tie again, it selects the variant with the highest variant calling quality
      - if there is a tie again, it selects an arbitrary variant
    
The output VCF should not contain any duplicate.

## Support for different variant callers

The number of allele calls is calculated equally for all supported variant callers based on the genotypes called.

### Strelka 

* Allele frequency
   - For SNVs: `FORMAT[alternate_base + "U"] / (FORMAT[alternate_base + "U"] + FORMAT[reference_base + "U"])`
   - For indels: `FORMAT["TIR"] / (FORMAT["TIR"] + FORMAT["TAR"])`
* Variant calling quality: `INFO["VQSR"]`

### Starling

* Allele frequency: `format["AD"][1] / (format["AD"][1] + format["AD"][0])`
* Variant calling quality: `format["GQX"]`

### Platypus

* Allele frequency: `INFO["TR"][0] / INFO["TC"]`
* Variant calling quality: `QUAL`

### Generic

* Allele frequency: `INFO["AF"]`
* Variant calling quality: `QUAL`


## Assumptions

* Multi-allelic variants are not supported
* VCF must be sorted by chromosome, position, reference and alternate bases (beware that some tools do not take into account reference and alternate bases for sorting). Otherwise, use the `--sort` option and the input file will be sorted before deduplicating.

## Requirements

* UNIX sort with the parallelization feature released with coreutils 8.6 (2010-10-15) is used.
* bgzip (samtools >= 1.2)

## Using the VCF dedupper

The VCF dedupper comes in two flavors:
1. A command line script
2. A python module that can be used programmatically

### Script

```
vcf_dedupper --input-vcf $INPUT_VCF --output-vcf $OUTPUT_VCF --variant-caller $(one of "generic", "strelka", "starling", "platypus" or "duplication_finder") --selection-method $(one of "af", "quality", "allele_calls", "arbitrary")
```

The `duplication_finder` mode writes to the output only those variants that are duplicated without removing any duplicates. This mode is intended to facilitate the analysis of the root cause of duplications.

The `generic` mode relies only on annotations defined in the VCF specification. This mode is intended to be used with any type of valid VCF.

The `equality-mode` may be changed to use a less restrictive definition of duplicates using only chromosome, position and reference. By default it uses chromosome, position, reference and alternate.

The `sample-idx` indicates as a 0-based index for multisample VCFs the sample to which collision criteria will be applied. Beware that for Strelka VCFs the tumor sample is in the second sample (i.e.: index = 1), this might change.

The `sample-name` indicates the name of the sample to which collision criteria will be applied. This parameter overrides `sample-idx`. If the sample name does not exist an error is raised.

Sorting the input vcf can be enabled with flag `--sort`. Sort uses the UNIX sort, it can be parallelized using the parameter `--sort-threads` and while the default temporary folder is that of the output VCF this can be customised using `--sort-temp-folder`. The maximum amount of memory used by sort is set to 80% of the available memory by default, but this can be customised using `--sort-mem-percentage`.

### Python module

```
from vcf_dedup.runner import VcfDedupRunner


config = {
        "input_vcf" : input_vcf,
        "output_vcf" : output_vcf,
        "variant_caller": variant_caller,
        "selection_method": selection_method,
        "equality_mode": equality_mode,
        "sample_idx": sample_idx,
        "sort_vcf": False,
        "sort_threads": 1,
        "temp_folder": "",
        "sort_mem_percentage": 90,
        "verbose": True
        "log_file": "vcf-dedup.log"
    }
# Calls the VCF dedupper
runner = VcfDedupRunner(config)
runner.process_vcf()
```

### Output

* If the parameter `output-vcf` is provided this file will be created.
* If the parameter `output-vcf` is not provided the resulting VCF will be written to the standard output.
* The duplicated variants will be written in `${output-vcf}.duplicated.vcf` if `output-vcf` is provided, otherwise they will be written to the standard error.
* If the parameter `--verbose` is provided logs will be written to the standard error or to the log file (`--log-file`), if provided. Beware that logs and the duplicated variants might be both written to the standard error.

### Sample commands

For strelka it is agreed to run the following command:
```
vcf_dedupper --input-vcf $file --output-vcf ${file}.dedupped.vcf --variant-caller strelka --selection-method af --equality-mode 1 --sample-idx 1  --sort 2> ${file}.duplications.vcf  | bgzip > ${file}.dedupped.vcf.gz
```

For Starling it is agreed to run the following commands:
```
vcf_dedupper --input-vcf $file --output-vcf ${file}.dedupped.vcf --variant-caller starling --selection-method quality --equality-mode 1 --sort 2> ${file}.duplications.vcf  | bgzip > ${file}.dedupped.vcf.gz
```

For Platypus it is agreed to run the following command:
```
vcf_dedupper --input-vcf $file --variant-caller platypus --selection-method allele_calls --equality-mode 1 --sort 2> ${file}.duplications.vcf  | bgzip > ${file}.dedupped.vcf.gz
```

For the generic variant caller it is agreed to run the following command:
```
vcf_dedupper --input-vcf $file --variant-caller generic --selection-method af|quality  --equality-mode 1 --sort 2> ${file}.duplications.vcf  | bgzip > ${file}.dedupped.vcf.gz
```

For the `duplication_finder` run the following command:
```
vcf_dedupper --input-vcf $file --output-vcf ${file}.duplications.vcf --variant-caller duplication_finder --selection-method af --sort --sort-threads 2 2> ${file}.log
```

## Known issues

* When the input VCF does not have header (i.e.: no line starting with '#') the first character in the first line of the output VCF is replaced with the character '#'.
* When running unit tests simultaneously the resulting VCFs contain variants from different tests. This can be avoided by running unit tests separately. This behaviour only affects the test environment.
* The time to compute a given VCF is determined by the number of variants and the number of samples. Even if for some execution modes (e.g.: generic mode) the sample information is not used this affects the time to compute. We have observed that processing VCFs with a huge amount of samples (i.e.: 6773) have very low performance (0.2 variants/second).
