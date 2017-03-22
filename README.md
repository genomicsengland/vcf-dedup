# VCF dedupper

VCF files may contain duplicated after normalisation or merging of VCFs. This module intends to remove those duplicates in a sensible way. This process is very specific to the variant caller, we only support Strelka and Starling variant callers. Platypus might be supported in future releases. The variant caller might be selected with the parameter `--variant-caller`.

Variant duplicates are identified by having the same genomic coordinates and optionally the same alternate allele. This is controlled with the parameter `--equality-mode`.

The currrent implementation takes into account:
1. The filtering information. Whenever there is one and only one variant non filtered (i.e.: PASS) This will be the resulting variant.
2. Other criteria are used when there is a collision in filtering status (i.e.: two or more non filtered variants or all filtered variants). This is set with `--selection-method`
    * Highest allele frequency (i.e.: ratio of supporting reads for the variant call)
    * Highest variant calling quality.
    
The output VCF should not contain any duplicate.

## Assumptions

* Multi-allelic variants are not supported
* VCF must be sorted


## Using the VCF dedupper

The VCF dedupper comes in two flavors:
1. A command line script
2. A python module that can be used programmatically

### Script

```
vcf_dedupper --input-vcf $INPUT_VCF --output-vcf $OUTPUT_VCF --variant-caller $(one of "strelka", "starling" or "duplication_finder") --selection-method $(one of "af", "quality", "arbitrary")
```

The `duplication_finder` mode writes to the output only those variants that are duplicated without removing any duplicates. This mode is intended to facilitate the analysis of the root cause of duplications.

The `equality-mode` may be changed to use a less restrictive definition of duplicates using only chromosome, position and reference. By default it uses chromosome, position, reference and alternate.

The `sample-idx` indicates as a 0-based index for multisample VCFs the sample to which collision criteria will be applied. Beware that for Strelka VCFs the tumor sample is in the second sample (i.e.: index = 1), this might change.

The `sample-name` indicates the name of the sample to which collision criteria will be applied. This parameter overrides `sample-idx`. If the sample name does not exist an error is raised.

### Python module

```
from vcf_dedup.runner import VcfDedupRunner


config = {
        "input_vcf" : input_vcf,
        "output_vcf" : output_vcf,
        "variant_caller": variant_caller,
        "selection_method": selection_method,
        "equality_mode": equality_mode,
        "sample_idx": sample_idx
    }
# Calls the VCF dedupper
runner = VcfDedupRunner(config)
runner.process_vcf()
```

### Output

If the parameter `output-vcf` is provided this file will be created.
If the parameter `output-vcf` is not provided the resulting VCF will be written to the standard output.

### Sample commands

For strelka it is agreed to run the following command:
```
vcf_dedupper --input-vcf $file --output-vcf ${file}.dedupped.vcf --variant-caller strelka --selection-method af --equality-mode 1 --sample-idx 1 2> ${file}.log
```

For Starling it is agreed to run the following commands:
```
vcf_dedupper --input-vcf $file --output-vcf ${file}.dedupped.vcf --variant-caller starling --selection-method quality --equality-mode 1 2> ${file}.log
```

For the `duplication_finder` run the following command:
```
vcf_dedupper --input-vcf $file --output-vcf ${file}.duplications.vcf --variant-caller duplication_finder --selection-method af 2> ${file}.log
```

