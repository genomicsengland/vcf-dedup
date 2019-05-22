from distutils.core import setup
from setuptools import find_packages

setup(
    name='vcf-dedup',
    version='0.4.8',
    description='Removal of duplicated variants from a VCF',
    packages=find_packages(),
    scripts=['scripts/vcf_dedupper'],
    url='https://github.com/genomicsengland/vcf-dedup',
    download_url='https://github.com/genomicsengland/vcf-dedup/archive/v0.4.8.tar.gz',
    license='Apache',
    author='Pablo Riesgo Ferreiro',
    author_email='pablo.riesgo-ferreiro@genomicsengland.co.uk',
    requires=['PyVCF', 'argparse'],
    install_requires=['PyVCF==0.6.8'],
    keywords=['VCF']
)
