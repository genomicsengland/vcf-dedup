from distutils.core import setup
from setuptools import find_packages

setup(
    name='vcf-dedup',
    version='0.3.1',
    packages=find_packages(),
    scripts=['scripts/vcf_dedupper'],
    url='',
    license='',
    author='priesgo',
    author_email='pablo.ferreiro@genomicsengland.co.uk',
    description='', requires=['PyVCF'],
    install_requires=[
        'PyVCF==0.6.8'
    ]
)
