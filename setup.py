#!/usr/bin/env python
"""
setup.py adapted from https://github.com/pypa/sampleproject

"""

from __future__ import absolute_import, print_function, division
from setuptools import setup, find_packages
from codecs import open
from os import path

# Fetch long_description from DESCRIPTION.rst.
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='binstarsolver',
    # Versions should comply with PEP440.
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.0.1',
    description='Estimate physical quantities of a binary star system from observed quantities.',
    long_description=long_description,
    url='https://github.com/ccd-utexas/binstarsolver',
    author='White Dwarf Research Group, Astronomy Dept, UT Austin',
    author_email='harrold@astro.as.utexas.edu',
    license='MIT',
    classifiers=[
        # Project maturity values: 3 - Alpha; 4 - Beta; 5 - Production/Stable
        # https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Scientists',
        'Topic :: Astronomy',
        'License :: OSI Approved :: MIT License',
        # Supported Python versions.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
    ],
    keywords='astronomy binary stars observing',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    # Run-time dependencies. Will be installed by pip.
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[
        'numpy>=1.8.2',
        'scipy>=0.14.0',
        'astropy>=0.4.2',
        'matplotlib>=1.4.0'],
    # Additional dependencies (e.g. development dependencies).
    # Install using: $ pip install -e .[dev,test]
    extras_require = {
        'dev': ['check-manifest'],
        'test': ['coverage', 'pytest'],
    },
    # Data files included in installed packages.
    package_data={
        # 'sample': ['package_data.dat'],
    },
    # Data files not included in installed packages.
    # 'package_data' is preferred approach.
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
    # Example: 'data_file' will be installed into '<sys.prefix>/my_data'
    data_files=[
        #('my_data', ['data/data_file'])
        ],
    # For executable scripts, use 'entry_points' rather than 'scripts'
    # for cross-platform portability.
    entry_points={
        #'console_scripts': [
        #    'sample=sample:main',
        #],
    },
)
