#!/usr/bin/env python
"""Setup script adapted from [1]_, [2]_.

References
----------
.. [1] https://github.com/pypa/sampleproject
.. [2] https://bitbucket.org/birkenfeld/sphinx-contrib/src/dae39f177383/napoleon/?at=default

"""

from __future__ import absolute_import, division, print_function
import os
import codecs
import setuptools

# Load descriptions.
# Note: As of 2015-01-25, check-manifest will not load docs/requirements_readthedocs.txt.
# To build docs on local machine, install dependencies manually.
fpath = os.path.abspath(os.path.dirname(__file__))
long_description = codecs.open(os.path.join(fpath, 'DESCRIPTION.rst'), encoding='utf-8').read()

setuptools.setup(
    name='binstarsolver',
    # Versions should comply with PEP440.
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.1.3',
    description='Estimate physical quantities of a binary star system from observed quantities.',
    long_description=long_description,
    url='http://binstarsolver.readthedocs.org',
    author='Samuel Harrold; White Dwarf Research Group, Astronomy Dept, UT Austin',
    author_email='samuel.harrold@gmail.com',
    license='MIT',
    classifiers=[
        # Project maturity values: 3 - Alpha; 4 - Beta; 5 - Production/Stable
        # https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: MIT License',
        # Supported Python versions.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
    ],
    keywords='astronomy astrophysics binary stars observing',
    packages=setuptools.find_packages(exclude=['contrib', 'docs', 'tests*']),
    # Run-time dependencies. Will be installed by pip.
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['numpy>=1.8.2', 'scipy>=0.14.0', 'matplotlib>=1.4.0'],
    # Additional dependencies (e.g. development dependencies).
    # Install using: $ pip install -e .[dev,test]
    extras_require = {
        'dev': ['check-manifest>=0.22'],
        'test': ['coverage>=3.7.1', 'pytest>=2.6.3']
    },
    # Data files included in installed packages.
    package_data={},
    # Data files not included in installed packages. 'package_data' is preferred approach.
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
    # Example: 'data_file' will be installed into '<sys.prefix>/my_data'
    data_files=[],
    # For executable scripts, use 'entry_points' rather than 'scripts'
    # for cross-platform portability.
    entry_points={
        'console_scripts': [
            'binstarsolver=binstarsolver.main:main'
        ],
    },
)
