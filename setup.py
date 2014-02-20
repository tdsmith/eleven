from __future__ import print_function
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import io
import codecs
import os
import sys

import eleven

here = os.path.abspath(os.path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding= kwargs.get('encoding', 'utf-8')
    sep= kwargs.get('sep', '\n')
    buf=[]
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.rst')

"""
class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_test(self):
        import unittest
        errcode = unittest.main()
        sys.exit(errcode)
"""

setup(
    name='eleven',
    version=eleven.__version__,
    url='http://github.com/tdsmith/eleven/',
    license='BSD',
    author='Tim D. Smith',
    author_email='eleven@tim-smith.us',
    description='A friendly implementation of the GeNorm multi-gene RT-qPCR normalization algorithm',
    long_description = long_description,
    packages = ['eleven'],
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        ],
    include_package_data=True,
    platforms='any',
    test_suite='eleven.test_eleven',
    install_requires=[
        'numpy>=1.7.0',
        'scipy>=0.12',
        'pandas>=0.12',
        ],
    # cmdclass={'test': PyTest},
)

