#!/usr/bin/env python                                                                                                                                                                                           

from distutils.core import setup, Extension
import glob
import os

# Get matfiles and images for testing
matfiles=glob.glob(os.path.join('tests/matfiles/*.mat'))
images=glob.glob(os.path.join('tests/data/*'))
docs=glob.glob(os.path.join('documentation/*'))

setup(
    name='pyphasesym',
    version='1.0',
    description='Python implementation of phasesym program',
    author='Abhijit Bendale',
    author_email='bendale@mit.edu',
    py_modules = ['pyphasesym','tests.test_pyphasesym_regression',
                  'tests.test_pyphasesym_unit'],
    data_files = [('documentation',['documentation/pyphasesym.pdf'])]
    )
    
