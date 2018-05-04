import numpy
from distutils.core import setup, Extension

setup(name='PEhmm2.c',version=1.0,ext_modules=[Extension('PEhmm2', ['PEhmm2.c'],include_dirs=[numpy.get_include()])])
