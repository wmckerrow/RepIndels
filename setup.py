import numpy
from distutils.core import setup, Extension

setup(name='PEhmm.c',version=1.0,ext_modules=[Extension('PEhmm', ['PEhmm.c'],include_dirs=[numpy.get_include()])])
