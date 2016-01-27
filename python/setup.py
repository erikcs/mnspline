#build locally with $python setup.py build_ext --inplace
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os
import numpy

os.environ["CC"] = "gcc-5" # gcc from brew

ext_modules = [Extension("pymnspline",
                        ["pymnspline.pyx", "../src/mnspline.c"],
                        language='c',
                        extra_compile_args=["-fopenmp"],
                        extra_link_args=["-fopenmp"],
                        include_dirs=[numpy.get_include()]
                        )]

setup(
     name = 'pymnspline',
     cmdclass = {'build_ext': build_ext},
     ext_modules = ext_modules,
     url='https://github.com/nuffe/mnspline'
)

