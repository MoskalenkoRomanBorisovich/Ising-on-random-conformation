from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import numpy

import mc_lib

ising = Extension("cy_ising", ["cy_ising.pyx"],
                include_dirs = [numpy.get_include(),
                                mc_lib.get_include()],
                language='c++',)


ising_cluster = Extension("cy_ising_cluster", ["cy_ising_cluster.pyx"],
                include_dirs = [numpy.get_include(),
                                mc_lib.get_include()],
                language='c++',)



ising_exact = Extension("exact_ising", ["exact_ising.pyx"],
                include_dirs = [numpy.get_include()],
                language='c++',)


setup(ext_modules=cythonize([ising, ising_cluster, ising_exact], build_dir='ising_builds'),
      cmdclass = {'build_ext': build_ext})
