#from setuptools import setup
import setuptools # I think this is needed for the following
from numpy.distutils.core import Extension,setup
#ext_prop = Extension(name = 'propagator',
#        sources = ['propagator.f90'],
#        f2py_options = '--fcompiler=gnu95 -llapack',
#        )
ext_test = Extension(name = 'pyspecdata.test_f90',
        sources = ['pyspecdata/test_f90.f90'])
#f2py_options = ['--fcompiler=gnu95','--compiler=mingw32','-lmsvcr71'])
ext_prop = Extension(name = 'pyspecdata.propagator',
        sources = ['pyspecdata/propagator.f90'],
        libraries = ['lapack','refblas'])

setup(
    name='pySpecData',
    author='J. M. Franck',
    version='0.1.0',
    packages=['pyspecdata'],
    license='LICENSE.md',
    description='object-oriented N-dimensional data processing with notebook functionality',
    long_description=open('README.rst').read(),
    install_requires=[
        "sympy",
        "numpy",
        "scipy",
        "matplotlib",
        "tables",
        ],
    ext_modules = [ext_test,ext_prop],
)
