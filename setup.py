#from setuptools import setup
import setuptools # I think this is needed for the following
from numpy.distutils.core import Extension,setup
import os
#ext_prop = Extension(name = 'propagator',
#        sources = ['propagator.f90'],
#        f2py_options = '--fcompiler=gnu95 -llapack',
#        )
ext_test = Extension(name = 'pyspecdata.test_module',
        sources = ['pyspecdata/test_f90.pyf','pyspecdata/test_f90.f90','pyspecdata/anothertest.f90','pyspecdata/lprmpt.c','pyspecdata/fortrancall.h'],
        define_macros = [('ADD_UNDERSCORE',None)],
        )
#f2py_options = ['--fcompiler=gnu95','--compiler=mingw32','-lmsvcr71'])

kwargs = {}
if os.name == 'nt':
    print "It looks like you're on windows, so I'm going to try to use a version of lapack (from http://netlib.org/) that was built with MinGW.  If this fails, you should get MinGW"
    kwargs.update(dict(library_dirs = ['.']))

try:
    ext_prop = Extension(name = 'pyspecdata.propagator',
            sources = ['pyspecdata/propagator.f90'],
            libraries = ['lapack','refblas'],
            **kwargs)
except:
    raise RuntimeError("There was as an error -- probabably becuase of lapack\nIf you want to no more, remove the try/except clause in setup.py")

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
        "mayavi",
        "PyQt4",
        ],
    ext_modules = [ext_test,ext_prop],
#    entry_points=dict(
#        notebook_info=["data_dir = pyspecdata:datadir ["+os.path.expanduser('~')+os.path.sep+'exp_data]']
#        )
)
