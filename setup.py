#from setuptools import setup
import setuptools # I think this is needed for the following
from numpy.distutils.core import Extension,setup
import subprocess
import sys
import os
ext_test = Extension(name = 'pyspecdata.test_module',
        sources = ['pyspecdata/test_f90.pyf','pyspecdata/test_f90.f90','pyspecdata/anothertest.f90','pyspecdata/lprmpt.c','pyspecdata/fortrancall.h'],
        define_macros = [('ADD_UNDERSCORE',None)],
        )
ext_modules = [ext_test]

kwargs = {}
if os.name == 'nt':
    print "It looks like you're on windows, so I'm going to build lapack (from http://netlib.org/) on MinGW."
    kwargs.update(dict(library_dirs = ['lapack-3.4.0']))
    os.chdir('lapack-3.4.0')
    subprocess.call(["cmd","makelibs.bat"])
    os.chdir('..')
try:
    ext_prop = Extension(name = 'pyspecdata.propagator',
            sources = ['pyspecdata/propagator.f90'],
            libraries = ['lapack','refblas'],
            **kwargs)
    ext_modules.append(ext_prop)
except:
    print "Lapack didn't build, so I'm not going to try to load the prop module."
    print "this is OK, but press enter to acknowledge"
    sys.stdin.readline(1)

tryagain = True
while tryagain == True:
    try:
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
            ext_modules = ext_modules,
            tryagain = False
        #    entry_points=dict(
        #        notebook_info=["data_dir = pyspecdata:datadir ["+os.path.expanduser('~')+os.path.sep+'exp_data]']
        #        )
        )
    except:
        if len(ext_modules) == 2:
            print "something went wrong, so I'm going to try to rebuild without lapack -- this means you won't be able to use the prop module"
            print "this is OK, but press enter to acknowledge"
            sys.stdin.readline(1)
            ext_modules = [ext_test]
            tryagain = True
        elif len(ext_modules) == 1:
            print "something STILL went wrong, so I'm rebuilding without the fortran test module"
            print "this is OK, but press enter to acknowledge"
            sys.stdin.readline(1)
            ext_modules = []
            tryagain = True
