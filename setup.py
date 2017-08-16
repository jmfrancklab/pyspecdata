#from setuptools import setup
import setuptools # I think this is needed for the following
from numpy.distutils.core import Extension,setup
from distutils.spawn import find_executable
import subprocess
import sys
import os

general_error = "I couldn't import {:s} -- go install it first!!\n(I'm doing this because dependency-based install of PyQt, and some others does not usually go well -- use your distro software (conda install ..., aptitude, etc) instead)\nIn fact, you probably want to install:\n\tpyqt, unxutils, matplotlib, mingw, and libpython"
try:
    import paramset_pyspecdata
except:
    raise RuntimeError("install the paramset_pyspecdata package first!\nIf using setup.py, run 'python setup_paramset.py install'")
try:
    import matplotlib
except:
    raise RuntimeError(general_error.format('matplotlib'))
ext_modules = []

setup(
    name='pySpecData',
    author='J. M. Franck',
    version='0.9.5.1.3',
    packages=setuptools.find_packages(exclude=['paramset_pyspecdata',]),
    license='LICENSE.md',
    author_email='jmfranck@notgiven.com',
    url='http://github.com/jmfranck/pyspecdata',
    description='object-oriented N-dimensional data processing with notebook functionality',
    long_description=open('README.rst').read(),
    install_requires=[
        "paramset_pyspecdata",
        "sympy",
        "numpy",
        "scipy",
        "h5py",
        "matplotlib",
        "tables",
        ],
    ext_modules = ext_modules,
    entry_points=dict(console_scripts=
        ['scons_continuous=pyspecdata.latexscripts:repeat_scons',
        'update_notebook_pythonscripts=pyspecdata.latexscripts:main',
        'pdflatex_notebook_wrapper=pyspecdata.latexscripts:wraplatex',
        'pdflatex_notebook_view_wrapper=pyspecdata.latexscripts:wrapviewer']
        ),
)
tryagain = False

## Later, I should probably use the setuptools equivalent of install_data to do both this and the lapack stuff
#print "\n\nNow that everything else is set up, I'm going to check your notebook and data directories, possibly asking you to set them."
#print "\n--> The notebook directory is the root directory where you store tex files for your notebook."
#print "\n--> The data directory is the root directory where you store your data."
#print "\nFor the notebook and data directory, I'm going to assume that there is one root directory with, possibly many, subfolders underneath.  But, you must choose directories that exist.  Note that you can use ~ for your home directory, even on windows."
#import pyspecdata.datadir as d
#print "Trying to grab the data directory:"
#d.getDATADIR()
#print "Trying to grab the notebook directory:"
#d.get_notebook_dir()
