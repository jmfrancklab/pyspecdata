#from setuptools import setup
import setuptools # I think this is needed for the following
from numpy.distutils.core import Extension,setup
from distutils.spawn import find_executable
import subprocess
import sys
import os

general_error = "I couldn't import {:s} -- go install it first!!\n(I'm doing this because dependency-based install of PyQt, and some others does not usually go well -- use your distro software (conda install ..., aptitude, etc) instead)\nIn fact, you probably want to install:\n\tpyqt, unxutils, matplotlib, mingw, and libpython\nAlso, it's recommended to start by running ``python setup.py config --fcompiler=gfortran develop``"
try:
    import matplotlib
    #import PyQt5
except:
    raise RuntimeError(general_error.format('matplotlib'))
ext_modules = []
exec(compile(open('pyspecdata/version.py', "rb").read(), 'pyspecdata/version.py', 'exec'))

ext_modules.append(Extension(name = 'pyspecdata._nnls',
        sources = ['nnls/nnls.pyf','nnls/nnls.f','nnls/nnls_regularized.f90','nnls/nnls_regularized_loop.f90'],
        define_macros = [('ADD_UNDERSCORE',None)],
        extra_compile_args = ['-g'],# debug flags
        ))

setup(
    name='pySpecData',
    author='J. M. Franck',
    version=__version__,
    packages=setuptools.find_packages(),
    license='LICENSE.md',
    author_email='jmfranck@notgiven.com',
    url='http://github.com/jmfranck/pyspecdata',
    description='object-oriented N-dimensional data processing with notebook functionality',
    long_description=open('README.rst',encoding='utf-8').read(),
    install_requires=[
        "sympy",
        "numpy",
        "scipy",
        "h5py",
        "matplotlib",
        "pillow",
        ],
    ext_modules = ext_modules,
    entry_points=dict(console_scripts=
        ['scons_continuous=pyspecdata.latexscripts:repeat_scons',
        'update_notebook_pythonscripts=pyspecdata.latexscripts:main',
        'pdflatex_notebook_wrapper=pyspecdata.latexscripts:wraplatex',
        'pdflatex_notebook_view_wrapper=pyspecdata.latexscripts:wrapviewer',
        'pyspecdata_dataconfig=pyspecdata.genconfig:genconfig']
        ),
)
tryagain = False

print("You can now run pyspecdata_dataconfig to generate a template configuration file (which will show up in your home directory).")
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
