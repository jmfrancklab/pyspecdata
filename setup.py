#from setuptools import setup
import setuptools # I think this is needed for the following
from numpy.distutils.core import Extension,setup
from distutils.spawn import find_executable
import subprocess
import sys
import os
setup(
    name='paramset_pyspecdata',
    author='J. M. Franck',
    version='0.1.0',
    packages=['paramset_pyspecdata'],
    license='LICENSE.md',
    description='helper for pyspecdata',
    long_description="just a module for storing global variables -- needed to change what's imported for notebook vs. graphical display",
)

try:
    import PyQt4.QtCore
    import PyQt4.QtGui
except:
    raise RuntimeError("I couldn't import PyQt -- go install it first!!\n(I'm doing this because dependency-based install of PyQt does not usually go well -- use your distro software (conda install ..., aptitude, etc) instead)\nIn fact, you probably want to install:\n\tpyqt, mayavi, unxutils, and libpython")
try:
    import mayavi
except:
    raise RuntimeError("I couldn't import MayaVi -- go install it first!!\n(I'm doing this because dependency-based install of MayaVi does not usually go well -- use your distro (conda install ..., aptitude, etc) instead)")
ext_test = Extension(name = 'pyspecdata.test_module',
        sources = ['pyspecdata/test_f90.pyf','pyspecdata/test_f90.f90','pyspecdata/anothertest.f90','pyspecdata/lprmpt.c','pyspecdata/fortrancall.h'],
        define_macros = [('ADD_UNDERSCORE',None)],
        )
ext_modules = [ext_test]

libraries = ['lapack']
lapack_success = False
skip_prop = False

if os.name == 'nt':
    if find_executable('make') is None:
        print "It looks like you're on windows, but I can't find make, so I'm skipping lapack"
        print "- If you have python(x,y), make should be installed, so this is a problem with your path"
        print "- If you have anaconda, you can install make with `conda install unxutils`"
        print "\nPress enter to acknowledge or Cntrl-C to break, fix this and try again"
        sys.stdin.readline(1)
        skip_prop = True
    else:
        libraries.append('refblas')
        print "It looks like you're on windows, so I'm going to build lapack (from http://netlib.org/) on MinGW."
        gfortran_name = find_executable('gfortran')
        if gfortran_name is None:
            raise ValueError("Install gfortran!\nFor Anaconda on windows, that means you probably need to run 'conda install libpython'")
        target = os.path.dirname(os.path.dirname(gfortran_name)) + os.sep + 'lib' + os.sep + 'gcc' # switched to gfortran because if you have two gcc's (eg. haskell from pandoc compile), this breaks
        if not os.path.exists(target):
            target = os.path.dirname(os.path.dirname(find_executable('gfortran'))) + os.sep + 'libs' # works for anaconda -- switched to gfortran because if you have two gcc's (eg. haskell from pandoc compile), this breaks
        if 'liblapack.a' in os.listdir(target) and 'librefblas.a' in os.listdir(target):
            print "I see liblapack and librefblas in",target,"so I'm not going to rebuild"
            lapack_success = True
        else:
            os.chdir('lapack-3.4.0')
            my_env = os.environ.copy()
            ar_directory = os.path.dirname(find_executable('ar'))
            my_env['PATH'] += ';'+ar_directory
            print "new path is",my_env['PATH']
            subprocess.Popen("cmd /c makelibs.bat",env = my_env).wait()
            # possibly replace with start
            # possibly look at this: https://gcc.gnu.org/wiki/LAPACK%20on%20Windows
            print "trying to move the built libraries to "+target
            os.rename('liblapack.a',target + os.sep +'liblapack.a')
            os.rename('librefblas.a',target + os.sep +'librefblas.a')
            print "press enter..."
            sys.stdin.readline(1)
            os.chdir('..')
if not skip_prop:
    ext_prop = Extension(name = 'pyspecdata.propagator',
            sources = ['pyspecdata/propagator.f90'],
            libraries = libraries)
    ext_modules.append(ext_prop)

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
    except (RuntimeError, TypeError, NameError),e:
        if not lapack_success and len(ext_modules) == 2:
            print "something went wrong, so I'm going to try to rebuild without lapack -- this means you won't be able to use the prop module"
            print "this is OK, but press enter to acknowledge"
            sys.stdin.readline(1)
            ext_modules = [ext_test]
            tryagain = True
        elif not lapack_success and len(ext_modules) == 1:
            print "Error message was:"
            print e
            print "something STILL went wrong, so I'm rebuilding without the fortran test module"
            print "this is OK, but press enter to acknowledge"
            sys.stdin.readline(1)
            ext_modules = []
            tryagain = True
        else:
            print "I failed with no extension modules or with lapack_success"
            raise

# Later, I should probably use the setuptools equivalent of install_data to do both this and the lapack stuff
print "\n\nNow that everything else is set up, I'm going to check your notebook and data directories, possibly asking you to set them."
import pyspecdata.datadir as d
print "Trying to grab the data directory:"
d.getDATADIR()
print "Trying to grab the notebook directory:"
d.get_notebook_dir()
