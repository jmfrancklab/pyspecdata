import setuptools  # I think this is needed for the following
from numpy.distutils.core import Extension, setup
from pyspecdata.version import __version__
from distutils.spawn import find_executable
import os

ext_modules = []
exec(
    compile(
        open("pyspecdata/version.py", "rb").read(),
        "pyspecdata/version.py",
        "exec",
    )
)

ext_modules.append(
    Extension(
        name="pyspecdata._nnls",
        sources=[
            "nnls/nnls.pyf",
            "nnls/nnls.f",
            "nnls/nnls_regularized.f90",
            "nnls/nnls_regularized_loop.f90",
        ],
        define_macros=[("ADD_UNDERSCORE", None)],
        # extra_compile_args = ['-g'],# debug flags
        extra_f77_compile_args=[
            "-fallow-argument-mismatch"
        ],  # seems to be required on linux, but doesn't work on windows
    )
)
setup(
    name="pySpecData",
    author="J. M. Franck",
    version=__version__,
    packages=setuptools.find_packages(),
    license="LICENSE.md",
    author_email="jmfranck@notgiven.com",
    url="http://github.com/jmfranck/pyspecdata",
    description=(
        "object-oriented N-dimensional data processing with notebook"
        " functionality"
    ),
    long_description=open("README.rst", encoding="utf-8").read(),
    install_requires=[
        "sympy",
        "numpy",
        "scipy",
        "h5py",
        "matplotlib>=3.8.0",  # sharex is a problem in 3.5-6
        "pillow",
        "pint",
        "lmfit>=1.1",  # we recently found that at least 1.0.3 generates
        #                output parameters that are unchanged, but still
        #                returns a "success" condition
    ],
    ext_modules=ext_modules,
    entry_points=dict(
        console_scripts=[
            "scons_continuous=pyspecdata.latexscripts:repeat_scons",
            "update_notebook_pythonscripts=pyspecdata.latexscripts:main",
            "pdflatex_notebook_wrapper=pyspecdata.latexscripts:wraplatex",
            ("pdflatex_notebook_view_wrapper="
             "pyspecdata.latexscripts:wrapviewer"),
            "pyspecdata_dataconfig=pyspecdata.latexscripts:genconfig",
            ("pyspecdata_register_dir="
             "pyspecdata.datadir:register_directory"),
        ]
    ),
)

print(
    "You can now run pyspecdata_dataconfig to generate a template"
    " configuration file (which will show up in your home directory)."
)
