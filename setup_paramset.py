import setuptools # I think this is needed for the following
from numpy.distutils.core import Extension,setup
setup(
    name='paramset_pyspecdata',
    author='J. M. Franck',
    version='0.1.0',
    packages=['paramset_pyspecdata'],
    license='LICENSE.md',
    description='helper for pyspecdata',
    long_description="just a module for storing global variables -- needed to change what's imported for notebook vs. graphical display",
)
