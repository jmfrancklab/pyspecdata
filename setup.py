from setuptools import setup

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
)
