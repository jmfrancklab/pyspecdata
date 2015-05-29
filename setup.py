from distutils.core import setup

setup(
    name='pySpecData',
    author='J. M. Franck',
    version='0.2.0',
    packages=['pyspecdata','pyspecdata.propagate'],
    license='LICENSE.txt',
    description='object-oriented N-dimensional data processing with notebook functionality',
    long_description=open('README.txt').read(),
    install_requires=[
        "sympy",
        "numpy",
        "scipy",
        "matplotlib",
        "pytables",
        ],
)
