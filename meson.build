
project('pyspecdata', 'fortran', version: '0.9.5.4', meson_version: '>=0.63.0')

py = import('python').find_installation('python3')

# Retrieve NumPy include directory
numpy_include_dir = run_command(
  py, '-c', 'import numpy; print(numpy.get_include())'
).stdout().strip()

# Define custom target to generate _nnlsmodule.c from .pyf file
f2py_wrapper = custom_target(
  '_nnlsmodule.c',
  input: 'nnls/nnls.pyf',
  output: '_nnlsmodule.c',
  command: [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'pyspecdata._nnls', '--quiet'],
  capture: true,
)

# Define the Fortran extension module
nnls_extension = py.extension_module(
  'pyspecdata._nnls',
  sources: [
    f2py_wrapper,
    'nnls/nnls.f',
    'nnls/nnls_regularized.f90',
    'nnls/nnls_regularized_loop.f90',
  ],
  include_directories: include_directories(numpy_include_dir),
  install: true,
)

# Install the Python sources
py.install_sources('pyspecdata')
