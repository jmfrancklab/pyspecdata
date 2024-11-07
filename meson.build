project('pyspecdata', ['fortran', 'c'], version: '0.9.5.5.2', meson_version: '>=0.63.0')

py = import('python').find_installation('python3')

# Retrieve NumPy include directory
numpy_include_dir = run_command(
  py, '-c', 'import numpy; print(numpy.get_include())'
).stdout().strip()
# Find the path to numpy/f2py/src where fortranobject.h resides
f2py_include = run_command(py, '-c', 'import numpy.f2py; import os; print(os.path.join(os.path.dirname(numpy.f2py.__file__), "src"))', check: true).stdout().strip()

python_inc_dir = run_command(py, ['-c', 'from sysconfig import get_paths; print(get_paths()["include"])']).stdout().strip()

# Add linker arguments for linking against the Fortran runtime
add_project_link_arguments('-lgfortran', language: 'fortran')

# Custom target to copy fortranobject.h into the build directory
fortran_header = custom_target(
  'copy_fortranobject',
  input: f2py_include + '/fortranobject.h',
  output: 'fortranobject.h',
  command: ['cp', '@INPUT@', meson.project_build_root()],
)

# Custom target to compile fortranobject.c into an object file
fortranobject_o = custom_target(
  'fortranobject.o',
  input: join_paths(f2py_include, 'fortranobject.c'),
  output: 'fortranobject.o',
  command: ['gcc', '-fPIC', '-c', '@INPUT@', '-o', '@OUTPUT@', '-I' + f2py_include, '-I'+python_inc_dir],
  build_by_default: true
)

# Define custom target to generate _nnlsmodule.c from .pyf file
f2py_wrapper = custom_target(
  '_nnlsmodule.c',
  input: 'nnls/nnls.pyf',
  output: '_nnlsmodule.c',
  command: [py, '-m', 'numpy.f2py', '@INPUT@'],
  capture: false,
)

# Define the Fortran extension module
nnls_extension = py.extension_module(
  '_nnls',
  sources: [
    f2py_wrapper,
    fortran_header,
    fortranobject_o,
    'nnls/nnls.f',
    'nnls/nnls_regularized.f90',
    'nnls/nnls_regularized_loop.f90',
  ],
  include_directories: include_directories(numpy_include_dir),
  install: true,
  cpp_args: ['-DADD_UNDERSCORE', '-I'+numpy_include_dir],  # Equivalent to define_macros=[("ADD_UNDERSCORE", None)]
  fortran_args: ['-fallow-argument-mismatch'],  # Equivalent to extra_f77_compile_args
)

subdir('pyspecdata')
