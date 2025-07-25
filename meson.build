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
  command: [py, '-c', # python copy b/c cp doesn't work on windows
    'import shutil, sys; shutil.copy(sys.argv[1], sys.argv[2])',
    '@INPUT@', meson.project_build_root()],
)

# Custom target to compile fortranobject.c into an object file
fortranobject_o = custom_target(
  'fortranobject.o',
  input: join_paths(f2py_include, 'fortranobject.c'),
  output: 'fortranobject.o',
  command: ['gcc', '-fPIC', '-c', '@INPUT@', '-o', '@OUTPUT@', '-I' + f2py_include, '-I'+python_inc_dir, '-I'+numpy_include_dir],
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

fortran_args = []
if host_machine.system() == 'windows'
    # Windows-specific flags can go here if needed
    # fortran_args += ['some-windows-specific-flag']
else
    # Assume Linux or other Unix-like systems
    fortran_args += ['-fallow-argument-mismatch'] # Equivalent to extra_f77_compile_args
endif

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
    'nnls/venk_brd.f90',
    'lapack-3.4.0/SRC/dgesv.f',
    'lapack-3.4.0/SRC/dgetrf.f',
    'lapack-3.4.0/SRC/dgetrs.f',
    'lapack-3.4.0/SRC/dgetf2.f',
    'lapack-3.4.0/SRC/dlaswp.f',
    'lapack-3.4.0/SRC/ilaenv.f',
    'lapack-3.4.0/SRC/iparmq.f',
    'lapack-3.4.0/SRC/ieeeck.f',
    'lapack-3.4.0/INSTALL/dlamch.f',
    'lapack-3.4.0/BLAS/SRC/dgemm.f',
    'lapack-3.4.0/BLAS/SRC/dgemv.f',
    'lapack-3.4.0/BLAS/SRC/dtrsm.f',
    'lapack-3.4.0/BLAS/SRC/dger.f',
    'lapack-3.4.0/BLAS/SRC/dscal.f',
    'lapack-3.4.0/BLAS/SRC/dswap.f',
    'lapack-3.4.0/BLAS/SRC/idamax.f',
    'lapack-3.4.0/BLAS/SRC/daxpy.f',
    'lapack-3.4.0/BLAS/SRC/lsame.f',
    'lapack-3.4.0/BLAS/SRC/xerbla.f',
  ],
  include_directories: include_directories(numpy_include_dir),
  install: true,
  cpp_args: ['-DADD_UNDERSCORE', '-I'+numpy_include_dir],  # Equivalent to define_macros=[("ADD_UNDERSCORE", None)]
  fortran_args: fortran_args,
)

subdir('pyspecdata')

scripts = [
    ['update_notebook_pythonscripts', 'pyspecdata.latexscripts', 'main'],
    ['pdflatex_notebook_wrapper', 'pyspecdata.latexscripts', 'wraplatex'],
    ['pdflatex_notebook_view_wrapper', 'pyspecdata.latexscripts', 'wrapviewer'],
    ['pyspecdata_dataconfig', 'pyspecdata.latexscripts', 'genconfig'],
    ['pyspecdata_register_dir', 'pyspecdata.datadir', 'register_directory'],
]

foreach script : scripts
    configure_file(
        output: script[0],
        input: 'scripts/entry_point_template',
        configuration: {
            'python_executable': py.path(),
            'module': script[1],
            'function': script[2],
        },
        install: true,
        install_dir: get_option('bindir'),
    )
endforeach
