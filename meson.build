project('pyspecdata', ['fortran', 'c'], version: '0.9.5.5.2', meson_version: '>=0.63.0')

py = import('python').find_installation('python3')

# Retrieve NumPy include directory
numpy_include_dir = run_command(
  py, '-c', 'import numpy; print(numpy.get_include())'
).stdout().strip()
numpy_inc = run_command(py, '-c', 'import numpy,sys;sys.stdout.write(numpy.get_include())').stdout().strip() 
# Find the path to numpy/f2py/src where fortranobject.h resides
fortranobject_src = run_command(py, '-c', 'import numpy,os,sys;import os;sys.stdout.write(os.path.join(os.path.dirname(numpy.__file__),"f2py","src","fortranobject.c"))').stdout().strip()
numpy_f2py_src = run_command(py, '-c', 'import numpy,os,sys;sys.stdout.write(os.path.join(os.path.dirname(numpy.__file__),"f2py","src"))').stdout().strip()
cc = meson.get_compiler('c')
fc = meson.get_compiler('fortran')
python_dep = py.dependency()
fortran_sources = files(
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
  )
python_inc_dir = run_command(py, ['-c', 'from sysconfig import get_paths; print(get_paths()["include"])']).stdout().strip()

c_wrappers = custom_target(
    'f2py_wrappers',
    # interface blocks require the modules containing fit variables
    input: [
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
    output: ['_nnlsmodule.c', '_nnls-f2pywrappers2.f90'],
    command: [py, '-m', 'numpy.f2py', '@INPUT@', '-m', '_nnls', '--lower'],
    build_by_default: true,
)

fixed_lib = static_library(
    'nnls_fixed',
    fortran_sources,
    fortran_args: ['-O2','-g','-fopenmp','-std=gnu','-ffixed-form','-ffixed-line-length-none'],
    pic: true,
    install: false,
)

mod_dir = join_paths(meson.current_build_dir(), 'libnnls_fixed.a.p')

c_sources = files(
    fortranobject_src,
)

c_lib = static_library(
    'nnls_c',
    c_sources + [c_wrappers[0]],
    c_args: ['-DADD_UNDERSCORE', '-I' + numpy_inc, '-I' + numpy_f2py_src],
    dependencies: [python_dep],
    pic: true,
    install: false,
)

wrapper_lib = static_library(
    'nnls_wrapper',
    c_wrappers[1],
    fortran_args: ['-O2','-g','-std=gnu', '-I' + mod_dir],
    dependencies: [python_dep],
    pic: true,
    install: false,
)

py.extension_module(
    '_nnls',
    [],
    objects: [
        fixed_lib.extract_all_objects(),
        wrapper_lib.extract_all_objects(),
        c_lib.extract_all_objects(),
    ],
    c_args: ['-I'+numpy_inc, '-I'+numpy_f2py_src],
    subdir: 'nnls',
    install: true,
)

fortran_args = []
if host_machine.system() == 'windows'
    # Windows-specific flags can go here if needed
    # fortran_args += ['some-windows-specific-flag']
else
    # Assume Linux or other Unix-like systems
    fortran_args += ['-fallow-argument-mismatch'] # Equivalent to extra_f77_compile_args
endif

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
