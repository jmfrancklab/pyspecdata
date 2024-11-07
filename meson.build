project('pyspecdata', ['fortran', 'c'], version: '0.9.5.5.1', meson_version: '>=0.63.0')

py = import('python').find_installation('python3')

# Retrieve NumPy include directory
numpy_include_dir = run_command(
  py, '-c', 'import numpy; print(numpy.get_include())'
).stdout().strip()

# Find the path to numpy/f2py/src where fortranobject.h resides
f2py_include = run_command(py, '-c', 'import numpy.f2py; import os; print(os.path.join(os.path.dirname(numpy.f2py.__file__), "src"))', check: true).stdout().strip()

# Custom target to copy fortranobject.h into the build directory
fortran_header = custom_target(
  'copy_fortranobject',
  input: f2py_include + '/fortranobject.h',
  output: 'fortranobject.h',
  command: ['cp', '@INPUT@', meson.project_build_root()],
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
  'pyspecdata._nnls',
  sources: [
    f2py_wrapper,
    fortran_header,
    'nnls/nnls.f',
    'nnls/nnls_regularized.f90',
    'nnls/nnls_regularized_loop.f90',
  ],
  include_directories: include_directories(numpy_include_dir),
  install: true,
  cpp_args: ['-DADD_UNDERSCORE', '-I'+numpy_include_dir],  # Equivalent to define_macros=[("ADD_UNDERSCORE", None)]
  fortran_args: ['-fallow-argument-mismatch'],  # Equivalent to extra_f77_compile_args
)

py_sources = [
'pyspecdata/fl_dummy.py',
'pyspecdata/file_saving/hdf_save_dict_to_group.py',
'pyspecdata/axis_class.py',
'pyspecdata/datadir.py',
'pyspecdata/axis_manipulation/shear.py',
'pyspecdata/axis_manipulation/__init__.py',
'pyspecdata/axis_manipulation/secsy.py',
'pyspecdata/axis_manipulation/inhomog_coords.py',
'pyspecdata/axis_manipulation/register_axis.py',
'pyspecdata/fornotebook.py',
'pyspecdata/generate_fake_data.py',
'pyspecdata/hdf_utils.py',
'pyspecdata/general_functions.py',
'pyspecdata/core.py',
'pyspecdata/__init__.py',
'pyspecdata/placeholder.py',
'pyspecdata/plot_funcs/__init__.py',
'pyspecdata/plot_funcs/image.py',
'pyspecdata/plot_funcs/pcolormesh.py',
'pyspecdata/lmfitdata.py',
'pyspecdata/dict_utils.py',
'pyspecdata/version.py',
'pyspecdata/matrix_math/svd.py',
'pyspecdata/matrix_math/__init__.py',
'pyspecdata/matrix_math/dot.py',
'pyspecdata/matrix_math/nnls.py',
'pyspecdata/units.py',
'pyspecdata/mpl_utils.py',
'pyspecdata/DCCT_function.py',
'pyspecdata/figlist.py',
'pyspecdata/ndshape.py',
'pyspecdata/fourier/ft.py',
'pyspecdata/fourier/convolve.py',
'pyspecdata/fourier/shear.py',
'pyspecdata/fourier/ft_shift.py',
'pyspecdata/fourier/__init__.py',
'pyspecdata/fourier/ift.py',
'pyspecdata/fourier/_ft_conj.py',
'pyspecdata/fourier/ftshift.py',
'pyspecdata/load_files/bruker_nmr.py',
'pyspecdata/load_files/prospa.py',
'pyspecdata/load_files/open_subpath.py',
'pyspecdata/load_files/load_cary.py',
'pyspecdata/load_files/bruker_esr.py',
'pyspecdata/load_files/__init__.py',
'pyspecdata/load_files/acert.py',
'pyspecdata/nnls.py',
'pyspecdata/ipy.py',
'pyspecdata/latexscripts.py',
]
py.install_sources(py_sources, subdir: 'pyspecdata')
