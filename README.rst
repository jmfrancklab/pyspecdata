.. image:: https://zenodo.org/badge/24356894.svg
   :target: https://zenodo.org/badge/latestdoi/24356894
   
To learn more about pyspecdata, you can head over to the `documentation <http://jmfrancklab.github.io/pyspecdata>`_.

If you already know that you want to install,
see `Installation <#installation>`_

Please note this package is heavily utilized by three other packages that our lab manages on github:

-   `ODNP processing scripts <https://github.com/jmfrancklab/proc_scripts/>`_
-   `Classes for communicating with instruments <https://github.com/jmfrancklab/inst_notebooks/>`_
-   `SpinCore Extension <https://github.com/jmfrancklab/spincore_apps/>`_

We have somewhat recently added fast compiled Fortran functions for things like
2D ILT (Tikhonov regularization with basis set compression) for NMR (Nuclear Magnetic Resonance),
so please read the install instructions
carefully!

===========
pySpecData
===========

Object-oriented Python package for processing spectral data -- or in general, *n*-dimensional data with labeled axes (i.e. *N*-dimensional gridded data or "nddata").
It depends on numpy, which provides very fast manipulations of *N*-dimensional gridded arrays ("ndarray").
This package has some overlap with xarray,
but it doesn't attempt to mimic pandas notation,
shooting instead for very compact notation for natural slicing, etc.
It mainly focuses on making it easier to *quickly write good code
for processing spectroscopy data*.
In particular, it takes care of various features related to fourier
transformation, error propagation, and direct products in multidimensional data with
little to no interaction from the user.

If you are working in a lab developing new spectroscopic methodologies, then this package is definitely for you.
If you deal with multi-dimensional data of some other form, then it's likely for you.
Features include:

Features
========

* Labeled axes allow one to manipulate datasets (potentially with different dimensions) without having to explicitly keep track of what the different dimensions correspond to.  Code becomes more legible.  Also, tiling, direct product, and gridding functions become obsolete.

* Fourier transformation with automatic manipulation of axes.

* Automatic error propagation.

* Commands like ``plot(data)`` will generate a plot with automatically labeled
  axes, errors, and units.
  All of this information is also written to HDF5 when the data is saved.

* Simplified curve fitting that takes advantage of labeled axes and Python's symbolic algebra package (sympy).

* The code is written so that it can be integrated into a nicely formatted PDF lab notebook.

    * The same code can be run on the command line (to generate pop-up plot windows) and embedded into a LaTeX document.

    * Extension to other output formats, such as HTML or markdown, should be relatively straightforward.

* In a multimedia environment like jupyter, you don't need a separate plot
  command.  The code can automatically choose a plotting style appropriate to
  the code (eventually, the general preferences for this can just be configured
  at the beginning of the jupyter notebook).

More detailed web documentation will be coming soon.

NMR/ESR specific
================

Because it was written primarily for NMR and ESR data, it also includes:

* Routines for reading commercial raw data (*e.g.* Bruker, Kea) into nddata
  objects with all relevant information.

* The object-oriented features make it much easier to process raw phase-cycled
  data and to simultaneously view multiple (potentially interfering) coherence
  pathways.

* Contains functions for baseline correction, peak integration, *etc.*

* (Not yet in packaged version) A basic compiled routine for propagating
  density matrices that can be used to predict the response to shaped pulses.

Version Notes
=============

Note that the current version is intended just for collaborators, *etc.*
(Though, if you do really want to use it for interesting science,
we are happy to work with you to make it work for your purposes.)
A public-use version 1.0.0, to be accompanied by useful demonstrations, is planned within a year.
*(Note that the email currently linked to the PyPI account is infrequently checked --if you have interest in this software, please find J. Franck's website and contact by that email.)*

History/Roadmap
---------------

(Current version in bold) 

0.9.5
    First version distributed on pypi.python.org.

0.9.5.1
    - 0.9.5.1.1
      Some important debugging, and also added `pyspecdata.ipy` → executing the following at the top of a jupyter notebook:

        .. code-block:: python

            %pylab inline
            %load_ext pyspecdata.ipy

      will cause nddata to "display" as labeled plots.

    - 0.9.5.1.2
      added ability to load power saturation 2D data from Bruker

    - 0.9.5.1.3
      XEpr data loaded with dBm units rather than W units

      added ``to_ppm`` function for Bruker files

    - 0.9.5.1.4
      Improved internal logging, and started to remove gratuitous dependencies,
      ``%load_ext pyspecdata.ipy`` includes
      ``%pylab inline``, so that only

        .. code-block:: python

            %load_ext pyspecdata.ipy

        is required for jupyter.
  
    - 0.9.5.1.6

        - Removed several legacy modules, and added docstrings for the remaining modules.

        - Begin phasing out earlier `CustomError` class.

        - Make `numpy` pretty printing available from the `general_functions` module.

        - Add xelatex support to the notebook wrapper.

        - Start to move file search routines away from demanding a single "data directory."

        - Improved support for 2D Bruker XEPR

        - Made it possible to call standard trig functions with `nddata` as an argument.
    - 0.9.5.1.7
        - ILT (Tikhonov regularization) with SVD Kernel compression
          (1 and 2 dimensions)
        - ``smoosh`` and ``chunk`` deal with axes properly

0.9.5.3 **Current Version**
    upgrade to Python 3 and begin to flesh out documentation

0.9.5.4
    - 0.9.5.4.1
      - ``to_ppm`` should only be a method of inherited class
      - 1.5 and 2.5 D ILT
1.0
    We are working on four major upgrades relative to the 0.9 sequence:

    - Axes as objects rather than a set of separate attributes of nddata.
    - Remove dependence on pytables in favor of h5py.
    - Replace figure lists with “plotting contexts,” which will still
      enable PDF vs. GUI plotting, but will better integrated with Qt and
      more object-oriented in nature
    - Comma-separated indexing to work correctly with all indexing types.
      (0.9.5 requires sequential brackets rather than comma-separated
      indexing for some combined range selections.)

1.0.2
    GUI for setting configuration directories.

    Means for dealing with non-linearly spaced data in image plots
    (0.9.5 auto-detects log spacing in 1D plots,
    but pretends that image plots are linear -- we will implement linear spline
    interpolation algorithm)

1.0.3
    Bruker DSP phase correction for raw data from newer versions of Topspin that is in sync with the code from nmrglue.

1.0.4
    Package a make-less copy of lapack to allow a cross-platform build of density matrix propagation routines.

1.1.0
    Integrate with ACERT NLSL Python package for simulation and fitting of ESR spectra.

1.2.0
    Implement a version of figure list that can be interfaced with Qt.

Installation
============

On Windows with `Anaconda 3.X <https://www.anaconda.com/blog/individual-edition-2020-11>`_,
just run
``conda install -y -c anaconda numpy scipy sympy pyqt pytables matplotlib h5py libpython``
followed by ``conda install -c msys2 m2w64-toolchain`` (the libpython and m2w64-toolchain are only required if you are a developer).
Then (if not a developer) install either via pip (`pip install pyspecdata`) or (if you want to be able to develop or modify the code) follow the `installation for developers <#installation-for-developers>`_ below.

On CentOS7, we've tested
``yum install python-matplotlib python-matplotlib-qt4 python-devel sympy h5py python-tables scipy``
(after running ``yum install epel-release`` to install the EPEL distribution)

On Mac, your python distribution needs to have a working Fortran compiler, since some of the modules use Fortran.

More generally,
these instructions are based on the fact that it's *Highly Recommended* 
that you install the following packages using a good package-management system (conda or linux package manager), rather than relying on `pip` or `setuptools` to install them:

* numpy

* scipy

* sympy

* pyqt

* pytables (in future work, we hope to eliminate dependence on this package)

* matplotlib

* h5py

* The python libraries, and a Fortran compiler.  Under anaconda, these are supplied by `libpython` and `mingw`, respectively.

(If you don't install these packages with your system `pip` will try to install them, and there is a good chance it will fail -- it's known not to work great with several of these; `setuptools` should error out and tell you to install the packages.)

*mayavi*: Mayavi can be used (and gives very nice graphics), but frequently lags behind common Python distros.
Therefore, this package was written so that it doesn't depend on mayavi.
Rather, you can just import ``mayavi.mlab`` and pass it to any figure list that you initialize:
``figlist_var(mlab = mayavi.mlab)``

Installation for developers
---------------------------

Once these are installed,
to install from github, just ``git clone https://github.com/jmfranck/pyspecdata.git`` then move to the directory where setup.py lives,
and do
``python setup.py develop``.
Make sure that this terminates with a successful message, and without any compilation errors.

*Important note for conda on Windows 10:*
For reasons that we don't understand, the Fortran compiler can give odd errors, depending on which terminal you are using to install.
This appears to be Windows' fault, rather than conda's (?).
We highly recommend trying both the Anaconda prompt, as well as the standard dos prompt (press start: type `cmd`) if you experience errors related to compilation.


Notes on compilation of compiled extensions
-------------------------------------------

We recently added a compiled extension that performs non-negative least-squares for regularization (DOSY/Relaxometry/etc.)

Under linux or mac, you should have a gcc and gfortran compiler installed, and should make sure you have libpython for this to work.

Under anaconda on windows, we have run into some trouble sometimes where it gives you an error 127.
We recommend using the normal dos command prompt (cmd) to install pyspecdata, and make sure that your path is set such that
``where gcc`` yields a gcc.exe (NOT .bat) file and ``where python`` yields the anaconda python executable.
(Recent versions of mingw appear to put .bat files in a preferential location
in the path, and these .bat files seem to mess everything up, including
compatibility with the git bash prompt.)

Further installation notes
--------------------------

Upon upgrading from Python 2.X to 3.X, we made some notes in
`conda_upgrade.md <conda_upgrade.md>`_;
this includes some useful (but possibly dated) instructions on how to
implement different environments in anaconda,
how to deal with AppLocker permissions, and Windows permissions generally,
if you run into any of these issues.

Open an issue!
--------------

If you have issues with installing or using pyspecdata, don't hesitate to open
an issue on this page!
