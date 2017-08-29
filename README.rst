If you already know that you want to install -- you can skip to the quick-start_.

===========
pySpecData
===========

Object-oriented Python package for processing spectral data -- or in general, *n*-dimensional data with labeled axes (i.e. *N*-dimensional gridded data or "nddata").
It depends on numpy, which provides very fast manipulations of *N*-dimensional gridded arrays ("ndarray").

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
    - 0.9.5.1.7
      ``to_ppm`` should only be a method of inherited class

    - **0.9.5.1.3**
      XEpr data loaded with dBm units rather than W units

      added ``to_ppm`` function for Bruker files

    - 0.9.5.1.2
      added ability to load power saturation 2D data from Bruker

    - 0.9.5.1.1
      Some important debugging, and also added `pyspecdata.ipy` → executing the following at the top of a jupyter notebook:

        .. code-block:: python

            %pylab inline
            %load_ext pyspecdata.ipy

      will cause nddata to "display" as labeled plots.

0.9.5.2
    Comma-separated indexing to work correctly with all indexing types.
    (0.9.5 requires sequential brackets rather than comma-separated indexing for some combined range selections.)

0.9.5.4
    GUI for setting configuration directories.

    Means for dealing with non-linearly spaced data in image plots
    (0.9.5 auto-detects log spacing in 1D plots,
    but pretends that image plots are linear -- we will implement linear spline
    interpolation algorithm)

0.9.5.5
    Bruker DSP phase correction for raw data from newer versions of Topspin that is in sync with the code from nmrglue.

0.9.5.8
    Package a make-less copy of lapack to allow a cross-platform build of density matrix propagation routines.

1.1.0
    Integrate with ACERT NLSL Python package for simulation and fitting of ESR spectra.

1.2.0
    Implement a version of figure list that can be interfaced with Qt.


Installation Notes
==================

*Highly Recommended:* 
Install the following packages using a good package-management system (conda or linux package manager), rather than relying on `pip` or `setuptools` to install them:

* numpy

* scipy

* sympy

* pyqt

* pytables

* matplotlib

* h5py

For example, on Windows with `Anaconda 2.7`_.
-- just run
``conda install numpy scipy sympy pyqt pytables matplotlib h5py``.

(If you don't install these packages with your system `pip` will try to install them, and there is a good chance it will fail -- it's known not to work great with several of these; `setuptools` should error out and tell you to install the packages.)

*mayavi*: Mayavi can be used (and gives very nice graphics), but frequently lags behind common Python distros.
Therefore, this package was written so that it doesn't depend on mayavi.
Rather, you can just import ``mayavi.mlab`` and pass it to any figure list that you initialize:
``figlist_var(mlab = mayavi.mlab)``

For compiled extensions
-----------------------

All compiled extensions are currently stripped out, but will be slowly
    added back in.

If you are installing from github (or generally using setuptools -- *i.e.* ``python setup.py install``
or ``python setup.py develop``).

If you are on windows, you will need some additional packages to enable compilation:

* libpython

* unxutils

* mingw

The last two are specific to Windows, and provide things like the ``gcc`` and ``gfortran`` compiler, as well as ``make``.

Installation for developers
---------------------------

(Once these are installed,
to install from github, just ``git clone https://github.com/jmfranck/pyspecdata.git`` then move to the directory where setup.py lives,
and do
``python setup_paramset.py develop``
followed by
``python setup.py develop``)

Quick-Start
===========

To get started with this code:

1. Install a good Python 2.7 distribution

   * On Windows or MacOS: `Anaconda 2.7 <https://www.continuum.io/downloads>`_.  When installing select "install for all users."

2. Install libraries that pyspecdata depends on. (If you're interested in why you need to do this first, see installation notes below.)

   * On Windows or MacOS: in the Anaconda Prompt, run ``conda install numpy scipy sympy pyqt pytables matplotlib h5py``.

   * For Mac, you can also use homebrew.
     Note that, in the current version python is renamed to `python2`,
     and `pip` to `pip2`.
     Most packages can just be installed with `pip2` under homebrew.
     If you want HDF5 functionality, you will need to run `brew tap homebrew/science` followed by `brew install hdf5`.

   * On Linux, just use your package manager (``aptitude``, ``yum``, *etc.*) to install these libraries.

3. Install `pyspecdata`: ``pip install pyspecdata``

4. Set up directories -- create a file in your home directory
   called ``_pyspecdata`` (Windows  -- note the underscore)
   or ``.pyspecdata`` (Mac or Linux).
   Here is an example -- you can copy and paste it as a starting point:

   ::

        [General]
        data_directory = c:/Users/yourusername/exp_data
        notebook_directory = c:/Users/yourusername/notebook

   Note that any backslashes are substituted with forward slashes.
   Also note that you will
   need to change the directories to refer to real directories that already
   exist or that you create on your hard drive (see below).
   Note that on Windows, you can use notebook, *etc.* to create this file,
   but it cannot have a .txt, *etc.* `extension <http://www.wikihow.com/Change-a-File-Extension>`_.

   * Where is my "home directory"? (Where do I put the `_pyspecdata` file?)

       * On Windows, your home directory is likely something like
         ``C:\Users\yourusername``.
         You can access your home directory by opening any file folder window, and
         starting to type your name in the address bar -- it's the first folder that shows up
         underneath.

       * On MacOS and Linux, it's the directory indicated by ``~``.  On Linux,
         this typically expands to ``/home/yourusername``.

   * What are these directories? → You can either create them or point to existing directories.

       * ``data_directory`` must be set.  It is a directory, anywhere on the
         hard drive, where you store all your raw experimental data.  It must
         contain at least one subdirectory -- each subdirectory stores
         different "experiment types," typically acquired on different instruments
         (*e.g.* you might have subdirectories named ``400MHz_NMR``,
         ``500MHz_NMR``, ``95GHz_ESR``, and ``Xband_ESR``).

           * Data is assumed to be **unpacked** (*i.e.* as it is on the spectrometer -- not in .zip or .tgz files)

           * If you're setting up a lab, you might want to separately sync each different
             experiment type folders using `seafile <https://www.seafile.com/en/home/>`_.

             Or you can sync the whole data directory with dropbox.

       * If set, the ``notebook_directory`` is intended to contain latex
         files with embedded python code, as well as some processed
         output.

   * *Do not* use quotes to surround the directory name.  Even if it contains
     spaces, do not use quotes, and do not escape spaces with backslashes.

   * Note that on Windows, your desktop folder is typically in ``C:\Users\yourusername\Desktop``

   * Why do I need to do this?

       * Setting this configuration allows you to move code between different
         computers (*e.g.* a spectrometer computer, a desktop, and a laptop),
         and re-use the same code, even though the locations of the files are
         changing.  This should work even across different operating systems.

       * It specifically enables functions like ``find_file(...)``,
         ``get_datadir(...)``, *etc.* that can search the data directory for a
         file name matching some basic criteria.
         You should always use these to load your data,
         and *never* use the absolute path.

       * The GUI tool that will allow you to set up ``_pyspecdata`` by pointing
         and clicking has not yet been set up.
