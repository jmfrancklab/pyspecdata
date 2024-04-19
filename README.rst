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

`Roadmap`_.

.. _Roadmap: changelog.rst

Installation
============

**Important note:**
the package ships Fortran-based extensions that are used to provide fast ILT methods.
We believe this is a useful feature.
Unfortunately,
while the instructions below work for most cases,
not everyone's system is set up equally well for Fortran compilation.
If you experience difficulties, please don't hesitate to reach out to us at jmfranck [at] syr.edu;
we would be happy for the opportunity to test distribution on new platforms!
In all situations, note that this is a development library that works very well
in our hands -- we are happy to hear from you and work with you to try to
broaden its applicability!

On **Windows** with `Anaconda 3.X <https://www.anaconda.com/blog/individual-edition-2020-11>`_,
just run ``conda install -y -c anaconda numpy scipy sympy pyqt pytables matplotlib h5py libpython`` followed by ``conda install -y m2w64-toolchain`` (the libpython and m2w64-toolchain are for building compiled extensions such as the ILT).
Then follow the `installation for developers <#installation-for-developers>`_ below. We have a package on pip, but it currently lags behind the github repo.

On **CentOS7**, we've tested
``yum install python-matplotlib python-matplotlib-qt4 python-devel sympy h5py python-tables scipy``
(after running ``yum install epel-release`` to install the EPEL distribution).  Then follow the `installation for developers <#installation-for-developers>`_ below. 

On **Debian** (should also work for **Ubuntu**),
we've tested
``sudo apt-get install -y python3 python3-matplotlib libpython3.7 python3-dev python3-sympy python3-h5py python3-tables python3-scipy python3-setuptools gfortran pip``.  Then follow the `installation for developers <#installation-for-developers>`_ below. 

On **MacOS**, if you want to install as a developer your python distribution needs to have a working Fortran compiler, since some of the modules use Fortran.
We have tested ``conda install -c conda-forge fortran-compiler``, followed by
``conda install -y -c anaconda numpy scipy sympy pyqt pytables matplotlib h5py``.
However *due to a problem with more recent versions of MacOS/xcode*, you need to modify ``setup.py`` to tell it where to find the system libraries.
At about line 27, you need to add something like following as a keyword arg for the `Extension` function:
``library_dirs = ["/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"],``
(we recommmend just using keyword completion to find a similar directory that exists).
(Feel free to contact us if you have issues with this or would like to test deployment on pip for a Mac).

**More generally,**
these instructions are based on the fact that it's *Highly Recommended* 
that you install the following packages using a good package-management system (conda or linux package manager), rather than relying on `pip` or `setuptools` to install them:

* numpy

* scipy

* sympy

* pyqt

* pytables (in future work, we hope to eliminate dependence on this package)

* matplotlib

* h5py

* lmfit  

* The python libraries, and a Fortran compiler.  Under anaconda, these are supplied by `libpython` and `mingw`, respectively.

* If you plan on building the documentation, you also want `sphinx_rtd_theme sphinx-gallery`

(If you don't install these packages with your system `pip` will try to install them, and there is a good chance it will fail -- it's known not to work great with several of these; `setuptools` should error out and tell you to install the packages.)

*mayavi*: Mayavi can be used (and gives very nice graphics), but frequently lags behind common Python distros.
Therefore, this package was written so that it doesn't depend on mayavi.
Rather, you can just import ``mayavi.mlab`` and pass it to any figure list that you initialize:
``figlist_var(mlab = mayavi.mlab)``

Installation for developers
---------------------------

To install pySpecData from github, just ``git clone https://github.com/jmfranck/pyspecdata.git``. Then switch over to the anaconda prompt and move to the directory where setup.py lives (root directory of repository),
and type
``python setup.py develop``.
Make sure that this terminates with a successful message, and without any compilation errors.  In particular:

- If it gives an error about permissions (will happen for a system-wide anaconda install), you need to load the anaconda prompt as admin (right click and run as administrator).
- Near the end (above EXT compiler optimization) it should tell you that you can run `pyspecdata_dataconfig`.  You should do this, unless you've installed pyspecdata before on the computer you are working at.

Important notes for conda on Windows:

- **Warning** Before running the installation for developers, you must
  first check that the output of ``conda info`` on your git bash terminal
  matches the output of your anaconda prompt.
- For reasons that we don't understand, the Fortran compiler can give odd
  errors, depending on which terminal you are using to install.  This
  appears to be Windows' fault, rather than conda's (?).  We highly
  recommend trying both the Anaconda prompt, as well as the standard dos
  prompt (press start: type `cmd`) if you experience errors related to
  compilation.
- If you want to build the documentation, run:
  `conda install -y -c conda-forge sphinx_rtd_theme sphinx-gallery`

Data File Management
====================

pySpecData is designed to run the same script on different computers,
where the required data files might be stored in different paths
on the different computers.

The basic strategy is that you enter information on how to find your
files in the `_pyspecdata` config file (typically this is only required once,
at setup),
then the `find_file` and `search_filename` functions can use this info
to find your files.

Setting up your _pyspecdata configuration file
----------------------------------------------

Part of the pySpecData package is the datadir module, allowing the user to run the same code on 
different machines - even thought the location of the raw spectral data might change. 
This is controlled by the ``~/.pyspecdata`` (unix-like) or ``~/_pyspecdata`` (windows) config file,
which looks like the following.

::

    [General]
    data_directory = /home/jmfranck/exp_data
    qesr conversion = 162.66
    qesr diameter = 0.704
    qesr q = 4700

    [ExpTypes]
    odnp_nmr_comp/odnp = /home/jmfranck/exp_data/NMR_comp/ODNP

    [mode]
    figures = standard

    [RcloneRemotes]
    nmr_comp/odnp = jmf_teams:General/exp_data/NMR_comp/ODNP/

The ``General`` section points to the directory with the datasets of interest whether that is the
direct path to the drive with the datasets or if you prefer Rclone, this ``data_directory``
points to your local folder of datasets.
(This is also a good spot to include, *e.g.* proportionality constants for
QESR, which we have done here, and which are utilized in the `proc_scripts`
repo.)

The ``ExpTypes`` section gives the various locations to 
folders containing the appropriate data sets - either pointing to the
cloud storage or pointing to the local directory your rclone adds files to.
So when you call ``odnp_nmr_comp/odnp`` this will point
to the actual location, ``/home/jmfranck/exp_data/NMR_comp/ODNP``

Note that it's possible to point the different `exp_type` directly to shared drives,
pySpecData also offers a (we think superior) method that downloads local copies
of files on-demand using `rclone <https://rclone.org/>`_.
Obviously, you need to install rclone and add it to your path to do this (see next subsection).
Rclone is an amazing tool that can be configured to talk to virtually any type of cloud storage
(Google Drive accounts, OneDrive and SharePoint accounts, etc.)

Inside the ``RcloneRemote`` section, each key/variable points to a properly configured remote that
was set up with `rclone <https://rclone.org/>`_--
e.g., ``jmf_teams`` here is a properly configured  remote that shows up
in response to the shell command ``rclone config``.
*Note:* as you require datasets from other folders you will need to make new folders locally to match
for Rclone.
You will receive error messages that guide you to do this, and you should follow them.
For example, if you required a dataset from ``exp_data/francklab_esr/alex`` you
will need to go into your local ``exp_data`` folder and add a new folder called ``francklab_esr/alex``

Setting up Rclone
-----------------

To get set up with Rclone, download Rclone and follow the documentation which should include
running the command ``rclone config`` enabling you to set up the location and name of the cloud
drive you wish to pull from.
The documentation of rclone is pretty straightforward and can walk
you through this. 
If you are at an academic institution, we highly recommend asking your IT
department for a protocol for connecting rclone to your cloud storage of
choice.

Notes on compilation of compiled extensions
===========================================

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
==============

If you have issues with installing or using pyspecdata, don't hesitate to open
an issue on this page!
