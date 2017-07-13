===========
pySpecData
===========

Installation
============

*Highly Recommended:* 
If installing with `pip` install the following packages with your system, rather than relying on `pip` to install them:

* numpy

* scipy

* sympy

* pyqt

* pytables

* matplotlib

* h5py

For example, on windows with `Anaconda 2.7 <https://www.continuum.io/downloads>`_.
-- just run
``conda install numpy scipy sympy pyqt pytables matplotlib h5py``.

(If you don't install these packages with your system `pip` will try to install them, and there is a good chance it will fail -- it's known not to work great with several of these).

*mayavi*: Mayavi can be used (and gives very nice graphics), but can be difficult to install.
This package doesn't depend on mayavi.  Rather, you import it and pass it to the figure list that you initialize:
``figlist_var(mlab = mayavi.mlab)``

Version Notes
=============

Note that version is currently 0.9.4 -- currently intended just for collaborators, *etc.*
A general-use version 1.0.0 is planned within a year.
*(Note that the email currently linked to the PyPI account is infrequently checked --if you have interest in this software, please find J. Franck's website and contact by that email.)*

Object-oriented Python package for processing spectral data -- or in general, *n*-dimensional data with labeled axes (i.e. *n*-Dimensional gridded data like an HDF SDS).  If you are working in a lab developing new spectroscopic methodologies, then this package is definitely for you.  If you deal with multi-dimensional data of some other form, then it's likely for you.

* Labeled axes allow one to manipulate datasets (potentially with different dimensions) without having to explicitly keep track of what the different dimensions correspond to.  Code becomes more legible.  Also, tiling, direct product, and griding functions become obsolete.

* Fourier transformation with automatic manipulation of axes.

* Automatic error propagation.

* Reading and writing to HDF5.

* The code is written so that it can be integrated into a LaTeX lab notebook.  The same code that generates pop-up windows with plots from the command line can be embedded into a Latex document. Extension to other output formats, such as HTML or markdown, should be relatively straightforward.

More detailed web documentation will be coming soon.

NMR/ESR specific
================

Because it was written primarily for NMR data, it also includes:

* Routines for reading commercial raw data into objects with all relevant information.

* (Not yet in packaged version) A basic compiled routine for propagating density matrices that can be used to predict the response to shaped pulses.

