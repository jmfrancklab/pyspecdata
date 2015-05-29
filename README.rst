===========
pySpecData
===========

Object-oriented module for processing spectral data -- or in general, *n*-dimensional data with labeled axes (i.e. n-Dimensional gridded data like an HDF SDS).  If you are working in a lab developing new spectroscopic methodologies, then this lab is definitely for you.  If you deal with multi-dimensional data of some other form, then it's likely for you.

Please see `the homepage on github <http://jfranck.github.com>`_.

* Labeled axes allow one to manipulate datasets (potentiall with different dimensions) without having to explicitly keep track of what the different dimensions correspond to.  Code becomes more legible.  Also, tiling, direct product, and griding functions become obsolete.

* Fourier transformation with automatic manipulation of axes.

* Automatic error propagation.

* Reading and writing to HDF5.

* The code is written so that it can be integrated into a latex lab notebook.  The same code that generates pop-up windows with plots from the command line can be embedded into a Latex document. Extension to other output formats, such as HTML or markdown, should be relatively straightforward.

NMR/ESR specific
=========
Because it was written primarily for NMR data, it also includes:

* Routines for reading commercial raw data into objects with all relevant information.

* A basic compiled routine for propagating density matrices that can be used to predict the response to shaped pulses.

