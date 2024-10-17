.. pySpecData documentation master file, created by
   sphinx-quickstart on Sat Mar 12 19:53:01 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pySpecData: compact spectral data processing!
=============================================

pySpecData allows you to deal with multi-dimensional spectroscopy data in an object-oriented fashion.
:ref:`See pySpecData examples here <sphx_glr_auto_examples>`,
and for a simple example applied to 1D data, see
:ref:`here <sphx_glr_auto_examples_basic_example.py>`.

Please note this package is heavily utilized by two other packages that our lab manages on github:

*   `Classes for communicating with instruments <https://github.com/jmfrancklab/FLInst/>`_.

    *   Controls USB instrumentation like oscilloscopes and power meters connected via a USB/GPIB converter.
    *   Includes our Custom SpinCore NMR extension
    *   (Note that the previous two used be separate repositories -- they have been combined to improve maintenance).

*   `ODNP processing scripts <https://github.com/jmfrancklab/proc_scripts/>`_.

The Basics
==========

Our goal is that after you put in a little effort to learn the new way of manipulating data with pySpecData, you can then make the code for processing spectral data that is shorter and also more quickly legible.
PySpecData *automatically* handles the following issues, without any additional code:

*   relabeling axes after a Fourier transformation
*   propagation of errors
*   adding units to plots
*   calculating analytical Jacobians used during least-squares fitting  

To enable this, you work with a pySpecData `nddata` object
(which includes information about dimension names, axis values, errors, and the units)
rather than
working directly with traditional numpy `ndarray` objects.
(pySpecData is built on top of numpy.)

If you have ever worked with arrays in Matlab or Python before, you are
familiar with the additional code needed to convert between index numbers
and axis values.
Using the funny notation of pySpecData, you can do this automatically.

For example, say you have loaded an `nddata` object called `d` and you want to take the time-domain data and:

#.  Fourier transform,
#.  Select out the central 20 kHz in the frequency domain, and finally
#.  Inverse Fourier transform

... all while preserving the correct axes throughout.
That looks like this:

>>> d.ft('t2', shift=True)
>>> d = d['t2':(-10e3,10e-3)]
>>> d.ift('t2')

Note that most pySpecData methods operate *in-place* on the data;
it modifies the data inside the nddata, rather than making a new copy.
This is because we assume that we are progressively optimizing/filtering
our spectral data with each new line of code.
This allows us to quickly work through many operations (like the ft here)
without keeping many copies of a large dataset and with less typing for each operation.
(If you ever truly want to create a copy a dataset, just attach a `.C`)

..  describe updates allowing aliasing (nu2 vs. t2, etc.)

Note that

>>> plot(d)

Quickly generates a publication-quality plot;
it automatically plots the data on the correct axes and includes the units
and the name of dimension (*t2* in this example) and its units along the
*x* axis.

How do I generate an nddata object?
-----------------------------------

We primarily work in magnetic resonance,
so have written wrappers for a few different types of NMR
(nuclear magnetic resonance) and ESR (electron spin resonance)
file formats.
You can use the :func:`pyspecdata.find_file` function to automatically load them as nddata.

Additionally, we have written several classes that allow you to read
nddata objects directly from *e.g.* an oscilloscope.
These are available as separate repositories on github (please contact us for further info).

Finally, you can easily build nddata from standard arrays, as discussed in the section :doc:`about nddata objects <nddata>`.

What can I do with nddata objects?
------------------------------------------

To understand how to manipulate nddata objects, head over to the section
:doc:`about nddata objects <nddata>`.

Contents:
---------

There are many other features of pySpecData that govern the interaction
of ndata objects with plots and, *e.g.* that allow you to generate a nice
PDF laboratory notebook showing every step of your data processing.
These and further details are covered in the various sections of the documentation:

.. toctree::
    :maxdepth: 2

    
    nddata.rst
    fitdata.rst
    modules.rst
    notebook.rst
    figlist.rst
    units.rst
    examples.rst 

.. toctree::
    :maxdepth: 2
    :caption: Example Gallery

    auto_examples/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

