Examples
========
Here are some examples to help navigate the different functions and tools included in pyspecdata

.. note::
    *Alex*: if we can figure out what, specifically, we need to import from
    pylab to get the following to work, then we can probably get the figurelist
    code to work as-is -- as you see from my comments I can't just import the
    functions that I need from pylab -- I'm guessing there some other function
    that is called by sphinx that also needs to by imported

.. plot:: ../examples/matplotlib_test.py
    :include-source:

Here is an example of the use of fitdata and setting a initial guess then plotting simultaneously.

.. plot:: ../examples/plot_fit_fake_data.py
    :include-source:

.. toctree::
    :maxdepth: 2
    :caption: example of fitdata
    auto_examples/plot_fit_fake_data
