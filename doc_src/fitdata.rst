the fitdata class
=================

This is a child class of nddata used for fitting.

For old fitdata documentation
(obsolete, for development reference) see :ref:`fitdata_old`

creating new types of fitdata modules
-------------------------------------

There is a base class called “fitdata” that defines the basic routines
necessary for fitting. Currently, the fitdata class only supports
fitting along one dimension, so before constructing a new class, one
must first choose what dimension they will be fitting along.

.. todo::
    Of
    course, for multidimensional data, the fit will be repeated along the
    dimensions that are not the fit dimension. see how easy it would be to
    allow more than one dimension

To fit a new type of function, one simply creates a new type of class
that *inherits* from the fitdata class. We override all the methods that
have to do with the definition of the functional format. These are
defined in the first section, where we build up an example for fitting a
general :math:`T_1` recovery curve. This example should be used as a
starting point for making new fit classes. Then, we can make instances
of the new class, and use their methods (described in the subsequent
section) next.

.. todo::
    the option block :no-inherited-members: doesn't work -- not sure how to modify class.rst
    I put a template from stackexchange inside _templates
    -- see https://stackoverflow.com/questions/28147432/how-to-customize-sphinx-ext-autosummary-rst-template
    on how to use it

    then, I need to link to or include generated/pyspecdata.core.fitdata.rst

.. comment
    .. autosummary::
        :toctree: generated
        ~fitdata

.. currentmodule:: pyspecdata.core

.. autoclass:: fitdata
    :members:

