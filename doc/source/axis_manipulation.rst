Axis Manipulation Functions
===========================

.. autosummary::
    pyspecdata.axis_manipulation

.. currentmodule:: pyspecdata.axis_manipulation

.. rubric:: Member Functions
.. autosummary::
    :toctree: generated

    inhomog_coords.inhomog_coords
    shear.linear_shear
    register_axis.register_axis
    secsy.secsy_transform_manual

.. comment
    here, I have generated stub pages, and after a bit of experimentation, figure out that I do need to point to the actual source of the functions
    I can generate the above list with
    grep ~/notebook/pyspecdata/pyspecdata/axis_manipulation -rie "^def" | sed 's/.*\/\([^\/]*\)\.py:def *\([^(]*\).*/\1.\2/'

