HFSS-derived .fld data
======================

ACERT functions for processing HFSS-derived data.
These are intended to be paired with some IronPython scripts that can be run to export .fld files from within Ansys HFSS.

.. comment made the following file
    .. currentmodule:: pyspecdata
    .. autosummary::
        :toctree: generated
        (linebreak here)
        acert_hfss
    then ran
    sphinx-autogen -o generated hfss.rst
    with the import statements commented out
    then pulled the result in here and edited it
    by specifying members, I eliminate the need to comment out

.. currentmodule:: pyspecdata.acert_hfss

.. rubric:: Functions

Most likely, you will want to use :func:`load_fields <pyspecdata.acert_hfss.load_fields>`

.. autosummary::
    :toctree: generated

    construct_axes_from_positions
    contour_power
    gaussian_over
    load_fields
    load_hfss_scalar
    load_hfss_vectors
    w_index

