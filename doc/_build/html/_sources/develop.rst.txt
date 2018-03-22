.. _nddata_devel:

Overview of ND-Data for Development
===================================

Here we give an overview of the attributes that are used to store
the various nddata components:

:dimlabels:
    An ordered list of the dimension names.
    **In general, the order of other list objects must be given
    dimlabels.**

    If the data is zero-dimensional, this should be set to an
    empty list, in which case `ndshape(...).zero_dimensional`
    will return `True`.

    A dimension name of ``'INDEX'`` is special, and is used to
    explicitly indicate that the axes are unlabeled.
    (see :func:`nddata.__init__`)

:axis_coords:
    The axis labels (*x-*, *y-*, *z-*, and higher-dimensional
    coordinates).  This can be set to:

    - `None` or an empty list, meaning that none of the axes
       have been explicitly labeled.
    - A list of numpy arrays, which correspond to the axis
      labels.

      - One or more of these can be set to `None`, meaning that
        the axis isn't explicitly labeled.
      - Note that the numpy array axis is explicitly allowed to
        be a structured array → see :func:`chunk_auto`.

:axis_coords_error:
    The errors associated with the axis labels.  This can be set
    to the same datatypes as the axis labels themselves.  In
    general, error should only be supplied when coordinates are
    also supplied.

:axis_coords_units:
    The units associated with the axis labels.  Currently, this
    is just a text string, but will be upgraded to a fancier type
    of object that supports conversion.

:data:
    The data itself.
    The number of dimensions must match the length of dimlabels,
    or untowards things will happen!

:data_error:
    The same size as `data` → stores the uncertainty associated
    with data.

    There some future plans to extend this to a
    higher-dimensional object that could store correlated errors.

:data_units:
    The units of the data → again, a string, plans on upgrading.

:other_info:
    This is a dictionary of data that stores any other
    information.
    For example, information pertaining to the state of the axes
    before/after Fourier transformation is stored here,
    as well the acquisition parameters that are read from a data
    file.
    The data is accessed through :func:`nddata.get_prop` and
    :func:`nddata.set_prop`.

:_nosave:
    Attributes that are not saved to file (functions, *etc.*).

Of these, only `data` and `dimlabels` *must* be set → if the
others are missing, they are assumed not to exist.

Note that these attributes are in general **not** intended to be
directly manipulated -- certainly not but the user, and (when
possible) not even by internal functions.
To help in manipulating the data,
various internal helper functions are supplied.

Relevant functions
------------------

:func:`nddata.axn` is used to get the index number associated
with the various dimensions.
For example, ``self.axis_coords_units[self.axn('t2')]`` would
give the units associated with the axis labels for the dimension
called 't2'.

:func:`nddata.mkd` (make dictionary) takes one of the list
objects noted above and transforms it into a dictionary,
where the keys are the elements of dimlabels.
This makes the objects easier to manipulate.  To return the final
result as a list, use :func:`nddata.fld` (flatten dictionary).

