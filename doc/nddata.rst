ND-Data
=======

.. _numpy: http://numpy.org

.. _nddata-summary-label:

.. pulled this from generated/pyspsecdata.nddata.rst
.. add :toctree: to autosummary so stubs are generated
.. remove __init__ method at the top

(This is an introduction.  For detailed info, see :class:`API documentation <pyspecdata.nddata>`.)

The nddata class is built on top of numpy_.
Numpy allows you to create multi-dimensional arrays of data.
The nddata class labels the dimensions (with a short text identifier)
and allows you to associate
axes, units, and errors with the data.
These are carried along with the nddata object,
and *are correctly transformed* when you perform arithmetic operations.
Having dimension labels also makes many common tasks easier.

Dimension labels
----------------

All dimension labels can have a display name (used for printing and plotting)
and one or more short names (used for writing code).

You don't need to keep track of the order of dimensions or align dimensions
during multiplication.
When you do arithmetic with two arrays,
pyspecdata will first reshape the arrays
so that dimensions with the same names are aligned with each other.
Furthermore,
if a dimension is present in one array of an arithmetic operation,
but not the other,
pyspecdata will simply tile the smaller array along the missing dimension(s).

To see how this works, compare the results of

>>> a = nddata(r_[0:4],'a')
>>> b = nddata(r_[0:4],'a')
>>> a*b
[0, 1, 4, 9]
                +/-None
        dimlabels=['a']
        axes={`a':[0, 1, 2, 3]
                        +/-None}

which are arrays of data organized along the *same* dimension,
and

>>> a = nddata(r_[0:4],'a')
>>> b = nddata(r_[0:4],'b')
>>> print a*b
[[0, 0, 0, 0]
 [0, 1, 2, 3]
 [0, 2, 4, 6]
 [0, 3, 6, 9]]
                +/-None
        dimlabels=['a', 'b']
        axes={`a':[0, 1, 2, 3]
                        +/-None,
                `b':[0, 1, 2, 3]
                        +/-None}

which are arrays of data organized along two *different* dimensions.

You can refer to a time dimension, such as `t1`, `t_1`, `t_direct`, *etc.*
    as `f1`, `f_1`, *etc.* in order to retrieve the Fourier transform.
You can set the pairs ...

.. note::
    error propagation for trig functions doesn't yet work;
    t--f not yet done

Item selection
--------------

Pyspecdata offers several types of nddata item selection,
which are outlined below.
A series of these can be combined in a single square bracket, in the usual way.

Numpy Index Slicing
~~~~~~~~~~~~~~~~~~~

Slice by Index
~~~~~~~~~~~~~~

Axis Value
~~~~~~~~~~

To pull a single point along a particular dimension,
you can just use the value of the axis.
The point nearest to the axis will be returned.

>>> d['t2':1.1]

Will return the point (or slice of data) where the t2 axis is closest to 1.1 s.

Axis Ranges
~~~~~~~~~~~

You can specify an inclusive range of numbers along an axis.
For example, to select from 0 to 100 μs along `t2`, you use:

>>> d['t2':(0,100e-6)]

Either value in parentheses can be `None`, in which case, all values to the end of the axis will be selected.

Logical
~~~~~~~

You can use functions that return logical values to select

>>> d[lambda x: abs(x-2)<5]

>>> d['t2',lambda x: abs(x-2)<5]


Fourier domain
~~~~~~~~~~~~~~

>>> d['f2':(-1e6,1e6)]

switches to the frequency dimension (shift by default -- shift should be a property)

Sub-Topics
----------

A selection of the methods noted below are broken down by sub-topic.

.. toctree::
    :maxdepth: 1

    axis_manipulation.rst
    fourier.rst

