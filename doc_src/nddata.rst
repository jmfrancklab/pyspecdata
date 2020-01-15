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
These attributes
*are correctly transformed* when you perform arithmetic operations
or Fourier transforms,
and are used to automatically format plots.

.. todo::
    also implement with blocks for working in the Fourier domain

Very importantly, most pyspecdata functions are designed to operate on the data *in place*, meaning that
rather than doing something like:

>>> data = data.mean() # numpy

you simply do:

>>> data.mean('axisname')

and ``data`` is modified from here out.
We do this because we are typically processing multidimensional datasets that consist of many points,
and we subject them to a series of steps to process them.
We don't want to use up memory with lots of copies of the data.
Also, this allows us to string together several operation, e.g.:

>>> data.ft('axis1').sum('axis2')['axis1':(0,3)].mean('axis3')

So, while this general setup is different than the standard numpy setup,
*etc.*, it should lead to you writing more efficient code, with less variables
to keep track of, and generally also leads to far more compact code,
and you should probably not try to bypass it by creating copies of your data.

In rare circumstances, you *need* to create a copy of your data
(I think there is only one case where this is true:
when you
need to process the same data in parallel in two different ways to generate two results,
and then perform math that uses both results).
For these cases, you can use an easy shortcut for the copy method: `C` (as in `data.C`).

**say something about slicing and copies here**

Note that, if all else fails,
you can always use numpy directly:
*The `data` attribute of an nddata object is just a standard numpy array.*
For example, if you have an nddata object called ``mydata``,
then ``mydata.data`` is just a **standard numpy array**.
For beginning users,
this can be a very helpful fall back.
But, note that -- in most cases -- it should be more beneficial if you don't
directly access the `data` attribute.

Below, we outline how you can use
dimension labels to make code more legible and make many common tasks easier.
Then, we note how slicing operations are different (and easier) for nddata than for standard numpy arrays.
Finally, we outline several classes of methods by sub-topic.

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
during an arithmetic operation,
if a dimension is present in one array
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

.. warning::
    error propagation for trig functions doesn't yet work;
    automatic t *vs.* f naming not yet implemented.

Item selection and slicing
--------------------------

Pyspecdata offers several different synataxes 
for item selection and slicing.
These fall into two main categories:
*numbered indexing and slicing* (which will be familiar to most users of python)
and *axis-based indexing*a.
These can be combined in a single square bracket, separated by commas.

Axis-Based Indexing
~~~~~~~~~~~~~~~~~~~

Axis-based Indexing
``````````````````

To pull a single point along a particular dimension,
you can just use the value of the axis.
The point nearest to the axis will be returned.

>>> d['t2':1.1]

Will return the point (or slice of data) where the t2 axis is closest to 1.1 s.

Axis Ranges
``````````````````

You can specify an inclusive range of numbers along an axis.
For example, to select from 0 to 100 μs along `t2`, you use:

>>> d['t2':(0,100e-6)]

Either value in parentheses can be `None`, in which case, all values to the end of the axis will be selected.

Logical
``````````````````

You can use functions that return logical values to select

>>> d[lambda x: abs(x-2)<5]

>>> d['t2',lambda x: abs(x-2)<5]

When this is done, nddata will check to see if slices along any dimension are uniformly missing.
If they are, the dataset will be trimmed to remove them.

When the deselected data are scattered throughout, a mask is used instead.

Axis-Based Indexing and Slicing in the Fourier Domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data can be manipulated not only in the direct domain,
but also in the Fourier conjugate domain.
This is achieved in one of two ways → by changing the name
of the axis,
and by explicitly specifying the domain.


Changing the name of the axis
`````````````````````````````

.. warning::
    this feature is planned, not yet implemented.

The property ``FT_regexp`` is a dictionary that controls
By default, it is initially set to ``{'\<t':'\<f','\<T':'\<F'}``,
meaning that t and f at the beginning of word will be converted into each other.
This means that the following command:

>>> d['f2':(-1e6,1e6)]

would switch something with a dimension named `t2` to the frequency dimension
(shift by default -- shift should be a property).


Explicitly Changing the Domain
```````````````````````````````

Alternatively, a slice of the form
``dimname:'t/f':(range)``
or
``dimname:'u/v':(range)``
can be used.
Here, 't' and 'f' always correspond to the "original" and "conjugate" domains, respectively.
Since, sometimes this might lead to confusion, 'u' and 'v' serve the same respective functions.

Therefore this command performs the same function as the last line of code:

>>> d['t2':'f':(-1e6,1e6)]

Frequency Slices Outside Original Range
```````````````````````````````````````

Numbered Indexing and Slicing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Basic Examples
--------------

* Apply a filter (fromaxis).
* Slicing.
* Aliasing of FT.
  

Methods by Sub-Topic
--------------------

It's important to note that, in contrast to standard numpy,
nddata routines are designed to be called as methods,
rather than independent functions.
Also, these methods **modify the data in-place** rather than returning a copy.
For example, after executing
``d.ft('t2')``, the object ``d`` will contain the *Fourier transformed* data.
There is not need to assign the result to a new variable.
Alternatively, the property ``C`` offers easy access to a copy:
``a = d.C.ft('t2')`` leaves ``d`` alone, and returns the FT as a new object
called ``a``.

This encourages a style where methods are chained together, *e.g.* ``d.ft('t2').mean('t1')``.

A selection of the methods noted below are broken down by sub-topic.

.. toctree::
    :maxdepth: 1

    axis_manipulation.rst
    fourier.rst

