class axis_collection(object):
    """A collection of :class:`axis` objects.

    Designed so that an instance of `axis_collection` is an attribute of `nddata` called `axes`,
    which behaves like a dictionary whose keys are the `dimlabels` of the `nddata` object,
    and whose values are :class:`axis` objects.

    Used to make sure that no axis names or aliases are duplicated.

    You can add axes to the collection in any of the following ways, where `example_nddata` is an nddata instance.
    (Remember that all `nddata` instances have an attribute `axes` of type `axis_collection`).

    building a new axis
        `example_nddata.axes['t2'] = ax_[0:1.2:100j]`
        or
        `example_nddata.axes['t2'] = ax_[0:1.2:0.01]`
        (uses the same notation as numpy `r_[…]`)

        this takes the place of `labels` or `setaxis` in old versions of pyspecdata.

    associating an existing axis
        `example_nddata.axes += existing_axis` `existing_axis` **must** have a
        name or alias that matches one of `example_nddata.dimlabels`.

    Attributes
    ----------
    dimlabels: list
        This is the same `dimlabels` attribute as the instance of the parent class.
    names_used: ChainMap
        Stores a list of all the names and aliases used by the `axis` objects
        that are contained in the collection,
        as well as the axes for any conjugate domains.
        since these need to be unique.

        This is simply
        `ChainMap(ax1.references,ax2.references,…,etc.)`
    """
    def __init__(self,dimlabels):
        self.dimlabels = dimlabel
    def rename(self,oldname,newname):
        """Rename an axis. If `oldname` is the preferred name of the axis,
        also go into dimlabels, and change the name
        (since dimlabels is the same list used by the `nddata` that
        contains the collection, it will update the dimlabel there as
        well)"""
        if oldname in self.dimlabels:
            idx = self.dimlabels.index(oldname)
            self.dimlabels[idx] = oldname
            return self
        else:
            raise ValueError("here add code to go look at the aliases")
class _ax_class_maker(object):
    def __getslice__(self,inp_slice):
        pass
ax_ = _ax_class_maker()
class nddata_axis(object):
    """The axis that gives the list of coordinates along a particular dimension.

    .. todo::
        There is no actual code here -- this is a proposal for the new axis class

    Internally uses the minimum number of variables to store information about the axis.

    Also includes information about the chosen location (alias) in infinitely periodic domains.
    This is analogous to the concept of using `fftshift` in matlab or traditional numpy,
    but more general.

    The `nddata_axis` has overloading routines to deal with the following operations like a standard numpy array
    (`example_axis` below is an instance of `nddata_axis`)

    indexing
        >>> retval = example_axis[1]

        returns the second value along the axis

    slicing
        >>> retval = example_axis[0:20:5]

        returns every fifth value from zero up to, but not including, 20

    nddata-like slicing
        >>> retval = example_axis[(0,5.5)]

        returns everything where the values of the axis coordinates are between 0 and 5.5 (inclusive)

        >>> retval = example_axis[(0,None)]

        returns everything where the values of the axis coordinates are 0 (inclusive) or above

    multiplication
        >>> retval = example_axis * b

        or

        >>> retval = b * example_axis

        if ``b`` is a numpy array
            will return another numpy array

        if ``b`` is an nddata
            will return another nddata
            -- note that this replaces the previous use of ``fromaxis``.

    addition + subtraction + division
        same rules as multiplication

    argument of a function
        >>> retval = exp(example_axis)

        (or ``sin``, ``cos``, *etc.*)
        returns another axis object.
        Note that this just amounts to setting the transf_func attribute, below.

        If ``self.multiplier`` is set to a complex number,
        specialized routines are always used
        (*e.g.* ``exp`` can be calculated more efficiently, *etc.*)

    interpolation (``@``)
        >>> retval = b @ example_axis

        Here, ``b`` must be an nddata,
        and ``example_axis`` must have a ``name`` matching one of the dimension labels of ``b``.

        ``retval`` will consist of ``b`` interpolated onto the new axis.

        Note that while ``@`` is typically used for matrix multiplication,
        it is NOT here.

    Attributes
    ----------
    size: long
        the length of the axis
    dx: float
        Step size multiplying the base array.
        For a non-uniform array,
        if possible,
        divide by the smallest step size,
        then multiply by a number that will
        `convert the resulting floats to integers <https://stackoverflow.com/questions/44587875/find-common-factor-to-convert-list-of-floats-to-list-of-integers>`_.
    start: float
        determines the starting value of the axis:
        >>> self.start+self.dx*r_[0:self.size]

    names: list of strings or sympy variables
        Names for this dimension that this axis is used to label, in order of preference.
        The first name is the "preferred" name,
        and all subsequent names are "aliases".
        For example,
        you might want to have a nicely named
        :math:`B_0` (stored as ``$B_0$`` or a sympy variable)
        to describe your axis
    domains: OrderedDict
        The keys correspond to a list of allowed transformations.
        Currently these are (future plans for ``(I)LT``, ``(I)NUS``, ``(I)Wavelet``)

        - ``'FT'``
        - ``'IFT'``

        These are the names of transformations that have previously been applied
        (or can be applied, though the list doesn't need to be comprehensive in that case)
        to the ``nddata`` object that the ``nddata_axis`` is being used to label.
        ``I...`` must **always** stand for "inverse"; on application of a transformation,
        the new ``axis`` object that is generated must have a `domains` attribute
        with the opposite (``I`` removed or added) transformation listed.

        The values are `axis` objects that label the data in the conjugate domains (after the transformation has been applied).

        For example,
        on application of the `nddata` :func:`nddata.ft` method,
        the data will be labeled with an axis that has a `domains` attribute with a key containing `IFT`.
        The value of that key will point to the `axis` object of the data *before* transformation,
        and will be used in the even of a call to :func:`nddata.ift`.

        This makes the `get_ft_props` and `set_ft_props` of older versions of `nddata` obsolete.
    multiplier: complex, default None
        this is *only* used in the event that 
        the axis is subjected to arithmetic involving a complex number
        it changes the way that the axis acts as an argument to various functions (especially `exp`)
    transf_func: function or (default) None
        **this and following attributes pertain only to non-uniform (non-linear) axes**
        a function that is applied to a uniformly spaced axis to achieve non-uniform spacing
        (for example, an axis with `log10` spacing).
        If this is set, the axis is constructed as

        >>> self.transf_func(self.start+self.dx*r_[0:self.size])

    uneven_steps: int or float (default non-existent)
        if this attribute exists, it must be an array of length `self.size`
        and determines the axis values as follows:
        `self.offset+self.dx*cumsum(self.uneven_steps)`
    uneven_step_coords:
        if `self.uneven_steps` exists, this stores the value of `cumsum(self.uneven_steps)`
    """
    def __init__(self, *args):
        """Either instantiates an empty array (`nddata_axis()`) for modification or (`nddata(inp_array)`) analyzes a numpy array for conversion.

        To construct with the analog of `r_[...]`, see the helper class `ax_[...]`;
        because of the nature of the axis storage `ax_[...]` will be far more efficient than passing `inp_array`.

        Parameters
        ----------
        inp_array: ndarray
            A standard numpy array to be converted to an axis.
            The initialization routine will check to see (in order):

            -   is the axis uniformly spaced?
            -   is the axis a uniformly spaced axis that is transformed by `log10`, `log`, `sin`, `cos`, *etc.*
        """
    def to_array(self):
        "returns the axis as a standard numpy array"
    @property
    def references():
        """returns OrderedDict of all names and aliases such that keys all point to the current instance (`self`)

        the idea is that this should be placed in a `ChainMap` object to be used by the :class:`axis_collection` class that contains the axis.
        """
