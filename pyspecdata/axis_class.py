"""Minimal axis utilities used by the test-suite.

This module provides a very small subset of the functionality described in
the documentation.  The implementation only aims to support the behaviour
exercised in the accompanying unit tests and is **not** a drop-in
replacement for the full featured version.
"""

from __future__ import annotations

from collections import ChainMap, OrderedDict
from typing import Iterable, List, MutableMapping

import numpy as np


class axis_collection:
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

    def __init__(self, dimlabels: Iterable[str]):
        self.dimlabels: List[str] = list(dimlabels)
        self._axes: MutableMapping[str, nddata_axis] = {}
        self.names_used: ChainMap = ChainMap()

    def __getitem__(self, key: str) -> "nddata_axis":
        return self._axes[key]

    def __setitem__(self, key: str, value: "nddata_axis") -> None:
        if key not in self.dimlabels:
            raise KeyError(key)
        if key in self._axes:
            old = self._axes[key]
            self.names_used.maps = [
                m for m in self.names_used.maps if not m or list(m.values())[0] is not old
            ]
        if not getattr(value, "names", []):
            value.names = [key]
        for name in value.names:
            if name in self.names_used:
                raise ValueError("duplicate axis name")
        self.names_used.maps.append({n: value for n in value.names})
        self._axes[key] = value

    def __iadd__(self, axis: "nddata_axis") -> "axis_collection":
        for name in axis.names:
            if name in self.dimlabels:
                self[name] = axis
                return self
        raise ValueError("axis name does not match any dimlabel")

    def rename(self, old: str, new: str) -> "axis_collection":
        """Rename an axis. If `oldname` is the preferred name of the axis,
        also go into dimlabels, and change the name
        (since dimlabels is the same list used by the `nddata` that
        contains the collection, it will update the dimlabel there as
        well)"""
        if old not in self.dimlabels:
            raise ValueError(f"{old!r} not a current dimlabel")
        idx = self.dimlabels.index(old)
        self.dimlabels[idx] = new
        if old in self._axes:
            ax = self._axes.pop(old)
            ax.names = [new if n == old else n for n in ax.names]
            self._axes[new] = ax
        return self


class _ax_class_maker:
    def __getitem__(self, slc):
        if not isinstance(slc, slice):  # pragma: no cover - defensive
            raise TypeError("axis builder expects slice syntax")
        start = 0.0 if slc.start is None else float(slc.start)
        stop = float(slc.stop)
        step = slc.step
        if isinstance(step, complex):
            arr = np.linspace(start, stop, int(abs(step)))
        else:
            arr = np.arange(start, stop, float(step))
        return nddata_axis(arr)


ax_ = _ax_class_maker()


class nddata_axis:
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

    __array_priority__ = 1000

    def __init__(self, arr: Iterable[float] = ()):  # type: ignore[assignment]
        """Either instantiates an empty array or analyzes ``arr`` for conversion.

        To construct with the analog of ``r_[...]``, see the helper class
        :data:`ax_`, which is more efficient than passing a full array.

        Parameters
        ----------
        arr
            Iterable of coordinate values used to populate the axis.  If
            omitted, an empty axis is created for later modification.
        """
        self.data = np.array(arr, dtype=float)
        self.names: List[str] = []

    # ------------------------------------------------------------------
    # basic helpers
    @property
    def size(self) -> int:
        return self.data.size

    def to_array(self) -> np.ndarray:
        return self.data

    @property
    def references(self) -> OrderedDict:
        """returns OrderedDict of all names and aliases such that keys all point to the current instance (`self`)

        the idea is that this should be placed in a `ChainMap` object to be used by the :class:`axis_collection` class that contains the axis.
        """
        return OrderedDict((name, self) for name in self.names)

    # ------------------------------------------------------------------
    # slicing and indexing
    def __getitem__(self, item):
        if isinstance(item, tuple):
            lo, hi = item
            lo = -np.inf if lo is None else lo
            hi = np.inf if hi is None else hi
            mask = (self.data >= lo) & (self.data <= hi)
            return nddata_axis(self.data[mask])
        if isinstance(item, slice):
            return nddata_axis(self.data[item])
        return self.data[item]

    # ------------------------------------------------------------------
    # arithmetic helpers
    def _reshape_for(self, other: "nddata"):
        idx = other.dimlabels.index(self.names[0])
        shape = [1] * other.data.ndim
        shape[idx] = self.size
        return self.data.reshape(shape)

    def __mul__(self, other):
        if isinstance(other, np.ndarray):
            return self.data * other
        if hasattr(other, "data") and hasattr(other, "dimlabels"):
            arr = self._reshape_for(other)
            result = other.__class__(other.data * arr, other.dimlabels)
            result.axes = axis_collection(result.dimlabels)
            for name in result.dimlabels:
                if name == self.names[0]:
                    result.axes[name] = self
                else:
                    result.axes[name] = other.axes[name]
            return result
        return self.data * other

    __rmul__ = __mul__

    # ------------------------------------------------------------------
    # numpy ufunc support
    def __array__(self, dtype=None):  # pragma: no cover - numpy protocol
        return np.asarray(self.data, dtype=dtype)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != "__call__":
            return NotImplemented
        arrays = [x.data if isinstance(x, nddata_axis) else x for x in inputs]
        result = getattr(ufunc, method)(*arrays, **kwargs)
        return nddata_axis(result)


# delayed import to avoid circular dependency when type checking
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from .core import nddata  # noqa: F401

