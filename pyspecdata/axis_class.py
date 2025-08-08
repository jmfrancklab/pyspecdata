"""Minimal axis utilities used by the test-suite.

This module provides a very small subset of the functionality described in
the documentation.  The implementation only aims to support the behaviour
exercised in the accompanying unit tests and is **not** a drop-in
replacement for the full featured version.
"""

from __future__ import annotations

from collections import ChainMap
from typing import Iterable, List, MutableMapping

import numpy as np


class axis_collection:
    """Container that stores :class:`nddata_axis` instances keyed by name."""

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
    """Representation of coordinates along a single dimension."""

    __array_priority__ = 1000

    def __init__(self, arr: Iterable[float]):
        self.data = np.array(arr, dtype=float)
        self.names: List[str] = []

    # ------------------------------------------------------------------
    # basic helpers
    @property
    def size(self) -> int:
        return self.data.size

    def to_array(self) -> np.ndarray:
        return self.data

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

