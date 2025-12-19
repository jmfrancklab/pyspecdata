import numpy as np
from collections import OrderedDict
from copy import deepcopy


class axis_collection(object):
    """Manage a set of :class:`nddata_axis` objects keyed by dimension labels."""

    def __init__(self, dimlabels=None, on_update=None):
        if dimlabels is None:
            dimlabels = []
        self._dimlabels = list(dimlabels)
        self._fixed_labels = len(self._dimlabels) > 0
        self._axes = {}
        self._on_update = on_update

    @property
    def dimlabels(self):
        return self._dimlabels

    @dimlabels.setter
    def dimlabels(self, value):
        self._fixed_labels = True
        self._dimlabels = list(value)

    def _trigger_update(self):
        if self._on_update is not None:
            self._on_update(self)

    def _check_alias_conflict(self, axis, skip_label=None):
        for label, existing in self._axes.items():
            if skip_label is not None and label == skip_label:
                continue
            for name in existing.names:
                if name in axis.names:
                    raise ValueError(
                        "Axis aliases must be unique; %s is already used" % name
                    )

    def _resolve_label(self, axis, provided_label=None):
        preferred = axis.names[0] if len(axis.names) > 0 else None
        if provided_label is None:
            if len(self.dimlabels) == 0:
                return preferred
            for candidate in axis.names:
                if candidate in self.dimlabels:
                    return candidate
            return preferred
        return provided_label

    def __getitem__(self, key):
        return self._axes[key]

    def __setitem__(self, key, axis):
        if not isinstance(axis, nddata_axis):
            raise TypeError("Only nddata_axis instances can be stored")
        if len(axis.names) == 0:
            raise ValueError("Axis must have at least one name")
        self._check_alias_conflict(axis, skip_label=key if key in self else None)
        label = self._resolve_label(axis, provided_label=key)
        if self._fixed_labels and label not in self.dimlabels:
            raise ValueError(
                "Axis names %s do not match any existing dimlabel"
                % repr(axis.names)
            )
        if label not in self.dimlabels:
            self.dimlabels.append(label)
        axis.names = [label] + [x for x in axis.names if x != label]
        self._axes[label] = axis
        self._trigger_update()

    def __contains__(self, key):
        return key in self._axes

    def __len__(self):
        return len(self._axes)

    def __iadd__(self, axis):
        self.append(axis)
        return self

    def append(self, axis):
        if not isinstance(axis, nddata_axis):
            raise TypeError("Only nddata_axis instances can be appended")
        self._check_alias_conflict(axis)
        label = self._resolve_label(axis)
        if self._fixed_labels and label not in self.dimlabels:
            raise ValueError(
                "Axis names %s do not match any existing dimlabel"
                % repr(axis.names)
            )
        if label not in self.dimlabels:
            self.dimlabels.append(label)
        axis.names = [label] + [x for x in axis.names if x != label]
        self._axes[label] = axis
        self._trigger_update()

    def rename(self, oldname, newname):
        if oldname not in self._axes:
            raise ValueError("Cannot rename missing axis %s" % oldname)
        axis = self._axes.pop(oldname)
        axis.names = [newname] + [x for x in axis.names if x != newname]
        if oldname in self._dimlabels:
            idx = self._dimlabels.index(oldname)
            self._dimlabels[idx] = newname
        else:
            self._dimlabels.append(newname)
        self._axes[newname] = axis
        self._trigger_update()
        return self

    @property
    def units(self):
        return {lbl: self._axes[lbl].units for lbl in self._axes}


class _ax_class_maker(object):
    def __getitem__(self, args):
        if not isinstance(args, tuple) or len(args) == 0:
            raise ValueError("You must provide an axis name before the slice")
        name = args[0]
        if not isinstance(name, str):
            raise ValueError("Axis name must come before the slice")
        if len(args) == 1:
            raise ValueError("You must provide coordinates for the axis")
        coordinates = args[1]
        return nddata_axis(name, coordinates)


ax_ = _ax_class_maker()


class nddata_axis(object):
    """Axis object that stores coordinates and metadata."""

    __array_priority__ = 1000

    def __init__(self, name=None, values=None):
        self.domains = OrderedDict()
        self.multiplier = None
        self.transf_func = None
        self.units = None
        self.names = []
        if name is not None:
            self.names = [name]
        if values is None:
            self.data = np.array([])
            return
        self.data = self._generate_data(values)

    def _generate_data(self, values):
        if isinstance(values, slice):
            start = 0.0 if values.start is None else values.start
            stop = values.stop
            step = values.step
            if isinstance(step, complex):
                count = int(abs(step.imag))
                return np.linspace(start, stop, count)
            if step is None:
                return np.arange(start, stop)
            return np.arange(start, stop, step)
        return np.array(values)

    def copy(self):
        duplicate = nddata_axis(
            self.names[0] if len(self.names) > 0 else None, self.data.copy()
        )
        duplicate.names = list(self.names)
        duplicate.units = self.units
        duplicate.transf_func = self.transf_func
        duplicate.multiplier = self.multiplier
        duplicate.domains = deepcopy(self.domains)
        return duplicate

    @property
    def size(self):
        return len(self.data)

    @property
    def start(self):
        if self.size == 0:
            return None
        return self.data[0]

    @property
    def dx(self):
        if self.size < 2:
            return None
        return self.data[1] - self.data[0]

    def to_array(self):
        return np.array(self.data)

    def __array__(self, dtype=None):
        if dtype is None:
            return np.array(self.data)
        return np.array(self.data, dtype=dtype)

    def __len__(self):
        return len(self.data)

    def _is_nddata(self, obj):
        classname = obj.__class__.__name__
        return classname == "nddata"

    def _operate_with_nddata(self, dataset, operation):
        if len(self.names) == 0:
            raise ValueError("Axis needs a name to operate on nddata")
        label = self.names[0]
        if label not in dataset.dimlabels:
            raise ValueError(
                "Axis name %s is not present on the dataset" % label
            )
        if dataset.data.shape[dataset.axn(label)] != self.size:
            raise ValueError(
                "Axis length does not match the dataset dimension for %s" % label
            )
        newdata = dataset.copy()
        scale_shape = [1] * newdata.data.ndim
        scale_shape[newdata.axn(label)] = self.size
        scale = self.data.reshape(scale_shape)
        newdata.data = operation(newdata.data, scale)
        newdata.setaxis(label, self.copy())
        return newdata

    def __getitem__(self, key):
        if isinstance(key, tuple):
            if len(key) != 2:
                raise ValueError("Slices must be (min,max) when using a tuple")
            lower, upper = key
            mask = np.ones(self.data.shape, dtype=bool)
            if lower is not None:
                mask = np.logical_and(mask, self.data >= lower)
            if upper is not None:
                mask = np.logical_and(mask, self.data <= upper)
            return nddata_axis(
                self.names[0] if len(self.names) > 0 else None, self.data[mask]
            )
        if isinstance(key, slice):
            return nddata_axis(
                self.names[0] if len(self.names) > 0 else None, self.data[key]
            )
        return self.data[key]

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != "__call__":
            return NotImplemented
        arrays = []
        for item in inputs:
            if item is self:
                arrays.append(self.data)
            elif isinstance(item, nddata_axis):
                arrays.append(item.data)
            else:
                arrays.append(item)
        result = ufunc(*arrays, **kwargs)
        if isinstance(result, tuple):
            return tuple(
                nddata_axis(self.names[0], r)
                if isinstance(r, np.ndarray) and r.shape == self.data.shape
                else r
                for r in result
            )
        if isinstance(result, np.ndarray) and result.shape == self.data.shape:
            transformed = nddata_axis(
                self.names[0] if len(self.names) > 0 else None, result
            )
            transformed.units = self.units
            transformed.transf_func = ufunc
            return transformed
        return result

    def __add__(self, other):
        if isinstance(other, nddata_axis):
            if len(self.names) > 0 and len(other.names) > 0:
                if self.names[0] == other.names[0]:
                    if self.size != other.size:
                        raise ValueError("Axis lengths must match for addition")
                    result = nddata_axis(self.names[0], self.data + other.data)
                    result.units = self.units
                    return result
            return self.data[:, None] + other.data[None, :]
        return self.data + other

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, nddata_axis):
            if len(self.names) > 0 and len(other.names) > 0:
                if self.names[0] == other.names[0]:
                    if self.size != other.size:
                        raise ValueError("Axis lengths must match for subtraction")
                    result = nddata_axis(self.names[0], self.data - other.data)
                    result.units = self.units
                    return result
            return self.data[:, None] - other.data[None, :]
        return self.data - other

    def __rsub__(self, other):
        if isinstance(other, nddata_axis):
            return other.__sub__(self)
        return other - self.data

    def __mul__(self, other):
        if isinstance(other, nddata_axis):
            if len(self.names) > 0 and len(other.names) > 0:
                if self.names[0] == other.names[0]:
                    if self.size != other.size:
                        raise ValueError("Axis lengths must match for multiplication")
                    result = nddata_axis(self.names[0], self.data * other.data)
                    result.units = self.units
                    return result
            return self.data[:, None] * other.data[None, :]
        if self._is_nddata(other):
            return self._operate_with_nddata(other, np.multiply)
        if isinstance(other, np.ndarray):
            return self.data * other
        if np.isscalar(other):
            result = self.copy()
            result.data = self.data * other
            return result
        return NotImplemented

    def __rmul__(self, other):
        if self._is_nddata(other):
            return self._operate_with_nddata(other, np.multiply)
        if isinstance(other, np.ndarray):
            return other * self.data
        if np.isscalar(other):
            result = self.copy()
            result.data = other * self.data
            return result
        return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, nddata_axis):
            if len(self.names) > 0 and len(other.names) > 0:
                if self.names[0] == other.names[0]:
                    if self.size != other.size:
                        raise ValueError("Axis lengths must match for division")
                    result = nddata_axis(self.names[0], self.data / other.data)
                    result.units = self.units
                    return result
            return self.data[:, None] / other.data[None, :]
        if np.isscalar(other):
            result = self.copy()
            result.data = self.data / other
            return result
        if isinstance(other, np.ndarray):
            return self.data / other
        return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, np.ndarray):
            return other / self.data
        if np.isscalar(other):
            return other / self.data
        return NotImplemented

    def __pow__(self, power):
        result = self.copy()
        result.data = self.data ** power
        return result
