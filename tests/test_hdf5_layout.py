import sys
import types

import numpy as np
import h5py
import tables

from conftest import load_module

# ensure optional compiled dependency is stubbed
sys.modules.setdefault('_nnls', types.ModuleType('_nnls'))
# numpy.core.rec was removed in newer numpy, but nddata.hdf5_write expects it
rec = types.ModuleType('rec')
try:
    rec.fromarrays = np.core.records.fromarrays
except Exception:
    rec.fromarrays = np.core.records.fromarrays
sys.modules['numpy.core.rec'] = rec
np.core.rec = rec

# load nddata using helper to avoid requiring full dependencies
core = load_module('core', use_real_pint=True)
nddata = core.nddata


class DummyGroup(dict):
    """Minimal dict subclass mimicking an h5py group/dataset.

    Accessing ``attrs['key']`` simply returns ``self['key']``.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.attrs = self


def dummy_class_from_step1(obj):
    """Recursively convert dictionaries to :class:`DummyGroup`."""

    if isinstance(obj, dict):
        return DummyGroup({k: dummy_class_from_step1(v) for k, v in obj.items()})
    return obj


def _generate_nddata():
    a = nddata(np.arange(6).reshape(2, 3), ["t", "f"])
    a.labels({"t": np.linspace(0.0, 1.0, 2), "f": np.array([10, 20, 30])})
    a.set_units("t", "s")
    a.set_units("f", "Hz")
    a.set_units("V")
    a.name("test_nd")
    a.meta = {"level1": {"level2": 5, "level2list": [1, 2, 3]}}
    return a


def _check_layout(g, a):
    assert "data" in g
    dimlabels = [x[0].decode("utf-8") for x in g.attrs["dimlabels"]]
    assert dimlabels == ["t", "f"]

    axes_group = g["axes"]
    for name, coords, unit in [
        ("t", a.getaxis("t"), "s"),
        ("f", a.getaxis("f"), "Hz"),
    ]:
        assert name in axes_group
        ds = axes_group[name]
        np.testing.assert_allclose(ds["data"], coords)
        axis_unit = ds.attrs["axis_coords_units"]
        if isinstance(axis_unit, (bytes, np.bytes_)):
            axis_unit = axis_unit.decode()
        assert axis_unit == unit

    meta = g["meta"]
    assert "level1" in meta
    assert meta["level1"].attrs["level2"] == 5
    np.testing.assert_array_equal(
        meta["level1"].attrs["level2list"]["LISTELEMENTS"], [1, 2, 3]
    )


def test_hdf5_layout(tmp_path):
    a = _generate_nddata()
    a.hdf5_write("sample.h5", directory=str(tmp_path))
    with h5py.File(tmp_path / "sample.h5", "r") as f:
        assert "test_nd" in f
        _check_layout(f["test_nd"], a)


def test_hdf5_layout_from_state():
    a = _generate_nddata()
    g = dummy_class_from_step1(a.__getstate__())
    _check_layout(g, a)
