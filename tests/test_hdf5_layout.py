import sys
import types
import numpy as np
import h5py

from conftest import load_module

# ensure optional compiled dependency is stubbed
sys.modules.setdefault("_nnls", types.ModuleType("_nnls"))
# numpy.core.rec was removed in newer numpy, but nddata.hdf5_write expects it
rec = types.ModuleType("rec")
try:
    rec.fromarrays = np.core.records.fromarrays
except Exception:
    rec.fromarrays = np.core.records.fromarrays
sys.modules["numpy.core.rec"] = rec
np.core.rec = rec

# load nddata using helper to avoid requiring full dependencies
core = load_module("core", use_real_pint=True)
nddata = core.nddata


def _generate_nddata():
    a = nddata(np.arange(6).reshape(2, 3), ["t", "f"])
    a.labels({"t": np.linspace(0.0, 1.0, 2), "f": np.array([10, 20, 30])})
    a.set_units("t", "s")
    a.set_units("f", "Hz")
    a.set_units("V")
    a.name("test_nd")
    a.other_info.update({"level1": {"level2": 5, "level2list": [1, 2, 3]}})
    return a


def _check_layout(g, a):
    def get_value(obj, key):
        if hasattr(obj, "attrs") and key in obj.attrs:
            return obj.attrs[key]
        return obj[key]

    assert "data" in g
    raw_labels = get_value(g, "dimlabels")
    dimlabels = []
    for lbl in raw_labels:
        if not isinstance(lbl, (bytes, np.bytes_, str)) and hasattr(
            lbl, "__getitem__"
        ):
            lbl = lbl[0]
        if isinstance(lbl, (bytes, np.bytes_)):
            lbl = lbl.decode("utf-8")
        dimlabels.append(lbl)
    assert dimlabels == ["t", "f"]

    axes_group = g["axes"]
    for name, coords, unit in [
        ("t", a.getaxis("t"), "s"),
        ("f", a.getaxis("f"), "Hz"),
    ]:
        assert name in axes_group
        ds = axes_group[name]
        np.testing.assert_allclose(np.array(ds["data"]), coords)
        axis_unit = get_value(ds, "axis_coords_units")
        if isinstance(axis_unit, (bytes, np.bytes_)):
            axis_unit = axis_unit.decode()
        assert axis_unit == unit

    other_info = g["other_info"]
    assert "level1" in other_info
    level1 = other_info["level1"]
    assert get_value(level1, "level2") == 5
    lvl2list = get_value(level1, "level2list")
    if isinstance(lvl2list, dict):
        arr = lvl2list["LISTELEMENTS"]
    else:
        arr = np.array(lvl2list["LISTELEMENTS"])
    np.testing.assert_array_equal(arr, [1, 2, 3])


def test_hdf5_layout(tmp_path):
    a = _generate_nddata()
    a.hdf5_write("sample.h5", directory=str(tmp_path))
    with h5py.File(tmp_path / "sample.h5", "r") as f:
        assert "test_nd" in f
        _check_layout(f["test_nd"], a)


def test_hdf5_layout_from_state():
    a = _generate_nddata()
    g = a.__getstate__()
    _check_layout(g, a)
