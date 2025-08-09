import sys
import types
import numpy as np
import h5py
import pickle

from conftest import load_module

# ensure optional compiled dependency is stubbed
sys.modules.setdefault("_nnls", types.ModuleType("_nnls"))

# load nddata using helper to avoid requiring full dependencies
core = load_module("core", use_real_pint=True)
nddata = core.nddata
nddata_hdf5 = core.nddata_hdf5


class DummyGroup(dict):
    """Minimal dict subclass mimicking an h5py group/dataset.

    Accessing ``attrs['key']`` simply returns ``self['key']``.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.attrs = self
        for k, v in self.items():
            if type(v) is dict:
                self[k] = DummyGroup(v)


def _generate_nddata():
    data = np.arange(6).reshape(2, 3) + 1j * np.arange(1, 7).reshape(2, 3)
    a = nddata(data, ["t", "f"])
    a.labels({"t": np.linspace(0.0, 1.0, 2), "f": np.array([10, 20, 30])})
    a.set_units("t", "s")
    a.set_units("f", "Hz")
    a.set_units("V")
    a.set_error(np.full(a.data.shape, 0.1))
    a.set_error("t", np.array([0.01, 0.02]))
    a.set_error("f", np.array([0.1, 0.2, 0.3]))
    a.name("test_nd")
    a.other_info.update({"level1": {"level2": 5, "level2list": [1, 2, 3]}})
    return a


def _check_layout(g, a):
    assert "data" in g
    raw_labels = g.attrs["dimlabels"]
    if hasattr(raw_labels, "astype"):
        raw_labels = raw_labels.astype("S")
    dimlabels = []
    for lbl in raw_labels:
        if isinstance(lbl, (bytes, np.bytes_)):
            lbl = lbl.decode("utf-8")
        dimlabels.append(lbl)
    assert dimlabels == ["t", "f"]

    axes_group = g["axes"]
    shape = tuple(len(axes_group[name]["data"]) for name in dimlabels)
    data_ds = g["data"]
    if hasattr(data_ds, "keys") and "data.r" in data_ds:
        data_vals = (
            np.array(data_ds["data.r"]) + 1j * np.array(data_ds["data.i"])
        ).reshape(shape)
    else:
        data_vals = np.array(data_ds["data"]).reshape(shape)
    error_vals = np.array(data_ds["error"]).reshape(shape)
    np.testing.assert_allclose(data_vals.real, a.data.real)
    np.testing.assert_allclose(data_vals.imag, a.data.imag)
    np.testing.assert_allclose(error_vals, a.get_error())

    for name in dimlabels:
        coords = a.getaxis(name)
        unit = "s" if name == "t" else "Hz"
        assert name in axes_group
        ds = axes_group[name]
        np.testing.assert_allclose(ds["data"], coords)
        np.testing.assert_allclose(ds["error"], a.get_error(name))
        axis_unit = ds.attrs["axis_coords_units"]
        if isinstance(axis_unit, (bytes, np.bytes_)):
            axis_unit = axis_unit.decode()
        assert axis_unit == unit

    other_info = g["other_info"]
    assert "level1" in other_info
    lvl1 = other_info["level1"]
    assert lvl1.attrs["level2"] == 5
    if "level2list" in lvl1.attrs:
        lvl2list = lvl1.attrs["level2list"]
        if (
            isinstance(lvl2list, dict) and "LISTELEMENTS" in lvl2list
        ):
            lvl2list = lvl2list["LISTELEMENTS"]
        elif getattr(lvl2list, "dtype", None) is not None and "LISTELEMENTS" in getattr(lvl2list.dtype, "names", []) :
            lvl2list = lvl2list["LISTELEMENTS"]
    else:
        lvl2list = lvl1["level2list"].attrs["LISTELEMENTS"]
    np.testing.assert_array_equal(lvl2list, [1, 2, 3])


def _check_loaded(b, a):
    assert list(b.dimlabels) == list(a.dimlabels)
    np.testing.assert_allclose(b.data.real, a.data.real)
    np.testing.assert_allclose(b.data.imag, a.data.imag)
    np.testing.assert_allclose(b.get_error(), a.get_error())
    for name in a.dimlabels:
        np.testing.assert_allclose(b.getaxis(name), a.getaxis(name))
        np.testing.assert_allclose(b.get_error(name), a.get_error(name))
        assert b.get_units(name) == a.get_units(name)
    assert b.get_units() == a.get_units()
    assert b.name() == a.name()
    assert b.other_info == a.other_info


def test_hdf5_layout(tmp_path):
    a = _generate_nddata()
    a.hdf5_write("sample.h5", directory=str(tmp_path))
    with h5py.File(tmp_path / "sample.h5", "r") as f:
        assert "test_nd" in f
        _check_layout(f["test_nd"], a)


def test_state_layout():
    a = _generate_nddata()
    g = DummyGroup(a.__getstate__())
    _check_layout(g, a)


def test_nddata_hdf5_roundtrip(tmp_path):
    a = _generate_nddata()
    a.hdf5_write("sample.h5", directory=str(tmp_path))
    b = nddata_hdf5("sample.h5/test_nd", directory=str(tmp_path))
    _check_loaded(b, a)


def test_nddata_pickle_roundtrip(tmp_path):
    a = _generate_nddata()
    with open(tmp_path / "sample.pkl", "wb") as f:
        pickle.dump(a, f)
    with open(tmp_path / "sample.pkl", "rb") as f:
        b = pickle.load(f)
    _check_loaded(b, a)


def test_state_hdf_dict_roundtrip(tmp_path):
    a = _generate_nddata()
    state = a.__getstate__()
    hdf_mod = load_module("file_saving.hdf_save_dict_to_group")
    hdf_save_dict_to_group = hdf_mod.hdf_save_dict_to_group
    with h5py.File(tmp_path / "state.h5", "w") as f:
        g = f.create_group("test_nd")
        hdf_save_dict_to_group(g, state)
    with h5py.File(tmp_path / "state.h5", "r") as f:
        _check_layout(f["test_nd"], a)
