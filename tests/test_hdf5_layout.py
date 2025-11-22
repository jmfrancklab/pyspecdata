import sys
import types
import numpy as np
import pytest
import h5py
import pickle
import os

# tests/conftest.py
from conftest import load_module

# ensure optional compiled dependency is stubbed
sys.modules.setdefault("_nnls", types.ModuleType("_nnls"))

# load nddata using helper to avoid requiring full dependencies
core = load_module("core", use_real_pint=True, use_real_h5py=True)
nddata = core.nddata
nddata_hdf5 = core.nddata_hdf5
hmod = load_module(
    "file_saving.hdf_save_dict_to_group",
    use_real_pint=True,
    use_real_h5py=True,
)
hdf_save_dict_to_group = hmod.hdf_save_dict_to_group
decode_list = hmod.decode_list
encode_list = hmod.encode_list
hdf_load_dict_from_group = hmod.hdf_load_dict_from_group


class DummyDataset:
    def __init__(self, data):
        self.data = data
        self.attrs = {}

    def __getitem__(self, key):
        if key == ():
            return self.data
        return self.data[key]

    def __len__(self):
        return len(self.data)


class DummyGroup(dict):
    """Minimal dict subclass mimicking an h5py group/dataset."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.attrs = {}


def _add_group_creators(target):
    """Attach minimal creation helpers to a :class:`DummyGroup` instance.

    These helpers mimic the ``h5py`` API closely enough for ``encode_list`` to
    build list nodes without requiring a full HDF5 file.
    """

    def create_group(name):
        subgroup = DummyGroup()
        _add_group_creators(subgroup)
        target[name] = subgroup
        return subgroup

    def create_dataset(name, data=None, dtype=None):
        dataset = DummyDataset(data)
        target[name] = dataset
        return dataset

    target.create_group = create_group
    target.create_dataset = create_dataset


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


def _generate_nddata_noerr():
    data = np.arange(6).reshape(2, 3) + 1j * np.arange(1, 7).reshape(2, 3)
    a = nddata(data, ["t", "f"])
    a.labels({"t": np.linspace(0.0, 1.0, 2), "f": np.array([10, 20, 30])})
    a.set_units("t", "s")
    a.set_units("f", "Hz")
    a.set_units("V")
    a.name("test_nd")
    a.other_info.update({"level1": {"level2": 5, "level2list": [1, 2, 3]}})
    return a


def _check_layout(g, a, haserr=True):
    assert "data" in g
    print("g:", dir(g))
    print("g.attrs:", g.attrs.keys())
    print("g.attrs['dimlabels']:", g.attrs["dimlabels"])
    dimlabels = [j[0].decode("utf-8") for j in g.attrs["dimlabels"]]
    assert dimlabels[0] == "t"
    assert dimlabels[1] == "f"
    if "data_units" not in g["data"].attrs:
        raise KeyError(
            f"dir(g['data']):{dir(g['data'])}\n"
            f"g.attrs.keys():{g['data'].attrs.keys()}"
        )
    if a.get_units() is not None:
        assert g["data"].attrs["data_units"].decode("utf-8") == a.get_units()
    shape = tuple(len(g["axes"][name]["data"]) for name in dimlabels)
    data_vals = np.array(g["data"]["data"]).reshape(shape)
    np.testing.assert_allclose(data_vals.real, a.data.real)
    np.testing.assert_allclose(data_vals.imag, a.data.imag)
    if haserr:
        error_vals = np.array(g["data"]["error"]).reshape(shape)
        np.testing.assert_allclose(error_vals, a.get_error())
    for name in dimlabels:
        np.testing.assert_allclose(g["axes"][name]["data"], a.getaxis(name))
        if haserr:
            np.testing.assert_allclose(
                g["axes"][name]["error"], a.get_error(name)
            )
        assert g["axes"][name].attrs["axis_coords_units"].decode(
            "utf-8"
        ) == a.get_units(name)
    assert "level1" in g["other_info"]
    lvl1 = g["other_info"]["level1"]
    assert lvl1.attrs["level2"] == 5
    np.testing.assert_array_equal(
        lvl1.attrs["level2list"]["LISTELEMENTS"], [1, 2, 3]
    )


def _check_loaded(b, a):
    assert list(b.dimlabels) == list(a.dimlabels)
    np.testing.assert_allclose(b.data.real, a.data.real)
    np.testing.assert_allclose(b.data.imag, a.data.imag)
    if a.get_error() is not None:
        np.testing.assert_allclose(b.get_error(), a.get_error())
    assert a.data.shape == b.data.shape
    for name in a.dimlabels:
        np.testing.assert_allclose(b.getaxis(name), a.getaxis(name))
        if a.get_error(name) is not None:
            np.testing.assert_allclose(b.get_error(name), a.get_error(name))
        assert b.get_units(name) == a.get_units(name)
    assert b.name() == a.name()
    assert b.get_units() == a.get_units()
    assert b.other_info == a.other_info


def test_hdf5_layout(tmp_path):
    a = _generate_nddata()
    a._pytables_hack = True
    a.hdf5_write("sample.h5", directory=str(tmp_path))
    with h5py.File(tmp_path / "sample.h5", "r") as f:
        assert "test_nd" in f
        _check_layout(f["test_nd"], a)
    os.remove(os.path.join(str(tmp_path), "sample.h5"))


def test_hdf5_layout_noerr(tmp_path):
    a = _generate_nddata_noerr()
    a._pytables_hack = True
    a.hdf5_write("sample.h5", directory=str(tmp_path))
    with h5py.File(tmp_path / "sample.h5", "r") as f:
        assert "test_nd" in f
        _check_layout(f["test_nd"], a, haserr=False)
    os.remove(os.path.join(str(tmp_path), "sample.h5"))


def test_state_layout():
    a = _generate_nddata()
    a._pytables_hack = True
    g = DummyGroup()
    _add_group_creators(g)
    hdf_save_dict_to_group(g, a.__getstate__(), use_pytables_hack=True)
    _check_layout(g, a)


def test_state_layout_noerr():
    a = _generate_nddata_noerr()
    a._pytables_hack = True
    g = DummyGroup()
    _add_group_creators(g)
    hdf_save_dict_to_group(g, a.__getstate__(), use_pytables_hack=True)
    _check_layout(g, a, haserr=False)


@pytest.mark.parametrize("use_pytables_hack", [False, True])
def test_nddata_hdf5_roundtrip(tmp_path, use_pytables_hack):
    a = _generate_nddata()
    if use_pytables_hack:
        a._pytables_hack = True
    a.hdf5_write("sample.h5", directory=str(tmp_path))
    b = nddata_hdf5("sample.h5/test_nd", directory=str(tmp_path))
    _check_loaded(b, a)
    os.remove(os.path.join(str(tmp_path), "sample.h5"))


@pytest.mark.parametrize("use_pytables_hack", [False, True])
def test_nddata_hdf5_roundtrip_noerr(tmp_path, use_pytables_hack):
    a = _generate_nddata_noerr()
    if use_pytables_hack:
        a._pytables_hack = True
    a.hdf5_write("sample.h5", directory=str(tmp_path))
    b = nddata_hdf5("sample.h5/test_nd", directory=str(tmp_path))
    _check_loaded(b, a)
    os.remove(os.path.join(str(tmp_path), "sample.h5"))


@pytest.mark.parametrize("use_pytables_hack", [False, True])
def test_nddata_hdf5_roundtrip_ikkf(tmp_path, use_pytables_hack):
    a = _generate_nddata()
    if use_pytables_hack:
        a._pytables_hack = True
    a.set_prop("IKKF", ["REAL", "REAL", "REAL"])
    a.hdf5_write("sample.h5", directory=str(tmp_path))
    b = nddata_hdf5("sample.h5/test_nd", directory=str(tmp_path))
    _check_loaded(b, a)
    os.remove(os.path.join(str(tmp_path), "sample.h5"))


def test_nddata_hdf5_prop_tuple_roundtrip(tmp_path):
    a = _generate_nddata_noerr()
    # ensure tuple properties retain numeric types through HDF5 round trips
    a.set_prop("testprop", (3.0, "T"))
    a.set_prop("strint_mixture", ("05", 7))
    a.hdf5_write("sample.h5", directory=str(tmp_path))
    b = nddata_hdf5("sample.h5/test_nd", directory=str(tmp_path))
    prop_val = b.get_prop("testprop")
    assert isinstance(prop_val[0], float)
    assert prop_val[0] == 3.0
    assert prop_val[1] == "T"
    mixed_val = b.get_prop("strint_mixture")
    assert mixed_val[0] == "05"
    assert isinstance(mixed_val[1], int)
    assert mixed_val[1] == 7
    os.remove(os.path.join(str(tmp_path), "sample.h5"))


def test_nddata_hdf5_prop_complex_structures(tmp_path):
    a = _generate_nddata_noerr()
    complex_prop = {
        "mixed_list": [1, "two", {"inner_dict": {"label": "val"}}, (3.5, "units")],
        "tuple_wrapping": (["zero", 0], {"deep_list": [("keep", 1), "str"]}),
    }
    a.set_prop("complex_prop", complex_prop)
    a.hdf5_write("sample.h5", directory=str(tmp_path))
    with h5py.File(tmp_path / "sample.h5", "r") as f:
        assert "mixed_list" in f["test_nd"]["other_info"]["complex_prop"]
        assert (
            f["test_nd"]["other_info"]["complex_prop"]["mixed_list"].attrs[
                "LIST_NODE"
            ]
        )
        decoded_list = decode_list(
            f["test_nd"]["other_info"]["complex_prop"]["mixed_list"]
        )
        assert isinstance(decoded_list[3], tuple)
        assert decoded_list[3][0] == 3.5
        assert decoded_list[3][1] == "units"
    b = nddata_hdf5("sample.h5/test_nd", directory=str(tmp_path))
    assert b.get_prop("complex_prop") == complex_prop
    os.remove(os.path.join(str(tmp_path), "sample.h5"))


def test_encode_decode_list_invertible_no_hack():
    seq = [1, "two", (3.5, "u"), {"inner": 4}]
    root = DummyGroup()
    _add_group_creators(root)
    encoded = encode_list("seq", seq, False)
    hdf_save_dict_to_group(root, encoded, use_pytables_hack=False)
    reconstructed = decode_list(root["seq"])
    assert reconstructed == seq


def test_encode_decode_list_invertible_with_hack():
    seq = [1, 2, 3]
    root = DummyGroup()
    _add_group_creators(root)
    encoded = encode_list("seq", seq, True)
    hdf_save_dict_to_group(root, encoded, use_pytables_hack=True)
    reconstructed = decode_list(root.attrs["seq"])
    assert reconstructed == seq


def test_hdf_dict_preserves_lists_no_hack(tmp_path):
    data = {"mixed_list": [1, "two", (3.0, "u")], "scalar": 5}
    with h5py.File(tmp_path / "dict_no_hack.h5", "w") as f:
        grp = f.create_group("root")
        hdf_save_dict_to_group(grp, data, use_pytables_hack=False)
    with h5py.File(tmp_path / "dict_no_hack.h5", "r") as f:
        loaded = hdf_load_dict_from_group(f["root"])
    assert loaded["mixed_list"] == data["mixed_list"]
    assert loaded["scalar"] == data["scalar"]


def test_hdf_dict_preserves_lists_pytables_hack(tmp_path):
    data = {"uniform_list": [10, 20, 30], "scalar": 8}
    with h5py.File(tmp_path / "dict_hack.h5", "w") as f:
        grp = f.create_group("root")
        hdf_save_dict_to_group(grp, data, use_pytables_hack=True)
    with h5py.File(tmp_path / "dict_hack.h5", "r") as f:
        loaded = hdf_load_dict_from_group(f["root"])
    assert loaded["uniform_list"] == data["uniform_list"]
    assert loaded["scalar"] == data["scalar"]


def test_nddata_pickle_roundtrip(tmp_path):
    a = _generate_nddata()
    with open(tmp_path / "sample.pkl", "wb") as f:
        pickle.dump(a, f)
    with open(tmp_path / "sample.pkl", "rb") as f:
        b = pickle.load(f)
    _check_loaded(b, a)
    os.remove(tmp_path / "sample.pkl")


def test_state_hdf_dict_roundtrip(tmp_path):
    a = _generate_nddata()
    a._pytables_hack = True
    state = a.__getstate__()
    with h5py.File(tmp_path / "state.h5", "w") as f:
        g = f.create_group("test_nd")
        hdf_save_dict_to_group(g, state, use_pytables_hack=True)
    with h5py.File(tmp_path / "state.h5", "r") as f:
        _check_layout(f["test_nd"], a)
    os.remove(tmp_path / "state.h5")


def test_state_hdf_dict_roundtrip_noerr(tmp_path):
    a = _generate_nddata_noerr()
    a._pytables_hack = True
    state = a.__getstate__()
    with h5py.File(tmp_path / "state.h5", "w") as f:
        g = f.create_group("test_nd")
        hdf_save_dict_to_group(g, state, use_pytables_hack=True)
    with h5py.File(tmp_path / "state.h5", "r") as f:
        _check_layout(f["test_nd"], a, haserr=False)
    os.remove(tmp_path / "state.h5")
