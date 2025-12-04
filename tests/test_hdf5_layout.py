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

    def create_group(self, name):
        """Create a nested :class:`DummyGroup` and attach it under *name*."""

        subgroup = DummyGroup()
        self[name] = subgroup
        return subgroup

    def create_dataset(self, name, data=None, dtype=None):
        """Create a :class:`DummyDataset` and attach it under *name*."""

        dataset = DummyDataset(data)
        self[name] = dataset
        return dataset


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


def _check_state(state, a, haserr=True):
    expected_keys = {"data", "dimlabels", "axes", "other_info"}
    assert expected_keys.issubset(set(state.keys()))

    dimlabels = list(state["dimlabels"])
    assert dimlabels == list(a.dimlabels)

    np.testing.assert_allclose(state["data"]["data"], a.data)
    if haserr:
        np.testing.assert_allclose(state["data"]["error"], a.get_error())
    if a.get_units() is None:
        assert state["data"]["data_units"] is None
    else:
        assert state["data"]["data_units"] == a.get_units()

    for lbl in dimlabels:
        assert "NUMPY_DATA" in state["axes"][lbl]
        axis_array = state["axes"][lbl]["NUMPY_DATA"]
        if axis_array.dtype.names is None:
            np.testing.assert_allclose(axis_array, a.getaxis(lbl))
        else:
            np.testing.assert_allclose(axis_array["data"], a.getaxis(lbl))
            if haserr:
                np.testing.assert_allclose(
                    axis_array["error"], a.get_error(lbl)
                )
        if a.get_units(lbl) is None:
            assert "axis_coords_units" not in state["axes"][lbl]
        else:
            assert state["axes"][lbl]["axis_coords_units"] == a.get_units(lbl)

    assert state["other_info"] == a.other_info
    assert isinstance(state["other_info"]["level1"]["level2list"], list)


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
    state = a.__getstate__()
    _check_state(state, a)


def test_state_layout_noerr():
    a = _generate_nddata_noerr()
    state = a.__getstate__()
    _check_state(state, a, haserr=False)


def test_state_axis_error_structured_merge():
    a = _generate_nddata_noerr()
    structured_axis = np.core.records.fromarrays(
        [a.getaxis("t")], names="data"
    )
    a.axis_coords[0] = structured_axis
    a.set_error("t", np.array([0.01, 0.02]))
    state = a.__getstate__()
    assert state["axes"]["t"]["NUMPY_DATA"].dtype.names == ("data", "error")
    np.testing.assert_allclose(
        state["axes"]["t"]["NUMPY_DATA"]["data"], a.getaxis("t")["data"]
    )
    np.testing.assert_allclose(
        state["axes"]["t"]["NUMPY_DATA"]["error"], a.get_error("t")
    )


def test_state_axis_error_plain_to_structured():
    a = _generate_nddata()
    state = a.__getstate__()
    for lbl in a.dimlabels:
        assert state["axes"][lbl]["NUMPY_DATA"].dtype.names == (
            "data",
            "error",
        )
        np.testing.assert_allclose(
            state["axes"][lbl]["NUMPY_DATA"]["data"], a.getaxis(lbl)
        )
        np.testing.assert_allclose(
            state["axes"][lbl]["NUMPY_DATA"]["error"], a.get_error(lbl)
        )


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
        "mixed_list": [
            1,
            "two",
            {"inner_dict": {"label": "val"}},
            (3.5, "units"),
        ],
        "tuple_wrapping": (["zero", 0], {"deep_list": [("keep", 1), "str"]}),
    }
    a.set_prop("complex_prop", complex_prop)
    a.hdf5_write("sample.h5", directory=str(tmp_path))
    with h5py.File(tmp_path / "sample.h5", "r") as f:
        assert "mixed_list" in f["test_nd"]["other_info"]["complex_prop"]
        assert f["test_nd"]["other_info"]["complex_prop"]["mixed_list"].attrs[
            "LIST_NODE"
        ]
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
    encoded = encode_list("seq", seq, False)
    hdf_save_dict_to_group(root, encoded, use_pytables_hack=False)
    reconstructed = decode_list(root["seq"])
    assert reconstructed == seq


def test_encode_decode_list_invertible_with_hack():
    seq = [1, 2, 3]
    root = DummyGroup()
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


def test_legacy_axis_loading(tmp_path):
    data = np.arange(6).reshape(2, 3)
    with h5py.File(tmp_path / "legacy.h5", "w") as f:
        g = f.create_group("test_nd")
        g.attrs["dimlabels"] = np.array([(b"t",), (b"f",)])
        axes_group = g.create_group("axes")
        legacy_t = np.core.records.fromarrays(
            [np.linspace(0.0, 1.0, 2)], names="data"
        )
        axes_group.create_dataset("t", data=legacy_t)
        axes_group["t"].attrs["axis_coords_units"] = b"s"
        axes_group.create_dataset("f", data=np.array([10, 20, 30]))
        axes_group["f"].attrs["axis_coords_units"] = b"Hz"
        data_group = g.create_group("data")
        data_group.create_dataset("data", data=data)
        data_group.attrs["data_units"] = b"V"
    loaded = nddata_hdf5("legacy.h5/test_nd", directory=str(tmp_path))
    np.testing.assert_allclose(loaded.getaxis("t"), np.linspace(0.0, 1.0, 2))
    np.testing.assert_allclose(loaded.getaxis("f"), np.array([10, 20, 30]))
    assert loaded.get_units("t") == "s"
    assert loaded.get_units("f") == "Hz"
    assert loaded.get_units() == "V"
    np.testing.assert_allclose(loaded.data.real, data)
    os.remove(tmp_path / "legacy.h5")


def test_dataset_node_for_data_loading(tmp_path):
    data = np.arange(6).reshape(2, 3)
    with h5py.File(tmp_path / "dataset_data.h5", "w") as f:
        g = f.create_group("test_nd")
        g.attrs["dimlabels"] = np.array([(b"t",), (b"f",)])
        axes_group = g.create_group("axes")
        axes_group.create_dataset("t", data=np.linspace(0.0, 1.0, 2))
        axes_group["t"].attrs["axis_coords_units"] = b"s"
        axes_group.create_dataset("f", data=np.array([10, 20, 30]))
        axes_group["f"].attrs["axis_coords_units"] = b"Hz"
        g.create_dataset("data", data=data)
        g["data"].attrs["data_units"] = b"V"
    loaded = nddata_hdf5("dataset_data.h5/test_nd", directory=str(tmp_path))
    assert list(loaded.dimlabels) == ["t", "f"]
    np.testing.assert_allclose(loaded.data.real, data)
    np.testing.assert_allclose(loaded.getaxis("t"), np.linspace(0.0, 1.0, 2))
    np.testing.assert_allclose(loaded.getaxis("f"), np.array([10, 20, 30]))
    assert loaded.get_units("t") == "s"
    assert loaded.get_units("f") == "Hz"
    assert loaded.get_units() == "V"
    os.remove(tmp_path / "dataset_data.h5")


def test_structured_dataset_data_unwrap(tmp_path):
    data = np.arange(6).reshape(2, 3) + 1j * np.arange(1, 7).reshape(2, 3)
    structured_data = np.zeros(data.shape, dtype=[("data", "<c16")])
    structured_data["data"] = data

    with h5py.File(tmp_path / "structured_data.h5", "w") as f:
        g = f.create_group("test_nd")
        g.attrs["dimlabels"] = np.array([(b"t",), (b"f",)])
        axes_group = g.create_group("axes")
        axes_group.create_dataset("t", data=np.linspace(0.0, 1.0, 2))
        axes_group["t"].attrs["axis_coords_units"] = b"s"
        axes_group.create_dataset("f", data=np.array([10, 20, 30]))
        axes_group["f"].attrs["axis_coords_units"] = b"Hz"
        data_group = g.create_group("data")
        data_group.create_dataset("data", data=structured_data)
        data_group.attrs["data_units"] = b"V"

    loaded = nddata_hdf5("structured_data.h5/test_nd", directory=str(tmp_path))
    assert list(loaded.dimlabels) == ["t", "f"]
    assert loaded.data.dtype.names is None
    np.testing.assert_allclose(loaded.data, data)
    np.testing.assert_allclose(loaded.getaxis("t"), np.linspace(0.0, 1.0, 2))
    np.testing.assert_allclose(loaded.getaxis("f"), np.array([10, 20, 30]))
    assert loaded.get_units("t") == "s"
    assert loaded.get_units("f") == "Hz"
    assert loaded.get_units() == "V"
    os.remove(tmp_path / "structured_data.h5")


def test_pytables_metadata_attributes_ignored(tmp_path):
    data = np.arange(4.0).reshape(2, 2)
    with h5py.File(tmp_path / "pytables_meta.h5", "w") as f:
        g = f.create_group("test_nd")
        g.attrs["dimlabels"] = np.array([(b"t",), (b"f",)])
        axes_group = g.create_group("axes")
        t_ds = axes_group.create_dataset("t", data=np.linspace(0.0, 1.0, 2))
        t_ds.attrs["CLASS"] = b"ARRAY"
        t_ds.attrs["TITLE"] = b"time axis"
        t_ds.attrs["VERSION"] = b"1.0"
        f_ds = axes_group.create_dataset("f", data=np.array([10, 20]))
        f_ds.attrs["CLASS"] = b"ARRAY"
        data_group = g.create_group("data")
        main_ds = data_group.create_dataset("data", data=data)
        main_ds.attrs["CLASS"] = b"ARRAY"
        main_ds.attrs["TITLE"] = b"data"
        main_ds.attrs["VERSION"] = b"2.0"
        data_group.attrs["data_units"] = b"V"
        other_info_group = g.create_group("other_info")
        cp_ds = other_info_group.create_dataset(
            "coherence_pathway", data=np.array([1, -1])
        )
        cp_ds.attrs["CLASS"] = b"ARRAY"
        cp_ds.attrs["TITLE"] = b"pathway"
        cp_ds.attrs["VERSION"] = b"1.0"

    with h5py.File(tmp_path / "pytables_meta.h5", "r") as f:
        loaded_state = hdf_load_dict_from_group(f["test_nd"])
    assert "CLASS" not in loaded_state["axes"]["t"]
    assert "TITLE" not in loaded_state["axes"]["t"]
    assert "VERSION" not in loaded_state["axes"]["t"]
    assert "CLASS" not in loaded_state["axes"]["f"]
    assert "CLASS" not in loaded_state["data"]["data"]
    assert "TITLE" not in loaded_state["data"]["data"]
    assert "VERSION" not in loaded_state["data"]["data"]
    assert set(loaded_state["other_info"]["coherence_pathway"].keys()) == {
        "NUMPY_DATA"
    }

    loaded = nddata_hdf5("pytables_meta.h5/test_nd", directory=str(tmp_path))
    assert "coherence_pathway" in loaded.other_info
    assert "CLASS" not in loaded.other_info["coherence_pathway"]
    assert loaded.data.shape == (2, 2)
    np.testing.assert_allclose(loaded.getaxis("t"), np.linspace(0.0, 1.0, 2))
    np.testing.assert_allclose(loaded.getaxis("f"), np.array([10, 20]))
    os.remove(tmp_path / "pytables_meta.h5")


def test_attributes_of_main_tree_roundtrip(tmp_path):
    a = _generate_nddata_noerr()
    state = a.__getstate__()
    with h5py.File(tmp_path / "attrs.h5", "w") as f:
        g = f.create_group("test_nd")
        hdf_save_dict_to_group(g, state, use_pytables_hack=False)
    with h5py.File(tmp_path / "attrs.h5", "r") as f:
        axis_t_units = f["test_nd"]["axes"]["t"].attrs["axis_coords_units"]
        axis_f_units = f["test_nd"]["axes"]["f"].attrs["axis_coords_units"]
        data_units_attr = f["test_nd"]["data"].attrs["data_units"]
        if isinstance(axis_t_units, (bytes, np.bytes_)):
            assert axis_t_units == b"s"
        else:
            assert axis_t_units == "s"
        if isinstance(axis_f_units, (bytes, np.bytes_)):
            assert axis_f_units == b"Hz"
        else:
            assert axis_f_units == "Hz"
        if isinstance(data_units_attr, (bytes, np.bytes_)):
            assert data_units_attr == b"V"
        else:
            assert data_units_attr == "V"
    with h5py.File(tmp_path / "attrs.h5", "r") as f:
        loaded_state = hdf_load_dict_from_group(f["test_nd"])
    assert loaded_state["axes"]["t"]["axis_coords_units"] == "s"
    assert loaded_state["axes"]["f"]["axis_coords_units"] == "Hz"
    assert loaded_state["data"]["data_units"] == "V"
    os.remove(tmp_path / "attrs.h5")


def test_pytables_hack_dimlabels_loading(tmp_path):
    # construct a dimlabels attribute that matches the PyTables structured
    # array style so we can confirm it decodes to plain strings
    structured_labels = np.zeros(1, dtype=[("LISTELEMENTS", "S5")])
    structured_labels["LISTELEMENTS"][0] = np.bytes_("time")
    data = np.arange(3.0)
    with h5py.File(tmp_path / "pytables_dimlabels.h5", "w") as f:
        g = f.create_group("test_nd")
        g.attrs["dimlabels"] = structured_labels
        axes_group = g.create_group("axes")
        axes_group.create_dataset("time", data=np.array([0.0, 1.0, 2.0]))
        data_group = g.create_group("data")
        data_group.create_dataset("data", data=data)
    loaded = nddata_hdf5(
        "pytables_dimlabels.h5/test_nd", directory=str(tmp_path)
    )
    assert list(loaded.dimlabels) == ["time"]
    np.testing.assert_allclose(loaded.data.real, data)
    np.testing.assert_allclose(
        loaded.getaxis("time"), np.array([0.0, 1.0, 2.0])
    )
