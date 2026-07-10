import importlib
import struct
import sys
import types
from pathlib import Path

import numpy as np

from conftest import load_module

pkg_root = Path(__file__).resolve().parents[1] / "pyspecdata"
_module_names = [
    "pyspecdata",
    "pyspecdata.core",
    "pyspecdata.datadir",
    "pyspecdata.general_functions",
    "pyspecdata.load_files",
    "pyspecdata.load_files.bruker_esr",
    "pyspecdata.load_files.load_cary",
]
_missing = object()
_saved_modules = {
    name: sys.modules.get(name, _missing) for name in _module_names
}

pyspecdata_pkg = types.ModuleType("pyspecdata")
pyspecdata_pkg.__path__ = [str(pkg_root)]
sys.modules["pyspecdata"] = pyspecdata_pkg


class _DummyConfig:
    def get_setting(self, *args, **kwargs):
        return None

    def set_setting(self, *args, **kwargs):
        pass


datadir_stub = types.ModuleType("pyspecdata.datadir")
datadir_stub.pyspec_config = _DummyConfig()
datadir_stub.rclone_search = lambda *args, **kwargs: None
sys.modules["pyspecdata.datadir"] = datadir_stub

load_files_pkg = types.ModuleType("pyspecdata.load_files")
load_files_pkg.__path__ = [str(pkg_root / "load_files")]
sys.modules["pyspecdata.load_files"] = load_files_pkg

sys.modules.pop("pyspecdata.core", None)
sys.modules.pop("pyspecdata.load_files.bruker_esr", None)
sys.modules.pop("pyspecdata.load_files.load_cary", None)

load_module("general_functions")
core = load_module("core", use_real_pint=True, use_real_h5py=True)
<<<<<<< Updated upstream
||||||| Stash base
# {{{ loading only portions of pyspecdata that are actually used
lmfitdata_mod = load_module(
    "lmfitdata", use_real_pint=True, use_real_h5py=True
)
pyspecdata_pkg.nddata = core.nddata
pyspecdata_pkg.ndshape = core.ndshape
pyspecdata_pkg.lmfitdata = lmfitdata_mod.lmfitdata
# }}}
=======
# {{{ loading only portions of pyspecdata that are actually used
pyspecdata_pkg.nddata = core.nddata
pyspecdata_pkg.ndshape = core.ndshape
pyspecdata_pkg.lmfitdata = load_module(
    "lmfitdata", use_real_pint=True, use_real_h5py=True
).lmfitdata
# }}}
>>>>>>> Stashed changes
bruker_esr = importlib.import_module("pyspecdata.load_files.bruker_esr")
load_cary = importlib.import_module("pyspecdata.load_files.load_cary")

for _name, _module in _saved_modules.items():
    if _module is _missing:
        sys.modules.pop(_name, None)
    else:
        sys.modules[_name] = _module


def _write_xepr_descriptor(path, *body_lines):
    path.write_text("#DESC\n" + "\n".join(body_lines) + "\n", encoding="utf-8")


def _write_big_endian_array(path, values, dtype):
    arr = np.asarray(values, dtype=np.dtype(dtype))
    path.write_bytes(arr.tobytes())

# SINGLE_USE_EXCEPTION
def _write_fake_cary_file(path, x_values, y_values, spectrum_name="fake_uv"):
    x_values = np.asarray(x_values, dtype="<f4")
    y_values = np.asarray(y_values, dtype="<f4")
    assert x_values.shape == y_values.shape
    record_dtype = np.dtype([("x", "<f4"), ("y", "<f4")])
    records = np.empty(len(x_values), dtype=record_dtype)
    records["x"] = x_values
    records["y"] = y_values
    header_dtype = np.dtype(
        [
            ("softversion", "<u4"),
            ("tmp1", "<u4"),
            ("end_x", "<i4"),
            ("start_x", "<i4"),
            ("V_min", "<i4"),
            ("V_max", "<i4"),
            ("points", "<u4"),
            ("tmp2", "<i4"),
            ("tmp3", "<u2"),
            ("tmp4", "<u2"),
            ("tmp5", "<i4"),
            ("spectrum_number", "<i4"),
            ("tmp6", "<f4"),
            ("tmp7", "<i4"),
            ("tmp8", "<f4"),
            ("tmp9", "<i4"),
            ("x_mode", "<f4"),
            ("y_mode", "<f4"),
        ]
    )
    magic = b"\x11Varian"
    prefix = bytes([len(magic)]) + magic + b"\x00" * (61 - len(magic))
    tstore_type = b"TContinuumStore"
    header = np.array(
        [
            (
                1,
                0,
                int(x_values[-1]),
                int(x_values[0]),
                0,
                0,
                len(x_values),
                0,
                0,
                0,
                0,
                1,
                0.0,
                0,
                0.0,
                0,
                0.0,
                0.0,
            )
        ],
        dtype=header_dtype,
    ).tobytes()
    name_bytes = spectrum_name.encode("ascii")[:255]
    name_bytes = name_bytes + b"\x00" + b"\x00" * (255 - len(name_bytes))
    padding = b"\x00" * (795 - (4 + len(tstore_type) + 4 + len(header)))
    block = b"".join(
        [
            struct.pack("<I", len(tstore_type)),
            tstore_type,
            struct.pack("<I", 1051 + records.nbytes),
            header,
            padding,
            name_bytes,
            records.tobytes(),
        ]
    )
    path.write_bytes(prefix + block)


def test_xepr_dsc_entrypoint_reads_chunkable_harmonics(tmp_path):
    dsc = tmp_path / "fake_esr.DSC"
    dta = tmp_path / "fake_esr.DTA"
    _write_xepr_descriptor(
        dsc,
        "IKKF REAL",
        "XPTS 4",
        "XWID 3",
        "XMIN 1",
        "XNAM Field",
        "XUNI 'G'",
        "Enable1stHarm 1",
        "Enable1stHarm90 1",
    )
    _write_big_endian_array(dta, np.arange(8, dtype=float), ">f8")

    data = bruker_esr.xepr(str(dsc))
    data.chunk_auto("harmonic", "phase")

    np.testing.assert_allclose(data.getaxis("phase"), np.array([0, 90]))
    np.testing.assert_allclose(data.getaxis("harmonic"), np.array([1]))
    np.testing.assert_allclose(data.data.shape, (4, 1, 2))


def test_xepr_dta_entrypoint_reads_power_axis_from_ygf(tmp_path):
    dsc = tmp_path / "fake_power.DSC"
    dta = tmp_path / "fake_power.DTA"
    ygf = tmp_path / "fake_power.YGF"
    _write_xepr_descriptor(
        dsc,
        "IKKF REAL",
        "XPTS 4",
        "XWID 3",
        "XMIN 1",
        "XNAM Field",
        "XUNI 'G'",
        "YPTS 3",
        "YTYP IGD",
        "YNAM 'Microwave Power'",
        "YUNI 'W'",
        "Enable1stHarm 1",
        "Enable1stHarm90 1",
    )
    _write_big_endian_array(dta, np.arange(24, dtype=float), ">f8")
    _write_big_endian_array(ygf, np.array([1e-3, 1e-2, 1e-1]), ">f8")

    data = bruker_esr.xepr(str(dta))
    data.chunk_auto("harmonic", "phase")
    microwave_power = data.getaxis("Microwave Power")

    np.testing.assert_allclose(
        microwave_power,
        np.array([-30.0, -20.0, -10.0]),
    )
    assert microwave_power.flags.writeable
    assert data.get_units("Microwave Power") == "dBm"
    data.set_axis("Microwave Power", lambda x: x + 1.0)
    np.testing.assert_allclose(
        data.getaxis("Microwave Power"),
        np.array([-29.0, -19.0, -9.0]),
    )


def test_cary_uv_loader_reads_fake_binary_file(tmp_path):
    path = tmp_path / "fake_uv.DSW"
    _write_fake_cary_file(
        path,
        x_values=np.array([500.0, 600.0, 700.0]),
        y_values=np.array([0.1, 0.2, 0.3]),
        spectrum_name="fake_uv",
    )

    loaded = load_cary.load_cary(str(path))
    assert len(loaded) == 1
    data = next(iter(loaded.values()))

    np.testing.assert_allclose(
        data.getaxis(r"$\lambda$"),
        np.array([500.0, 600.0, 700.0]),
    )
    np.testing.assert_allclose(data.data, np.array([0.1, 0.2, 0.3]))
    assert data.get_units(r"$\lambda$") == "nm"
