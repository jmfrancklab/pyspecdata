import os
import re
import sys
import importlib.machinery
import types
from pathlib import Path
import pytest

# allow the tests to import the local package without requiring an editable install
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# provide a stub for the compiled nnls extension so the package can be imported
_nnls_stub_module = types.ModuleType("_nnls")

def _nnls_stub(*args, **kwargs):
    raise RuntimeError("_nnls stub should not be executed in these tests")

_nnls_stub_module.nnls = _nnls_stub
_nnls_stub_module.nnls_regularized = _nnls_stub
_nnls_stub_module.nnls_regularized_loop = _nnls_stub
sys.modules["_nnls"] = _nnls_stub_module

pyspecdata_stub = types.ModuleType("pyspecdata")
pyspecdata_stub.__path__ = [str(Path(__file__).resolve().parents[1] / "pyspecdata")]
pyspecdata_stub.__spec__ = importlib.machinery.ModuleSpec(
    "pyspecdata", loader=None, is_package=True
)
sys.modules.setdefault("pyspecdata", pyspecdata_stub)

for _name in ["bruker_nmr", "prospa", "bruker_esr", "acert", "load_cary"]:
    sys.modules.setdefault(
        "pyspecdata.load_files." + _name, types.ModuleType(_name)
    )

_acert_stub_module = sys.modules["pyspecdata.load_files.acert"]

def _acert_passthrough(data, **kwargs):
    return data

for _attr in [
    "postproc_eldor_old",
    "postproc_eldor_3d",
    "postproc_generic",
    "postproc_echo_T2",
    "postproc_B1_se",
    "postproc_cw",
]:
    setattr(_acert_stub_module, _attr, _acert_passthrough)

open_subpath_stub = types.ModuleType("open_subpath")

def _open_subpath(*args, **kwargs):
    return False

open_subpath_stub.open_subpath = _open_subpath
sys.modules.setdefault("pyspecdata.load_files.open_subpath", open_subpath_stub)

zenodo_stub = types.ModuleType("zenodo")

def _zenodo_download(*args, **kwargs):
    raise RuntimeError("zenodo download should not run during these tests")

zenodo_stub.zenodo_download = _zenodo_download
sys.modules.setdefault("pyspecdata.load_files.zenodo", zenodo_stub)

core_stub = types.ModuleType("core")
core_stub.nddata_hdf5 = object()
sys.modules.setdefault("pyspecdata.core", core_stub)
sys.modules.setdefault("h5py", types.ModuleType("h5py"))

from pyspecdata.load_files import search_filename
from pyspecdata import datadir


def test_search_filename_respects_anchors(tmp_path, monkeypatch):
    base_dir = tmp_path / "experiment"
    base_dir.mkdir()
    (base_dir / "alpha.txt").write_text("alpha")
    (base_dir / "beta.dat").write_text("beta")
    (base_dir / "gamma.dat").write_text("gamma")

    def fake_getdatadir(exp_type=None):
        assert exp_type == "test_exp"
        return str(base_dir) + os.path.sep

    recorded = {}

    def fake_rclone(pattern, exp_type, directory):
        recorded["pattern"] = pattern
        recorded["exp_type"] = exp_type
        recorded["directory"] = directory
        return "cmd"

    monkeypatch.setattr("pyspecdata.load_files.getDATADIR", fake_getdatadir)
    monkeypatch.setattr("pyspecdata.load_files.rclone_search", fake_rclone)

    # ensure the anchored regex still finds the intended file locally
    assert sorted(os.listdir(str(base_dir))) == ["alpha.txt", "beta.dat", "gamma.dat"]
    check_pattern = re.compile(r"gamma\.dat$")
    assert check_pattern.search("gamma.dat")
    results = search_filename(r"gamma\.dat$", "test_exp", print_result=False)
    assert results == [os.path.join(str(base_dir), "gamma.dat")]
    assert recorded == {}


def test_search_filename_passes_raw_regex_to_rclone(tmp_path, monkeypatch):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    recorded = {}

    def fake_getdatadir(exp_type=None):
        assert exp_type == "missing_exp"
        return str(empty_dir) + os.path.sep

    def fake_rclone(pattern, exp_type, directory):
        recorded["pattern"] = pattern
        recorded["exp_type"] = exp_type
        recorded["directory"] = directory
        return "cmd"

    monkeypatch.setattr("pyspecdata.load_files.getDATADIR", fake_getdatadir)
    monkeypatch.setattr("pyspecdata.load_files.rclone_search", fake_rclone)

    # when nothing matches locally the regex should reach rclone unchanged
    with pytest.raises(RuntimeError):
        search_filename(r"delta$", "missing_exp", print_result=False)

    assert recorded["pattern"] == r"delta$"
    assert recorded["exp_type"] == "missing_exp"
    assert recorded["directory"] == str(empty_dir) + os.path.sep


def test_rclone_search_uses_regex_mode(monkeypatch, tmp_path):
    exp_type = "remote_exp"
    exp_key = datadir.PureWindowsPath(exp_type).as_posix().casefold()
    datadir.pyspec_config.config_vars["RcloneRemotes"][exp_key] = "example:remote"
    captured = {}

    def fake_system(command):
        captured["cmd"] = command
        return 0

    monkeypatch.setattr(datadir.os, "system", fake_system)

    destination = tmp_path / "local_copy"
    destination.mkdir()
    pattern = r"omega\.dat$"

    # the generated rclone command should use regex mode without adding wildcards
    command = datadir.rclone_search(pattern, exp_type, str(destination))

    assert command == captured["cmd"]
    assert '--include "{{' + pattern + '}}"' in command
    assert '--include "*' not in command
    assert "example:remote" in command
    assert str(destination) in command
