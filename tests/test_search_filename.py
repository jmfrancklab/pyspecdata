import os
import sys
import types
import pytest

from conftest import load_module

# stub out the optional _nnls extension so importing the package under test
# does not try to load a compiled module that is not built in the CI image
sys.modules.setdefault("_nnls", types.SimpleNamespace())

# load the modules under test using the shared helper so optional dependencies
# are stubbed the same way as the rest of the suite
datadir = load_module("datadir")
load_files = load_module("load_files.__init__")


def test_search_filename_respects_anchors(tmp_path, monkeypatch):
    base_dir = tmp_path / "experiment"
    base_dir.mkdir()
    (base_dir / "alpha.txt").write_text("alpha")
    (base_dir / "beta.dat").write_text("beta")
    (base_dir / "gamma.dat").write_text("gamma")

    # make sure the function resolves the experiment directory to the
    # temporary tree created above
    def fake_getdatadir(exp_type=None):
        assert exp_type == "test_exp"
        return str(base_dir) + os.path.sep

    # record whether the rclone fallback was triggered
    recorded = []

    def fake_rclone(pattern, exp_type, directory):
        recorded.append((pattern, exp_type, directory))
        return "cmd"

    monkeypatch.setattr(load_files, "getDATADIR", fake_getdatadir)
    monkeypatch.setattr(load_files, "rclone_search", fake_rclone)

    results = load_files.search_filename(
        r"gamma\.dat$", "test_exp", print_result=False
    )

    expected_path = str(base_dir) + os.path.sep + "gamma.dat"
    assert results == [expected_path]
    assert recorded == []


def test_search_filename_passes_raw_regex_to_rclone(tmp_path, monkeypatch):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()

    # the search should look inside the empty directory to confirm nothing is
    # present before falling back to rclone
    def fake_getdatadir(exp_type=None):
        assert exp_type == "missing_exp"
        return str(empty_dir) + os.path.sep

    recorded = {}

    def fake_rclone(pattern, exp_type, directory):
        recorded["pattern"] = pattern
        recorded["exp_type"] = exp_type
        recorded["directory"] = directory
        return "cmd"

    monkeypatch.setattr(load_files, "getDATADIR", fake_getdatadir)
    monkeypatch.setattr(load_files, "rclone_search", fake_rclone)

    with pytest.raises(RuntimeError):
        load_files.search_filename(
            r"delta$", "missing_exp", print_result=False
        )

    assert recorded["pattern"] == r"delta$"
    assert recorded["exp_type"] == "missing_exp"
    assert recorded["directory"] == str(empty_dir) + os.path.sep


def test_rclone_search_uses_regex_mode(monkeypatch, tmp_path):
    exp_type = "remote_exp"
    exp_key = datadir.PureWindowsPath(exp_type).as_posix().casefold()
    datadir.pyspec_config.config_vars["RcloneRemotes"][
        exp_key
    ] = "example:remote"

    captured = {}

    def fake_system(command):
        captured["cmd"] = command
        return 0

    monkeypatch.setattr(datadir.os, "system", fake_system)

    destination = tmp_path / "local_copy"
    destination.mkdir()
    pattern = r"omega\.dat$"

    # the generated command should enable regex matching without adding
    # wildcard padding around the original pattern
    command = datadir.rclone_search(pattern, exp_type, str(destination))

    assert command == captured["cmd"]
    assert '--include "{{' + pattern + '}}"' in command
    assert '--include "*' not in command
    assert "example:remote" in command
    assert str(destination) in command
