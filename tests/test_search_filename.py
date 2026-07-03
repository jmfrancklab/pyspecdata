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


def test_load_indiv_file_passes_zenodo_to_xepr_without_companion_search(
    tmp_path, monkeypatch
):
    dsc = tmp_path / "sample.DSC"
    dsc.write_text("descriptor")
    calls = []
    expected = types.SimpleNamespace(dimlabels=[])

    def fake_search_filename(*_args, **_kwargs):
        raise AssertionError(
            "load_indiv_file should not search for XEPR companions"
        )

    monkeypatch.setattr(load_files, "search_filename", fake_search_filename)
    monkeypatch.setattr(load_files, "_check_signature", lambda _path: "TXT")
    monkeypatch.setattr(load_files, "_check_extension", lambda _path: "DSC")

    def fake_xepr(filename, dimname="", exp_type=None, zenodo=None):
        calls.append((filename, dimname, exp_type, zenodo))
        return expected

    monkeypatch.setattr(load_files.bruker_esr, "xepr", fake_xepr)

    result = load_files.load_indiv_file(
        str(dsc), exp_type="remote_exp", zenodo="21084153"
    )

    assert result is expected
    assert calls == [
        (str(dsc), "", "remote_exp", "21084153"),
    ]


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


def test_search_filename_skips_zenodo_when_file_is_local(
    tmp_path, monkeypatch, capsys
):
    base_dir = tmp_path / "experiment"
    base_dir.mkdir()
    (base_dir / "local.dat").write_text("local")

    def fake_getdatadir(exp_type=None):
        assert exp_type == "local_exp"
        return str(base_dir) + os.path.sep

    def fake_zenodo_download(*_args, **_kwargs):
        raise AssertionError("zenodo_download should not be called")

    monkeypatch.setattr(load_files, "getDATADIR", fake_getdatadir)
    monkeypatch.setattr(load_files, "zenodo_download", fake_zenodo_download)

    results = load_files.search_filename(
        r"local\.dat$", "local_exp", print_result=False, zenodo="21041480"
    )

    assert results == [str(base_dir) + os.path.sep + "local.dat"]
    assert capsys.readouterr().out == ""


def test_zenodo_download_prints_after_download(tmp_path, monkeypatch, capsys):
    zenodo = load_module("load_files.zenodo")
    base_dir = tmp_path / "experiment"
    base_dir.mkdir()
    file_key = "Pure_T177R1a_pR_210615.BSW"
    file_url = "https://zenodo.org/api/records/21041480/files/" + file_key

    def fake_getdatadir(exp_type=None):
        assert exp_type == "remote_exp"
        return str(base_dir) + os.path.sep

    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {
                "files": [
                    {
                        "key": file_key,
                        "links": {"self": file_url},
                    }
                ]
            }

    retrieved = {}

    def fake_get(url):
        retrieved["record_url"] = url
        return FakeResponse()

    def fake_urlretrieve(url, dest):
        retrieved["file_url"] = url
        retrieved["dest"] = dest
        with open(dest, "w") as fp:
            fp.write("remote")

    monkeypatch.setattr(zenodo, "getDATADIR", fake_getdatadir)
    monkeypatch.setattr(zenodo.requests, "get", fake_get)
    monkeypatch.setattr(zenodo.urllib.request, "urlretrieve", fake_urlretrieve)

    path = zenodo.zenodo_download(
        "21041480", r".*T177R1a_pR_210615.*", exp_type="remote_exp"
    )

    assert path == str(base_dir) + os.path.sep + file_key
    assert retrieved["record_url"] == "https://zenodo.org/api/records/21041480"
    assert retrieved["file_url"] == file_url
    assert retrieved["dest"] == path
    assert (
        capsys.readouterr().out
        == f"Downloaded from zenodo '{file_url}' to {path}\n"
    )


def test_search_filename_uses_zenodo_when_file_is_missing(
    tmp_path, monkeypatch
):
    base_dir = tmp_path / "experiment"
    base_dir.mkdir()

    def fake_getdatadir(exp_type=None):
        assert exp_type == "remote_exp"
        return str(base_dir) + os.path.sep

    calls = []

    def fake_zenodo_download(deposition, searchstring, exp_type=None):
        calls.append((deposition, searchstring, exp_type))
        downloaded = base_dir / "remote.dat"
        downloaded.write_text("remote")
        return str(downloaded)

    def fake_rclone(*_args, **_kwargs):
        raise AssertionError("rclone_search should not be called")

    monkeypatch.setattr(load_files, "getDATADIR", fake_getdatadir)
    monkeypatch.setattr(load_files, "zenodo_download", fake_zenodo_download)
    monkeypatch.setattr(load_files, "rclone_search", fake_rclone)

    results = load_files.search_filename(
        r"remote\.dat$",
        "remote_exp",
        print_result=False,
        zenodo="21041480",
    )

    assert results == [str(base_dir) + os.path.sep + "remote.dat"]
    assert calls == [("21041480", r"remote\.dat$", "remote_exp")]


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


def test_zenodo_upload_existing_deposition_uses_documented_file_api(
    tmp_path, monkeypatch
):
    zenodo = load_module("load_files.zenodo")
    local_path = tmp_path / "payload.dat"
    local_path.write_text("payload")
    settings = []
    captured = {}

    def fake_post(url, headers=None, data=None, files=None, **kwargs):
        captured["url"] = url
        captured["headers"] = headers
        captured["data"] = data
        captured["kwargs"] = kwargs
        captured["file_contents"] = files["file"].read()

        class FakeResponse:
            def raise_for_status(self):
                pass

            def json(self):
                return {"name": "payload.dat"}

        return FakeResponse()

    monkeypatch.setattr(zenodo, "_get_token", lambda: "TOKEN")
    monkeypatch.setattr(
        zenodo,
        "create_deposition",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("should use existing deposition id")
        ),
    )
    monkeypatch.setattr(zenodo.requests, "post", fake_post)
    monkeypatch.setattr(
        zenodo.pyspec_config,
        "get_setting",
        lambda *_args, **_kwargs: "0",
    )
    monkeypatch.setattr(
        zenodo.pyspec_config,
        "set_setting",
        lambda *args, **_kwargs: settings.append(args),
    )

    deposition_id = zenodo.zenodo_upload(
        str(local_path), deposition_id="21084153"
    )

    assert deposition_id == "21084153"
    assert captured["url"].endswith("/deposit/depositions/21084153/files")
    assert captured["headers"] == {"Authorization": "Bearer TOKEN"}
    assert captured["data"] == {"name": "payload.dat"}
    assert captured["kwargs"] == {}
    assert captured["file_contents"] == b"payload"
    assert ("zenodo", "upload_deposition1", "21084153") in settings


def test_pyspecdata_zenodo_cli_uses_existing_numeric_deposition_id(
    tmp_path, monkeypatch
):
    zenodo = load_module("load_files.zenodo")
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    local_path = data_dir / "one.dat"
    local_path.write_text("one")
    (tmp_path / "data_files.csv").write_text(
        "Filename,Path,exp_type\n"
        f"one.dat,{data_dir},example\n",
        encoding="utf-8",
    )
    calls = []

    def fake_upload(local_path, title=None, deposition_id=None):
        calls.append((local_path, title, deposition_id))
        return deposition_id

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(zenodo, "zenodo_upload", fake_upload)

    assert zenodo.cmd(["21084153"]) == 0
    assert calls == [(str(local_path), None, "21084153")]


def test_pyspecdata_zenodo_cli_uses_title_for_new_deposition(
    tmp_path, monkeypatch
):
    zenodo = load_module("load_files.zenodo")
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    first_path = data_dir / "one.dat"
    second_path = data_dir / "two.dat"
    first_path.write_text("one")
    second_path.write_text("two")
    (tmp_path / "data_files.csv").write_text(
        "Filename,Path,exp_type\n"
        f"one.dat,{data_dir},example\n"
        f"two.dat,{data_dir},example\n",
        encoding="utf-8",
    )
    calls = []

    def fake_upload(local_path, title=None, deposition_id=None):
        calls.append((local_path, title, deposition_id))
        if deposition_id is None:
            return "777777"
        return deposition_id

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(zenodo, "zenodo_upload", fake_upload)

    assert zenodo.cmd(["my new data"]) == 0
    assert calls == [
        (str(first_path), "my new data", None),
        (str(second_path), None, "777777"),
    ]


def test_pyspecdata_zenodo_cli_requires_one_argument():
    zenodo = load_module("load_files.zenodo")
    with pytest.raises(SystemExit):
        zenodo.cmd([])
    with pytest.raises(SystemExit):
        zenodo.cmd(["21084153", "extra"])


def test_pyspecdata_zenodo_cli_requires_data_files_csv(
    tmp_path, monkeypatch
):
    zenodo = load_module("load_files.zenodo")
    monkeypatch.chdir(tmp_path)

    with pytest.raises(FileNotFoundError, match="data_files.csv"):
        zenodo.cmd(["21084153"])
