import importlib
import sys
import types
from pathlib import Path

pkg_root = Path(__file__).resolve().parents[1] / "pyspecdata"
# create stub package to avoid heavy imports
pyspecdata_pkg = types.ModuleType("pyspecdata")
pyspecdata_pkg.__path__ = [str(pkg_root)]
sys.modules.setdefault("pyspecdata", pyspecdata_pkg)

core_stub = types.ModuleType("pyspecdata.core")
core_stub.nddata = object
sys.modules.setdefault("pyspecdata.core", core_stub)

gen_stub = types.ModuleType("pyspecdata.general_functions")
gen_stub.strm = lambda *args, **kwargs: ""
gen_stub.lsafen = lambda x: x
sys.modules.setdefault("pyspecdata.general_functions", gen_stub)

datadir_stub = types.ModuleType("pyspecdata.datadir")
datadir_stub.rclone_search = lambda *args, **kwargs: None
sys.modules.setdefault("pyspecdata.datadir", datadir_stub)

load_files_pkg = types.ModuleType("pyspecdata.load_files")
load_files_pkg.__path__ = [str(pkg_root / "load_files")]
sys.modules.setdefault("pyspecdata.load_files", load_files_pkg)

bruker_esr = importlib.import_module("pyspecdata.load_files.bruker_esr")
xepr_load_acqu = bruker_esr.xepr_load_acqu


def test_line_continuation(tmp_path):
    dsc = tmp_path / "test.dsc"
    dsc.write_text(
        "#BLOCK\n"
        "VAR firstline\\\n"
        "secondline\\nthird\\\n"
        "fourth\n"
    )
    result = xepr_load_acqu(str(dsc))
    expected = "firstlinesecondline\nthirdfourth"
    assert result["BLOCK"]["VAR"] == expected
