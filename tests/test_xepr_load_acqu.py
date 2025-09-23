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


def test_ignore_blank_line(tmp_path):
    dsc = tmp_path / "blank.dsc"
    dsc.write_text(
        "#BLOCK\nVAR1 one\n\nVAR2 two\n"
    )
    result = xepr_load_acqu(str(dsc))
    assert result["BLOCK"]["VAR1"] == "one"
    assert result["BLOCK"]["VAR2"] == "two"


def test_multiline_parameter(tmp_path):
    dsc = tmp_path / "ppg.dsc"
    dsc.write_text(
        "#BLOCK\n"
        "PpgText Here is text belonging to the pulse program\\n\\\n"
        "Yes, this text belongs to the pulse program\\n\\\n"
        "Yes, indeed it does\\n\\\n"
        "\\n\n"
        "AnotherParam 20\n"
    )
    result = xepr_load_acqu(str(dsc))
    expected = (
        "Here is text belonging to the pulse program\n"
        "Yes, this text belongs to the pulse program\n"
        "Yes, indeed it does\n"
        "\n\n"
    )
    assert result["BLOCK"]["PpgText"] == expected
    assert result["BLOCK"]["AnotherParam"] == 20


def test_space_separated_string(tmp_path):
    dsc = tmp_path / "spaced.dsc"
    dsc.write_text(
        "#BLOCK\nPlsSPELEXPSlct     Echo separate cycles\n",
    )
    result = xepr_load_acqu(str(dsc))
    assert result["BLOCK"]["PlsSPELEXPSlct"] == "Echo separate cycles"


def test_multiline_numeric_text(tmp_path):
    dsc = tmp_path / "numeric.dsc"
    dsc.write_text(
        "#BLOCK\n"
        "PpgText start 180 \\\n"
        "number 90\\nnext \\\n"
        "final\n"
        "AnotherParam 20\n"
    )
    result = xepr_load_acqu(str(dsc))
    expected = "start 180 number 90\nnext final"
    assert result["BLOCK"]["PpgText"] == expected
    assert result["BLOCK"]["AnotherParam"] == 20
