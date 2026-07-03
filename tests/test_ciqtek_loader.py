import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from pyspecdata.load_files import ciqtek, load_indiv_file

BNAME = r"$B_0$"
FIELDS = [3346.0, 3396.0, 3446.0]
G_AXIS = [2.03795, 2.00750, 1.97881]


def _line(name, delay, real_values, imag_values):
    return {
        "ImData": list(map(list, zip(FIELDS, imag_values))),
        "ReData": list(map(list, zip(FIELDS, real_values))),
        "freq": 9.54314,
        "gData": list(map(list, zip(G_AXIS, real_values))),
        "name": name,
        "params": {"delay": str(delay)} if delay is not None else {},
    }


def _payload(lines, experiment_type="CW EPR/2D Field-Delay Sweep"):
    return {
        "_dataVersion": "V2.0",
        "_loadExperimentData": False,
        "_modify": False,
        "code": 1,
        "createTime": "2026-02-19 17:46:56.127",
        "dataStore": {
            "lineDataList": lines,
            "xAxisName": "Field[G]",
            "yAxisName": "Intensity",
            "zAxisName": "",
        },
        "devicename": "EPR200_2.0",
        "filename": "ciqtek_test_file",
        "name": "CWEPR-test",
        "setting": {
            "frequency": 9.54314,
            "lineEdit_2DFieldDelaySweep_NoOfPoints": len(FIELDS),
        },
        "type": experiment_type,
    }


def _write_payload(path, payload):
    with open(path, "w", encoding="utf-8") as fp:
        json.dump(payload, fp)


def test_load_ciqtek_2d_complex_data_and_props(tmp_path):
    path = tmp_path / "field_delay.epr"
    _write_payload(
        path,
        _payload(
            [
                _line("Delay_0", 0, [0.0, 1.0, 0.2], [0.0, 0.2, 0.1]),
                _line("Delay_600", 600, [0.1, 0.8, 0.3], [0.1, 0.3, 0.0]),
                _line(
                    "Delay_600_repeat",
                    600,
                    [-0.1, 0.6, 0.4],
                    [0.0, 0.4, -0.1],
                ),
            ]
        ),
    )

    data = ciqtek.load_ciqtek(str(path))

    assert data.dimlabels == [BNAME, "delay"]
    np.testing.assert_allclose(data.getaxis(BNAME), FIELDS)
    np.testing.assert_allclose(data.getaxis("delay"), [0.0, 600.0, 600.0])
    np.testing.assert_allclose(
        data.data,
        np.array(
            [
                [0.0 + 0.0j, 0.1 + 0.1j, -0.1 + 0.0j],
                [1.0 + 0.2j, 0.8 + 0.3j, 0.6 + 0.4j],
                [0.2 + 0.1j, 0.3 + 0.0j, 0.4 - 0.1j],
            ]
        ),
    )
    assert data.get_units(BNAME) == "G"
    assert data.get_prop("name") == "CWEPR-test"
    assert data.get_prop("source_format") == "CIQTEK JSON EPR"
    assert data.get_prop("ciqtek_data_version") == "V2.0"
    assert data.get_prop("ciqtek_type") == "CW EPR/2D Field-Delay Sweep"
    assert data.get_prop("ciqtek_device") == "EPR200_2.0"
    assert data.get_prop("ciqtek_filename") == "ciqtek_test_file"
    assert data.get_prop("ciqtek_create_time") == "2026-02-19 17:46:56.127"
    assert data.get_prop("ciqtek_settings")["frequency"] == 9.54314
    assert data.get_prop("ciqtek_line_names") == [
        "Delay_0",
        "Delay_600",
        "Delay_600_repeat",
    ]
    assert data.get_prop("ciqtek_line_params") == [
        {"delay": "0"},
        {"delay": "600"},
        {"delay": "600"},
    ]
    assert data.get_prop("ciqtek_line_freqs") == [9.54314] * 3
    np.testing.assert_allclose(data.get_prop("ciqtek_g_axis"), G_AXIS)
    assert data.get_prop("ReData") is None
    assert data.get_prop("ImData") is None


def test_load_ciqtek_1d_field_sweep(tmp_path):
    path = tmp_path / "field_sweep.epr"
    _write_payload(
        path,
        _payload(
            [_line("FieldSweep", None, [0.0, 1.0, 0.2], [0.0, 0.2, 0.1])],
            experiment_type="CW EPR/1D Field Sweep",
        ),
    )

    data = ciqtek.load_ciqtek(str(path))

    assert data.dimlabels == [BNAME]
    np.testing.assert_allclose(data.getaxis(BNAME), FIELDS)
    np.testing.assert_allclose(data.data, [0.0 + 0.0j, 1.0 + 0.2j, 0.2 + 0.1j])
    assert data.get_prop("ciqtek_line_params") == [{}]


def test_load_indiv_file_dispatches_ciqtek_epr(tmp_path):
    path = tmp_path / "field_delay.epr"
    _write_payload(path, _payload([_line("Delay_0", 0, [0, 1, 0], [0, 0, 0])]))

    data = load_indiv_file(str(path))

    assert data.get_prop("source_format") == "CIQTEK JSON EPR"
    np.testing.assert_allclose(data.data, [0.0 + 0.0j, 1.0 + 0.0j, 0.0 + 0.0j])


def test_non_ciqtek_json_epr_is_not_claimed(tmp_path):
    path = tmp_path / "not_ciqtek.epr"
    path.write_text('{"hello": "world"}', encoding="utf-8")

    assert not ciqtek.is_ciqtek_file(str(path))
    with pytest.raises(RuntimeError):
        load_indiv_file(str(path))


def test_ciqtek_gallery_example_runs(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    example = repo_root / "examples/ESR/ciqtek_example.py"
    local_data = (
        Path.home()
        / "exp_data"
        / "ciqtek"
        / "384K_1mM_TEMPO_SFO5_pt5G_MAmp_25dB_100G_22scans_2D.epr"
    )
    if not local_data.exists():
        pytest.skip("local CIQTEK example data is not present")
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    env["MPLCONFIGDIR"] = str(tmp_path)
    env["PYTHONPATH"] = str(repo_root)
    subprocess.run(
        [sys.executable, str(example)],
        cwd=repo_root,
        env=env,
        check=True,
    )
