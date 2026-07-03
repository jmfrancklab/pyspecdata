"""Load CIQTEK JSON ``.epr`` files."""

import json
import logging

import numpy as np

from ..core import nddata
from ..general_functions import strm

logger = logging.getLogger("pyspecdata.load_files.ciqtek")

b0_texstr = r"$B_0$"


def _as_two_column_array(values, key, line_name):
    try:
        array = np.asarray(values, dtype=np.double)
    except Exception as exc:
        raise ValueError(
            strm("CIQTEK", key, "for", line_name, "is not numeric")
        ) from exc
    if array.ndim != 2 or array.shape[1] != 2:
        raise ValueError(
            strm(
                "CIQTEK",
                key,
                "for",
                line_name,
                "must be a two-column array",
            )
        )
    return array


def is_ciqtek_payload(payload):
    """Return true when *payload* has the CIQTEK JSON data layout."""
    if not isinstance(payload, dict):
        return False
    data_store = payload.get("dataStore")
    if not isinstance(data_store, dict):
        return False
    line_data = data_store.get("lineDataList")
    if not isinstance(line_data, list) or len(line_data) == 0:
        return False
    first_line = line_data[0]
    if not isinstance(first_line, dict):
        return False
    if "ReData" not in first_line or "ImData" not in first_line:
        return False
    try:
        _as_two_column_array(
            first_line["ReData"], "ReData", first_line.get("name", "line 0")
        )
        _as_two_column_array(
            first_line["ImData"], "ImData", first_line.get("name", "line 0")
        )
    except ValueError:
        return False
    return True


def is_ciqtek_file(filename):
    """Return true when *filename* is a CIQTEK JSON ``.epr`` file."""
    try:
        with open(filename, encoding="utf-8") as fp:
            payload = json.load(fp)
    except Exception:
        return False
    return is_ciqtek_payload(payload)


def _load_payload(filename):
    with open(filename, encoding="utf-8") as fp:
        payload = json.load(fp)
    if not is_ciqtek_payload(payload):
        raise ValueError(strm(filename, "does not look like a CIQTEK .epr"))
    return payload


def _common_param_axis(line_params):
    if len(line_params) == 0:
        return None, None
    common_keys = set(line_params[0].keys())
    for params in line_params[1:]:
        common_keys &= set(params.keys())
    if "delay" in common_keys:
        axis_name = "delay"
    elif len(common_keys) == 1:
        axis_name = next(iter(common_keys))
    else:
        return None, None
    raw_values = [params[axis_name] for params in line_params]
    try:
        axis_values = np.asarray(raw_values, dtype=np.double)
    except Exception:
        axis_values = np.asarray(raw_values)
    return axis_name, axis_values


def _load_lines(payload):
    line_data = payload["dataStore"]["lineDataList"]
    line_names = []
    line_params = []
    line_freqs = []
    signals = []
    field_axes = []
    g_axes = []
    have_g_axis = True
    for line_number, line in enumerate(line_data):
        if not isinstance(line, dict):
            raise ValueError(
                strm("CIQTEK line", line_number, "is not a dictionary")
            )
        line_name = line.get("name", "line_%d" % line_number)
        re_data = _as_two_column_array(line.get("ReData"), "ReData", line_name)
        im_data = _as_two_column_array(line.get("ImData"), "ImData", line_name)
        if re_data.shape != im_data.shape:
            raise ValueError(
                strm("CIQTEK ReData and ImData shapes differ for", line_name)
            )
        if not np.allclose(re_data[:, 0], im_data[:, 0]):
            raise ValueError(
                strm(
                    "CIQTEK ReData and ImData field axes differ for",
                    line_name,
                )
            )
        field_axes.append(re_data[:, 0])
        signals.append(re_data[:, 1] + 1j * im_data[:, 1])
        line_names.append(line_name)
        line_params.append(line.get("params", {}))
        line_freqs.append(line.get("freq"))
        if "gData" in line:
            g_data = _as_two_column_array(line["gData"], "gData", line_name)
            if g_data.shape[0] != re_data.shape[0]:
                raise ValueError(
                    strm("CIQTEK gData length differs for", line_name)
                )
            g_axes.append(g_data[:, 0])
        else:
            have_g_axis = False
    field_axis = field_axes[0]
    for line_name, this_axis in zip(line_names[1:], field_axes[1:]):
        if len(this_axis) != len(field_axis) or not np.allclose(
            this_axis, field_axis
        ):
            raise ValueError(
                strm(
                    "CIQTEK field axis for",
                    line_name,
                    "does not match the first line",
                )
            )
    g_axis = None
    if have_g_axis and len(g_axes) == len(line_data):
        g_axis = g_axes[0]
        if not all(
            len(this_axis) == len(g_axis) and np.allclose(this_axis, g_axis)
            for this_axis in g_axes[1:]
        ):
            g_axis = np.vstack(g_axes).T
    return line_names, line_params, line_freqs, field_axis, g_axis, signals


def load_ciqtek(filename):
    """Load a CIQTEK JSON ``.epr`` file as an :class:`nddata` object."""
    payload = _load_payload(filename)
    (
        line_names,
        line_params,
        line_freqs,
        field_axis,
        g_axis,
        signals,
    ) = _load_lines(payload)
    if len(signals) == 1:
        data = nddata(signals[0], [len(field_axis)], [b0_texstr])
        data.labels([b0_texstr], [field_axis])
    else:
        indirect_name, indirect_axis = _common_param_axis(line_params)
        if indirect_name is None:
            indirect_name = "indirect"
            indirect_axis = np.r_[0 : len(signals)]
        signal_array = np.vstack(signals).T
        data = nddata(
            signal_array,
            [len(field_axis), len(signals)],
            [b0_texstr, indirect_name],
        )
        data.labels([b0_texstr, indirect_name], [field_axis, indirect_axis])
    data.set_units(b0_texstr, "G")
    props = {
        "name": payload.get("name"),
        "source_format": "CIQTEK JSON EPR",
        "ciqtek_data_version": payload.get("_dataVersion"),
        "ciqtek_type": payload.get("type"),
        "ciqtek_device": payload.get("devicename"),
        "ciqtek_filename": payload.get("filename"),
        "ciqtek_create_time": payload.get("createTime"),
        "ciqtek_settings": payload.get("setting", {}),
        "ciqtek_line_names": line_names,
        "ciqtek_line_params": line_params,
        "ciqtek_line_freqs": line_freqs,
    }
    for prop_name, prop_value in props.items():
        if prop_value is not None:
            data.set_prop(prop_name, prop_value)
    if g_axis is not None:
        data.set_prop("ciqtek_g_axis", g_axis)
    return data
