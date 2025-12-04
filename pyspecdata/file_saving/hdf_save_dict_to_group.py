import h5py
import numpy as np
import logging

logger = logging.getLogger("pyspecdata.hdf_save_dict_to_group")


def encode_list(name, input_list, use_pytables_hack=False):
    """Return a dictionary description of a Python sequence.

    Lists are stored as subgroups marked with ``LIST_NODE`` and numbered
    ``ITEM`` entries.  Tuples are preserved by a ``LIST_CLASS`` attribute.  If
    ``use_pytables_hack`` is requested, the legacy ``LISTELEMENTS`` attribute
    is emitted instead.
    """

    tree = {}

    if use_pytables_hack:
        serialized_elements = []
        for entry in input_list:
            if isinstance(entry, str):
                serialized_elements.append(entry.encode("utf-8"))
            else:
                serialized_elements.append(entry)
        arr = np.array(serialized_elements)
        if arr.dtype.kind == "U":
            arr = np.array([x.encode("utf-8") for x in arr.flat]).reshape(
                arr.shape
            )
        rec = np.zeros(1, dtype=[("LISTELEMENTS", arr.dtype, arr.shape)])[0]
        rec["LISTELEMENTS"] = arr
        tree[name] = rec
        return tree

    if (
        isinstance(input_list, (np.void, np.ndarray))
        and getattr(input_list, "dtype", None) is not None
        and getattr(input_list.dtype, "names", None) == ("LISTELEMENTS",)
    ):
        tree[name] = input_list
        return tree

    list_group = {"LIST_NODE": True}
    if isinstance(input_list, tuple):
        list_group["LIST_CLASS"] = "tuple"
    for idx, entry in enumerate(input_list):
        item_name = "ITEM" + str(idx)
        if isinstance(entry, (list, tuple)) or (
            isinstance(entry, (np.void, np.ndarray))
            and getattr(entry, "dtype", None) is not None
            and getattr(entry.dtype, "names", None) == ("LISTELEMENTS",)
        ):
            encoded_entry = encode_list(item_name, entry, use_pytables_hack)
            if item_name in encoded_entry:
                list_group[item_name] = encoded_entry[item_name]
        elif issubclass(type(entry), np.ndarray):
            list_group[item_name] = entry
        elif issubclass(type(entry), dict):
            list_group[item_name] = entry
        else:
            if entry is None:
                continue
            if isinstance(entry, str):
                list_group[item_name] = entry
            elif np.isscalar(entry):
                list_group[item_name] = entry
            elif hasattr(entry, "__len__") and len(entry) > 0:
                converted = []
                for element in entry:
                    if isinstance(element, str):
                        converted.append(element.encode("utf-8"))
                    else:
                        converted.append(element)
                list_group[item_name] = converted
    tree[name] = list_group
    return tree


def decode_list(source):
    """Return a Python sequence from either list representation.

    ``source`` may be a subgroup marked with ``LIST_NODE`` or a legacy
    ``LISTELEMENTS`` numpy record.  Both forms are converted back to standard
    Python lists or tuples.
    """

    if (
        isinstance(source, (np.void, np.ndarray))
        and getattr(source, "dtype", None) is not None
    ):
        if source.dtype.names == ("LISTELEMENTS",):
            decoded_attr = []
            for element in source["LISTELEMENTS"].flat:
                if isinstance(element, bytes):
                    decoded_attr.append(element.decode("utf-8"))
                else:
                    decoded_attr.append(element)
            return (
                np.array(decoded_attr, dtype=object)
                .reshape(source["LISTELEMENTS"].shape)
                .tolist()
            )
    if (
        hasattr(source, "attrs")
        and "LIST_NODE" in source.attrs
        and source.attrs["LIST_NODE"]
    ):
        indices = sorted({
            int(name[4:])
            for name in list(source.keys()) + list(source.attrs.keys())
            if name.startswith("ITEM")
        })
        reconstructed = []
        for idx in indices:
            item_name = "ITEM" + str(idx)
            if item_name in source:
                if isinstance(source[item_name], h5py.Dataset):
                    reconstructed.append(source[item_name][()])
                elif (
                    hasattr(source[item_name], "attrs")
                    and "LIST_NODE" in source[item_name].attrs
                    and source[item_name].attrs["LIST_NODE"]
                ):
                    reconstructed.append(decode_list(source[item_name]))
                else:
                    reconstructed.append(
                        hdf_load_dict_from_group(
                            source[item_name], include_attrs=False
                        )
                    )
            else:
                if isinstance(source.attrs[item_name], bytes):
                    reconstructed.append(
                        source.attrs[item_name].decode("utf-8")
                    )
                elif hasattr(
                    source.attrs[item_name], "__len__"
                ) and not np.isscalar(source.attrs[item_name]):
                    decoded_attr = []
                    for element in source.attrs[item_name]:
                        if isinstance(element, bytes):
                            decoded_attr.append(element.decode("utf-8"))
                        else:
                            decoded_attr.append(element)
                    reconstructed.append(decoded_attr)
                else:
                    if isinstance(source.attrs[item_name], np.generic):
                        reconstructed.append(source.attrs[item_name].item())
                    else:
                        reconstructed.append(source.attrs[item_name])
        if "LIST_CLASS" in source.attrs and (
            source.attrs["LIST_CLASS"] == b"tuple"
            or source.attrs["LIST_CLASS"] == "tuple"
        ):
            return tuple(reconstructed)
        return reconstructed
    return source


def hdf_save_dict_to_group(
    group,
    data,
    use_pytables_hack=False,
):
    """
    Copied as-is from ACERT hfesr code
    All numpy arrays are datasets.
    """

    if (
        isinstance(data, (np.void, np.ndarray))
        and getattr(data, "dtype", None) is not None
        and getattr(data.dtype, "names", None) == ("LISTELEMENTS",)
    ):
        group.attrs["LISTELEMENTS"] = data["LISTELEMENTS"]
        return

    for k, v in data.items():
        if k == "dimlabels":
            encoded_labels = []
            for lbl in v:
                if isinstance(lbl, str):
                    encoded_labels.append((lbl.encode("utf-8"),))
                else:
                    encoded_labels.append((lbl,))
            group.attrs[k] = np.array(encoded_labels)
            continue
        if issubclass(type(v), np.ndarray):
            logger.debug("Dataset type %s" % str(v.dtype))
            logger.debug("Adding %s=%s as dataset" % (k, v))
            group.create_dataset(k, data=v, dtype=v.dtype)
        elif isinstance(v, (list, tuple)):
            encoded = encode_list(k, v, use_pytables_hack)
            if k in encoded:
                if (
                    isinstance(encoded[k], (np.void, np.ndarray))
                    and getattr(encoded[k], "dtype", None) is not None
                    and getattr(encoded[k].dtype, "names", None)
                    == ("LISTELEMENTS",)
                ):
                    group.attrs[k] = encoded[k]
                else:
                    subgroup = group.create_group(k)
                    hdf_save_dict_to_group(
                        subgroup,
                        encoded[k],
                        use_pytables_hack,
                    )
        elif issubclass(type(v), dict):
            if "NUMPY_DATA" in v:
                if use_pytables_hack:
                    subgroup = group.create_group(k)
                    if (
                        getattr(v["NUMPY_DATA"], "dtype", None) is not None
                        and getattr(v["NUMPY_DATA"].dtype, "names", None)
                        is not None
                    ):
                        for field_name in v["NUMPY_DATA"].dtype.names:
                            subgroup.create_dataset(
                                field_name,
                                data=v["NUMPY_DATA"][field_name],
                                dtype=v["NUMPY_DATA"][field_name].dtype,
                            )
                    else:
                        subgroup.create_dataset(
                            "data",
                            data=v["NUMPY_DATA"],
                            dtype=v["NUMPY_DATA"].dtype,
                        )
                    for attr_name in v:
                        if attr_name == "NUMPY_DATA":
                            continue
                        if v[attr_name] is None:
                            continue
                        if isinstance(v[attr_name], str):
                            subgroup.attrs[attr_name] = np.bytes_(v[attr_name])
                        elif np.isscalar(v[attr_name]):
                            subgroup.attrs[attr_name] = v[attr_name]
                        elif hasattr(v[attr_name], "__len__"):
                            converted = []
                            for element in v[attr_name]:
                                if isinstance(element, str):
                                    converted.append(np.bytes_(element))
                                else:
                                    converted.append(element)
                            subgroup.attrs[attr_name] = converted
                else:
                    dataset = group.create_dataset(
                        k, data=v["NUMPY_DATA"], dtype=v["NUMPY_DATA"].dtype
                    )
                    for attr_name in v:
                        if attr_name == "NUMPY_DATA":
                            continue
                        if v[attr_name] is None:
                            continue
                        if isinstance(v[attr_name], str):
                            dataset.attrs[attr_name] = np.bytes_(v[attr_name])
                        elif np.isscalar(v[attr_name]):
                            dataset.attrs[attr_name] = v[attr_name]
                        elif hasattr(v[attr_name], "__len__"):
                            converted = []
                            for element in v[attr_name]:
                                if isinstance(element, str):
                                    converted.append(np.bytes_(element))
                                else:
                                    converted.append(element)
                            dataset.attrs[attr_name] = converted
            else:
                subgroup = group.create_group(k)
                hdf_save_dict_to_group(
                    subgroup,
                    v,
                    use_pytables_hack,
                )
        else:
            if v is None:
                continue
            if isinstance(v, str):
                logger.debug("Adding %s=%s as string attribute" % (k, v))
                group.attrs[k] = np.bytes_(v)
            elif np.isscalar(v):
                logger.debug(
                    "Adding "
                    + repr(k)
                    + "="
                    + repr(v)
                    + " as scalar attribute"
                )
                group.attrs[k] = v
            elif hasattr(v, "__len__") and len(v) > 0:
                logger.debug(
                    "Adding " + repr(k) + "=" + repr(v) + " as list attribute"
                )
                group.attrs[k] = [
                    np.bytes_(x) if isinstance(x, str) else x for x in v
                ]


def hdf_load_dict_from_group(
    group,
    base_path=None,
    attr_values=None,
    attr_locations=None,
    include_attrs=True,
):
    """Recursively load an HDF5 group into a plain ``dict``."""
    if base_path is None:
        base_path = []
    if "LIST_NODE" in group.attrs and group.attrs["LIST_NODE"]:
        return decode_list(group)
    retval = {}
    for k, v in group.items():
        if isinstance(v, h5py.Dataset):
            retval[k] = {"NUMPY_DATA": v[()]}
            for attr_name in v.attrs:
                if hasattr(v.attrs[attr_name], "dtype") and getattr(
                    v.attrs[attr_name].dtype, "names", None
                ) == ("LISTELEMENTS",):
                    retval[k][attr_name] = decode_list(v.attrs[attr_name])
                elif isinstance(v.attrs[attr_name], bytes):
                    retval[k][attr_name] = v.attrs[attr_name].decode("utf-8")
                elif hasattr(
                    v.attrs[attr_name], "__len__"
                ) and not np.isscalar(v.attrs[attr_name]):
                    decoded_attr = []
                    for element in v.attrs[attr_name]:
                        if isinstance(element, bytes):
                            decoded_attr.append(element.decode("utf-8"))
                        else:
                            decoded_attr.append(element)
                    retval[k][attr_name] = decoded_attr
                else:
                    if isinstance(v.attrs[attr_name], np.generic):
                        retval[k][attr_name] = v.attrs[attr_name].item()
                    else:
                        retval[k][attr_name] = v.attrs[attr_name]
        elif "LIST_NODE" in v.attrs and v.attrs["LIST_NODE"]:
            retval[k] = decode_list(v)
        else:
            retval[k] = hdf_load_dict_from_group(
                v, base_path + [k], attr_values, attr_locations, include_attrs
            )
    attr_dict = {}
    for k in group.attrs.keys():
        if k == "dimlabels":
            decoded_labels = []
            for entry in group.attrs[k]:
                value = entry
                if isinstance(entry, np.ndarray):
                    if entry.size == 0:
                        continue
                    value = entry.flat[0]
                elif isinstance(entry, (list, tuple)):
                    if len(entry) == 0:
                        continue
                    value = entry[0]
                if isinstance(value, bytes):
                    value = value.decode("utf-8")
                decoded_labels.append((value,))
            retval[k] = decoded_labels
            attr_dict[k] = retval[k]
            continue
        if hasattr(group.attrs[k], "dtype") and getattr(
            group.attrs[k].dtype, "names", None
        ) == ("LISTELEMENTS",):
            retval[k] = decode_list(group.attrs[k])
        elif isinstance(group.attrs[k], bytes):
            retval[k] = group.attrs[k].decode("utf-8")
        elif hasattr(group.attrs[k], "__len__") and not np.isscalar(
            group.attrs[k]
        ):
            decoded_attr = []
            for element in group.attrs[k]:
                if isinstance(element, bytes):
                    decoded_attr.append(element.decode("utf-8"))
                else:
                    decoded_attr.append(element)
            retval[k] = decoded_attr
        else:
            if isinstance(group.attrs[k], np.generic):
                retval[k] = group.attrs[k].item()
            else:
                retval[k] = group.attrs[k]
        attr_dict[k] = retval[k]
    return retval
