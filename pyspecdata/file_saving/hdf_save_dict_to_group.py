import h5py
import numpy as np
import logging

logger = logging.getLogger("pyspecdata.hdf_save_dict_to_group")


def encode_list(name, input_list, use_pytables_hack=False):
    """Return a dictionary description of a Python sequence.

    The returned dictionary may include ``ATTRIBUTES_OF_MAIN_TREE`` entries to
    describe attributes that must be attached after the primary HDF tree is
    built.  This keeps the encoding logic independent from any particular HDF
    API while still conveying the desired layout.
    """

    attr_values = []
    attr_locations = []
    tree = {}

    # dimlabels must remain attributes, so they are requested via the special
    # post-tree metadata list
    if name == "dimlabels":
        attr_values.append({name: np.array(input_list)})
        attr_locations.append([])
    elif use_pytables_hack:
        # request the legacy LISTELEMENTS attribute on the current node
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
        attr_values.append({name: rec})
        attr_locations.append([])
    elif (
        isinstance(input_list, (np.void, np.ndarray))
        and getattr(input_list, "dtype", None) is not None
        and getattr(input_list.dtype, "names", None) == ("LISTELEMENTS",)
    ):
        attr_values.append({name: input_list})
        attr_locations.append([])
    else:
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
                encoded_entry = encode_list(
                    item_name, entry, use_pytables_hack
                )
                if item_name in encoded_entry:
                    list_group[item_name] = encoded_entry[item_name]
                if "ATTRIBUTES_OF_MAIN_TREE" in encoded_entry:
                    for idx_attr in range(
                        len(encoded_entry["ATTRIBUTES_OF_MAIN_TREE"])
                    ):
                        attr_values.append(
                            encoded_entry["ATTRIBUTES_OF_MAIN_TREE"][idx_attr]
                        )
                        attr_locations.append(
                            [name]
                            + encoded_entry[
                                "ATTRIBUTES_OF_MAIN_TREE_LOCATIONS"
                            ][idx_attr]
                        )
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

    if len(attr_values) > 0:
        tree["ATTRIBUTES_OF_MAIN_TREE"] = attr_values
        tree["ATTRIBUTES_OF_MAIN_TREE_LOCATIONS"] = attr_locations
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
        indices = sorted(
            {
                int(name[4:])
                for name in list(source.keys()) + list(source.attrs.keys())
                if name.startswith("ITEM")
            }
        )
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
                        hdf_load_dict_from_group(source[item_name])
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
    base_path=None,
    attr_values=None,
    attr_locations=None,
):
    """
    Copied as-is from ACERT hfesr code
    All numpy arrays are datasets.
    """

    if base_path is None:
        base_path = []
    top_call = attr_values is None and attr_locations is None
    if top_call:
        attr_values = []
        attr_locations = []

    if "ATTRIBUTES_OF_MAIN_TREE" in data:
        for idx in range(len(data["ATTRIBUTES_OF_MAIN_TREE"])):
            attr_values.append(data["ATTRIBUTES_OF_MAIN_TREE"][idx])
            attr_locations.append(
                base_path + data["ATTRIBUTES_OF_MAIN_TREE_LOCATIONS"][idx]
            )

    for k, v in data.items():
        if (
            k == "ATTRIBUTES_OF_MAIN_TREE"
            or k == "ATTRIBUTES_OF_MAIN_TREE_LOCATIONS"
        ):
            continue
        if issubclass(type(v), np.ndarray):
            logger.debug("Dataset type %s" % str(v.dtype))
            logger.debug("Adding %s=%s as dataset" % (k, v))
            group.create_dataset(k, data=v, dtype=v.dtype)
        elif isinstance(v, (list, tuple)):
            encoded = encode_list(k, v, use_pytables_hack)
            if k in encoded:
                subgroup = group.create_group(k)
                hdf_save_dict_to_group(
                    subgroup,
                    encoded[k],
                    use_pytables_hack,
                    base_path + [k],
                    attr_values,
                    attr_locations,
                )
            if "ATTRIBUTES_OF_MAIN_TREE" in encoded:
                for idx in range(len(encoded["ATTRIBUTES_OF_MAIN_TREE"])):
                    attr_values.append(encoded["ATTRIBUTES_OF_MAIN_TREE"][idx])
                    attr_locations.append(
                        base_path
                        + encoded["ATTRIBUTES_OF_MAIN_TREE_LOCATIONS"][idx]
                    )
        elif issubclass(type(v), dict):
            subgroup = group.create_group(k)
            hdf_save_dict_to_group(
                subgroup,
                v,
                use_pytables_hack,
                base_path + [k],
                attr_values,
                attr_locations,
            )
        else:
            if v is None:
                continue
            if isinstance(v, str):
                logger.debug("Adding %s=%s as string attribute" % (k, v))
                group.attrs[k] = v.encode("utf-8")
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
                    x.encode("utf-8") if isinstance(x, str) else x for x in v
                ]

    if top_call:
        for idx in range(len(attr_values)):
            target = group
            for path_item in attr_locations[idx]:
                target = target[path_item]
            for attr_name in attr_values[idx]:
                if isinstance(attr_values[idx][attr_name], str):
                    target.attrs[attr_name] = attr_values[idx][
                        attr_name
                    ].encode("utf-8")
                elif (
                    hasattr(attr_values[idx][attr_name], "__len__")
                    and not np.isscalar(attr_values[idx][attr_name])
                    and not isinstance(
                        attr_values[idx][attr_name], (np.void, np.ndarray)
                    )
                ):
                    converted = []
                    for element in attr_values[idx][attr_name]:
                        if isinstance(element, str):
                            converted.append(element.encode("utf-8"))
                        else:
                            converted.append(element)
                    target.attrs[attr_name] = converted
                else:
                    target.attrs[attr_name] = attr_values[idx][attr_name]


def hdf_load_dict_from_group(group):
    """Recursively load an HDF5 group into a plain ``dict``.

    This mirrors :func:`hdf_save_dict_to_group` and returns a dictionary
    representation of ``group`` suitable for :meth:`nddata.__setstate__`.
    """
    if "LIST_NODE" in group.attrs and group.attrs["LIST_NODE"]:
        return decode_list(group)
    retval = {}
    for k, v in group.items():
        if isinstance(v, h5py.Dataset):
            retval[k] = v[()]
        elif "LIST_NODE" in v.attrs and v.attrs["LIST_NODE"]:
            retval[k] = decode_list(v)
        else:
            retval[k] = hdf_load_dict_from_group(v)
    for k in group.attrs.keys():
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
    return retval
