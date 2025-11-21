from ..general_functions import *
from pylab import *
import h5py

logger = logging.getLogger("pyspecdata.hdf_save_dict_to_group")


def encode_list(name, seq, group, use_pytables_hack=False):
    """Write a Python sequence to HDF5.

    When ``use_pytables_hack`` is true, lists are written using the legacy
    ``LISTELEMENTS`` record attribute for compatibility with older layout
    expectations.  Otherwise, the sequence is expanded into a subgroup marked
    with ``LIST_NODE`` so the order can be reconstructed without the numpy
    record hack.
    """

    # dimlabels are always stored directly as an attribute so they mirror the
    # historical layout expected by the tests that inspect the HDF5 structure
    if name == "dimlabels":
        group.attrs[name] = np.array(seq)
        return

    if use_pytables_hack:
        if (
            isinstance(seq, (np.void, np.ndarray))
            and getattr(seq, "dtype", None) is not None
            and getattr(seq.dtype, "names", None) == ("LISTELEMENTS",)
        ):
            group.attrs[name] = seq
            return
        elements = list(seq)
        if len(elements) > 0 and isinstance(elements[0], str):
            elements = [x.encode("utf-8") for x in elements]
        rec = np.rec.fromarrays([elements], names="LISTELEMENTS")
        group.attrs[name] = rec
        return
    if (
        isinstance(seq, (np.void, np.ndarray))
        and getattr(seq, "dtype", None) is not None
        and getattr(seq.dtype, "names", None) == ("LISTELEMENTS",)
    ):
        group.attrs[name] = seq
        return
    target_group = group.create_group(name)
    target_group.attrs["LIST_NODE"] = True
    if isinstance(seq, tuple):
        target_group.attrs["LIST_CLASS"] = "tuple"
    for idx, entry in enumerate(seq):
        item_name = "ITEM" + str(idx)
        if isinstance(entry, (list, tuple)) or (
            isinstance(entry, (np.void, np.ndarray))
            and getattr(entry, "dtype", None) is not None
            and getattr(entry.dtype, "names", None) == ("LISTELEMENTS",)
        ):
            encode_list(item_name, entry, target_group, use_pytables_hack)
        elif issubclass(type(entry), np.ndarray):
            target_group.create_dataset(item_name, data=entry, dtype=entry.dtype)
        elif issubclass(type(entry), dict):
            subgroup = target_group.create_group(item_name)
            hdf_save_dict_to_group(subgroup, entry, use_pytables_hack)
        else:
            if entry is None:
                continue
            if isinstance(entry, str):
                target_group.attrs[item_name] = entry.encode("utf-8")
            elif np.isscalar(entry):
                target_group.attrs[item_name] = entry
            elif hasattr(entry, "__len__") and len(entry) > 0:
                target_group.attrs[item_name] = [
                    x.encode("utf-8") if isinstance(x, str) else x for x in entry
                ]


def decode_list(source):
    """Return a Python sequence from either list representation.

    ``source`` may be a subgroup marked with ``LIST_NODE`` or a legacy
    ``LISTELEMENTS`` numpy record.  Both forms are converted back to standard
    Python lists or tuples.
    """

    if isinstance(source, (np.void, np.ndarray)) and getattr(source, "dtype", None) is not None:
        if source.dtype.names == ("LISTELEMENTS",):
            decoded_attr = []
            for element in source["LISTELEMENTS"].flat:
                if isinstance(element, bytes):
                    decoded_attr.append(element.decode("utf-8"))
                else:
                    decoded_attr.append(element)
            return np.array(decoded_attr, dtype=object).reshape(
                source["LISTELEMENTS"].shape
            ).tolist()
    if hasattr(source, "attrs") and "LIST_NODE" in source.attrs and source.attrs["LIST_NODE"]:
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
                elif hasattr(source[item_name], "attrs") and "LIST_NODE" in source[
                    item_name
                ].attrs and source[item_name].attrs["LIST_NODE"]:
                    reconstructed.append(decode_list(source[item_name]))
                else:
                    reconstructed.append(hdf_load_dict_from_group(source[item_name]))
            else:
                if isinstance(source.attrs[item_name], bytes):
                    reconstructed.append(source.attrs[item_name].decode("utf-8"))
                elif hasattr(source.attrs[item_name], "__len__") and not np.isscalar(
                    source.attrs[item_name]
                ):
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


def hdf_save_dict_to_group(group, data, use_pytables_hack=False):
    """
    Copied as-is from ACERT hfesr code
    All numpy arrays are datasets.
    """
    for k, v in data.items():
        if issubclass(type(v), np.ndarray):
            logger.debug("Dataset type %s" % str(v.dtype))
            logger.debug("Adding %s=%s as dataset" % (k, v))
            group.create_dataset(k, data=v, dtype=v.dtype)
        elif isinstance(v, (list, tuple)) or (
            isinstance(v, (np.void, np.ndarray))
            and getattr(v, "dtype", None) is not None
            and getattr(v.dtype, "names", None) == ("LISTELEMENTS",)
        ):
            encode_list(k, v, group, use_pytables_hack)
        elif issubclass(type(v), dict):
            if set(v.keys()) == {"LISTELEMENTS"}:
                encode_list(
                    k,
                    np.rec.fromarrays([v["LISTELEMENTS"]], names="LISTELEMENTS")[0],
                    group,
                    use_pytables_hack,
                )
            else:
                subgroup = group.create_group(k)
                hdf_save_dict_to_group(subgroup, v, use_pytables_hack)
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
                # added encoding in following
                group.attrs[k] = [
                    x.encode("utf-8") if isinstance(x, str) else x for x in v
                ]

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
        if (
            hasattr(group.attrs[k], "dtype")
            and getattr(group.attrs[k].dtype, "names", None) == ("LISTELEMENTS",)
        ):
            retval[k] = decode_list(group.attrs[k])
        elif isinstance(group.attrs[k], bytes):
            retval[k] = group.attrs[k].decode("utf-8")
        elif hasattr(group.attrs[k], "__len__") and not np.isscalar(group.attrs[k]):
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
