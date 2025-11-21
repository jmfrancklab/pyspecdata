from ..general_functions import *
from pylab import *
import h5py

logger = logging.getLogger("pyspecdata.hdf_save_dict_to_group")


def encode_list(name, seq, group):
    """Write a Python sequence to HDF5.

    Lists and tuples are stored as subgroups marked with ``LIST_NODE`` so that
    ordering is preserved without relying on the pytables ``LISTELEMENTS``
    record hack.  Legacy ``LISTELEMENTS`` records are still written directly as
    attributes when present to maintain backward compatibility.
    """

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
            encode_list(item_name, entry, target_group)
        elif issubclass(type(entry), np.ndarray):
            target_group.create_dataset(item_name, data=entry, dtype=entry.dtype)
        elif issubclass(type(entry), dict):
            subgroup = target_group.create_group(item_name)
            hdf_save_dict_to_group(subgroup, entry)
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


def decode_list(group):
    """Return a Python sequence from an HDF5 group marked ``LIST_NODE``.

    This reverses :func:`encode_list` by reading ``ITEM`` entries in order and
    honoring the stored ``LIST_CLASS`` or legacy ``LISTCLASS`` markers for
    tuples.
    """

    indices = sorted(
        {
            int(name[4:])
            for name in list(group.keys()) + list(group.attrs.keys())
            if name.startswith("ITEM")
        }
    )
    reconstructed = []
    for idx in indices:
        item_name = "ITEM" + str(idx)
        if item_name in group:
            if isinstance(group[item_name], h5py.Dataset):
                reconstructed.append(group[item_name][()])
            elif "LIST_NODE" in group[item_name].attrs and group[item_name].attrs[
                "LIST_NODE"
            ]:
                reconstructed.append(decode_list(group[item_name]))
            else:
                reconstructed.append(hdf_load_dict_from_group(group[item_name]))
        else:
            if isinstance(group.attrs[item_name], bytes):
                reconstructed.append(group.attrs[item_name].decode("utf-8"))
            elif hasattr(group.attrs[item_name], "__len__") and not np.isscalar(
                group.attrs[item_name]
            ):
                decoded_attr = []
                for element in group.attrs[item_name]:
                    if isinstance(element, bytes):
                        decoded_attr.append(element.decode("utf-8"))
                    else:
                        decoded_attr.append(element)
                reconstructed.append(decoded_attr)
            else:
                if isinstance(group.attrs[item_name], np.generic):
                    reconstructed.append(group.attrs[item_name].item())
                else:
                    reconstructed.append(group.attrs[item_name])
    if (
        "LIST_CLASS" in group.attrs
        and (group.attrs["LIST_CLASS"] == b"tuple" or group.attrs["LIST_CLASS"] == "tuple")
    ) or (
        "LISTCLASS" in group.attrs
        and (group.attrs["LISTCLASS"] == b"tuple" or group.attrs["LISTCLASS"] == "tuple")
    ):
        return tuple(reconstructed)
    return reconstructed


def hdf_save_dict_to_group(group, data):
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
            encode_list(k, v, group)
        elif issubclass(type(v), dict):
            if set(v.keys()) == {"LISTELEMENTS"}:
                encode_list(k, np.rec.fromarrays([v["LISTELEMENTS"]], names="LISTELEMENTS")[0], group)
            else:
                subgroup = group.create_group(k)
                hdf_save_dict_to_group(subgroup, v)
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
        else:
            retval[k] = hdf_load_dict_from_group(v)
    for k in group.attrs.keys():
        if (
            hasattr(group.attrs[k], "dtype")
            and getattr(group.attrs[k].dtype, "names", None) == ("LISTELEMENTS",)
        ):
            decoded_attr = []
            for element in group.attrs[k]["LISTELEMENTS"].flat:
                if isinstance(element, bytes):
                    decoded_attr.append(element.decode("utf-8"))
                else:
                    decoded_attr.append(element)
            retval[k] = np.array(decoded_attr, dtype=object).reshape(
                group.attrs[k]["LISTELEMENTS"].shape
            ).tolist()
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
