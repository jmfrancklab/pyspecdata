from ..general_functions import *
from pylab import *
import h5py

logger = logging.getLogger("pyspecdata.hdf_save_dict_to_group")


def hdf_save_dict_to_group(group, data):
    """
    Copied as-is from ACERT hfesr code
    All numpy arrays are datasets.
    """

    def save_sequence(seq, target_group):
        """Write a Python list or tuple to an HDF5 group.

        The group is marked so that :func:`hdf_load_dict_from_group` knows to
        rebuild the ordered sequence, and the contents are written in order
        using ITEM* names to avoid the numpy record hack that was previously
        required for pytables.
        """

        target_group.attrs["LIST_NODE"] = True
        if isinstance(seq, tuple):
            target_group.attrs["LISTCLASS"] = "tuple"
        for idx, entry in enumerate(seq):
            item_name = "ITEM" + str(idx)
            if issubclass(type(entry), np.ndarray):
                target_group.create_dataset(item_name, data=entry, dtype=entry.dtype)
            elif issubclass(type(entry), dict):
                nested_dict = target_group.create_group(item_name)
                nested_dict.attrs["DICT_NODE"] = True
                hdf_save_dict_to_group(nested_dict, entry)
            elif isinstance(entry, (list, tuple)):
                nested_group = target_group.create_group(item_name)
                save_sequence(entry, nested_group)
            else:
                if entry is None:
                    continue
                if isinstance(entry, str):
                    target_group.attrs[item_name] = entry.encode("utf-8")
                elif np.isscalar(entry):
                    target_group.attrs[item_name] = entry
                elif hasattr(entry, "__len__") and len(entry) > 0:
                    target_group.attrs[item_name] = [
                        x.encode("utf-8") if isinstance(x, str) else x
                        for x in entry
                    ]
    for k, v in data.items():
        if issubclass(type(v), np.ndarray):
            logger.debug("Dataset type %s" % str(v.dtype))
            logger.debug("Adding %s=%s as dataset" % (k, v))
            group.create_dataset(k, data=v, dtype=v.dtype)
        elif isinstance(v, (list, tuple)):
            subgroup = group.create_group(k)
            save_sequence(v, subgroup)
        elif issubclass(type(v), dict):
            if set(v.keys()) == {"LISTELEMENTS"}:
                elements = v["LISTELEMENTS"]
                if len(elements) > 0 and isinstance(elements[0], str):
                    elements = [x.encode("utf-8") for x in elements]
                arr = np.rec.fromarrays([elements], names="LISTELEMENTS")
                group.attrs[k] = arr
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
        indices = set()
        for name in group.keys():
            if name.startswith("ITEM"):
                indices.add(int(name[4:]))
        for name in group.attrs.keys():
            if name.startswith("ITEM"):
                indices.add(int(name[4:]))
        items = []
        for idx in sorted(indices):
            item_name = "ITEM" + str(idx)
            if item_name in group:
                if isinstance(group[item_name], h5py.Dataset):
                    items.append(group[item_name][()])
                elif "LIST_NODE" in group[item_name].attrs and group[item_name].attrs["LIST_NODE"]:
                    items.append(hdf_load_dict_from_group(group[item_name]))
                else:
                    items.append(hdf_load_dict_from_group(group[item_name]))
            else:
                if isinstance(group.attrs[item_name], bytes):
                    items.append(group.attrs[item_name].decode("utf-8"))
                elif hasattr(group.attrs[item_name], "__len__") and not np.isscalar(group.attrs[item_name]):
                    decoded_attr = []
                    for element in group.attrs[item_name]:
                        if isinstance(element, bytes):
                            decoded_attr.append(element.decode("utf-8"))
                        else:
                            decoded_attr.append(element)
                    items.append(decoded_attr)
                else:
                    if isinstance(group.attrs[item_name], np.generic):
                        items.append(group.attrs[item_name].item())
                    else:
                        items.append(group.attrs[item_name])
        if "LISTCLASS" in group.attrs and (
            group.attrs["LISTCLASS"] == b"tuple"
            or group.attrs["LISTCLASS"] == "tuple"
        ):
            return tuple(items)
        return items
    retval = {}
    for k, v in group.items():
        if isinstance(v, h5py.Dataset):
            retval[k] = v[()]
        else:
            retval[k] = hdf_load_dict_from_group(v)
    for k in group.attrs.keys():
        if (
            hasattr(group.attrs[k], "dtype")
            and getattr(group.attrs[k].dtype, "names", None) is not None
            and "LISTELEMENTS" in group.attrs[k].dtype.names
        ):
            retval[k] = group.attrs[k]
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
