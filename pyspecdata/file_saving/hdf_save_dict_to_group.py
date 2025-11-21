from ..general_functions import *
from pylab import *
import h5py

logger = logging.getLogger("pyspecdata.hdf_save_dict_to_group")


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
    retval = {}
    for k, v in group.items():
        if isinstance(v, h5py.Dataset):
            retval[k] = v[()]
        else:
            retval[k] = hdf_load_dict_from_group(v)
    for k, v in group.attrs.items():
        retval[k] = v
    return retval
