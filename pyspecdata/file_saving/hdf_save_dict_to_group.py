import logging
import numpy as np

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
                arr = np.rec.fromarrays(
                    [v["LISTELEMENTS"]], names="LISTELEMENTS"
                )
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
                logger.debug("Adding %s=%s as scalar attribute" % (k, v))
                group.attrs[k] = v
            elif hasattr(v, "__len__") and len(v) > 0:
                logger.debug("Adding %s=%s as list attribute" % (k, v))
                # added encoding in following
                group.attrs[k] = [
                    x.encode("utf-8") if isinstance(x, str) else x for x in v
                ]
