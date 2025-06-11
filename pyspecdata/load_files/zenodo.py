import os
import urllib.request
from ..datadir import getDATADIR
import logging

__all__ = ["download"]


def download(url, exp_type=None):
    """Download a file from Zenodo given by ``url`` and place it in the
    directory associated with ``exp_type`` using :func:`getDATADIR`.

    Parameters
    ----------
    url : str
        Direct URL to the file on zenodo (``?download=1`` link).
    exp_type : str
        Experiment type used to determine where the file should be stored via
        :func:`getDATADIR`.
    Returns
    -------
    str
        Path to the downloaded file.
    """
    if exp_type is None:
        raise ValueError("must provide exp_type")
    dest_dir = getDATADIR(exp_type=exp_type)
    os.makedirs(dest_dir, exist_ok=True)
    fname = os.path.basename(url.split("?", 1)[0])
    dest = os.path.join(dest_dir, fname)
    urllib.request.urlretrieve(url, dest)
    logging.debug(f"downloading zenodo '{url}' to '{dest}'")
    return dest
