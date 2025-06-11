import os
import urllib.request
from ..datadir import getDATADIR, pyspec_config
import logging
import requests

__all__ = ["zenodo_download", "zenodo_upload"]


def zenodo_download(url, exp_type=None):
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


def zenodo_upload(local_path):
    # retrieve authentication info from the config file
    token_path = pyspec_config.get_setting("token_file", section="zenodo")

    with open(os.path.expanduser(token_path)) as fp:
        token = fp.read().strip()

    # create an empty deposition so we can upload files
    r = requests.post(
        "https://zenodo.org/api/deposit/depositions",
        params={"access_token": token},
        json={},
    )
    r.raise_for_status()
    deposition = r.json()

    with open(local_path, "rb") as fp:
        r = requests.post(
            "https://zenodo.org/api/deposit/depositions/"
            f"{deposition['id']}/files",
            params={"access_token": token},
            files={"file": fp},
        )

    r.raise_for_status()
    info = r.json()
    print(
        "Uploaded",
        info["filename"],
        "download URL:",
        info["links"]["download"],
    )
    print("View deposition at", deposition["links"]["html"])
