import os
import urllib.request
import re
import datetime
from ..datadir import getDATADIR, pyspec_config
import logging
import requests

__all__ = ["zenodo_download", "zenodo_upload", "create_deposition"]


def _get_token():
    """Return the API access token for Zenodo."""
    token_path = pyspec_config.get_setting("token_file", section="zenodo")
    with open(os.path.expanduser(token_path)) as fp:
        return fp.read().strip()


def zenodo_download(deposition, searchstring, exp_type=None):
    """Download the file from Zenodo ``deposition`` that matches
    ``searchstring`` and place it in the directory associated with
    ``exp_type`` using :func:`getDATADIR`.

    Parameters
    ----------
    deposition : str
        Deposition identifier on Zenodo.
    searchstring : str
        Regular expression used to search the file names inside the deposition.
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

    r = requests.get(f"https://zenodo.org/api/records/{deposition}")
    r.raise_for_status()
    files = r.json().get("files", [])
    pattern = re.compile(searchstring)
    matches = [f for f in files if pattern.search(f.get("key", ""))]
    if len(matches) == 0:
        raise ValueError(
            f"no files matching {searchstring!r} in deposition {deposition}\n"
            f"files: {[f['key'] for f in files]}"
        )
    elif len(matches) > 1:
        raise ValueError(
            f"multiple files match {searchstring!r} in deposition"
            f" {deposition}: {[f['key'] for f in matches]}"
        )
    fileinfo = matches[0]
    dest = os.path.join(dest_dir, fileinfo["key"])
    url = fileinfo["links"]["self"]
    urllib.request.urlretrieve(url, dest)
    logging.debug(f"downloading zenodo '{url}' to '{dest}'")
    return dest


def create_deposition(title):
    """Create a new Zenodo deposition using ``title``.

    The deposition will pre-reserve a DOI, set the upload type to
    ``dataset`` and mark today's date as both the publication date and the
    availability date.
    """

    token = _get_token()

    today = datetime.date.today().isoformat()
    metadata = {
        "title": title,
        # reserve a DOI without providing an existing one
        "prereserve_doi": True,
        "upload_type": "dataset",
        "publication_date": today,
        # ``type`` corresponds to the "Type" field on the website
        "dates": [{"start": today, "end": today, "type": "Available"}],
    }

    r = requests.post(
        "https://zenodo.org/api/deposit/depositions",
        params={"access_token": token},
        json={"metadata": metadata},
    )
    try:
        r.raise_for_status()
    except requests.HTTPError as exc:
        msg = f"{exc}\n{r.text}"
        logging.error("failed to create deposition: %s", msg)
        raise requests.HTTPError(msg) from exc
    return r.json()["id"]


def zenodo_upload(local_path, title=None, deposition_id=None):
    """Upload ``local_path`` to Zenodo.

    Parameters
    ----------
    local_path : str
        Path to the local file that will be uploaded.
    title : str, optional
        Title of the deposition when creating a new one.  ``title`` must be
        provided if ``deposition_id`` is ``None``.
    deposition_id : str, optional
        Existing deposition identifier.  If ``None`` a new deposition is
        created using ``title``.
    """

    token = _get_token()

    if deposition_id is None:
        if title is None:
            raise ValueError(
                "must provide title when creating a new deposition"
            )
        deposition_id = create_deposition(title)

    with open(local_path, "rb") as fp:
        r = requests.post(
            "https://zenodo.org/api/deposit/depositions/"
            f"{deposition_id}/files",
            params={"access_token": token},
            files={"file": fp},
        )

    r.raise_for_status()
    info = r.json()
    print("Uploaded", info["filename"])
    print("View deposition at", f"https://zenodo.org/uploads/{deposition_id}")

    # record the upload in the config file
    n_uploads = int(
        pyspec_config.get_setting(
            "upload_number", section="zenodo", default="0"
        )
    )
    n_uploads += 1
    pyspec_config.set_setting("zenodo", "upload_number", str(n_uploads))
    pyspec_config.set_setting(
        "zenodo", f"upload_deposition{n_uploads}", str(deposition_id)
    )
    pyspec_config.set_setting(
        "zenodo", f"upload_localfile{n_uploads}", local_path
    )
    return deposition_id
