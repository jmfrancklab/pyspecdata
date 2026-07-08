import os
import urllib.request
import urllib.parse
import re
import datetime
import argparse
import csv
from ..datadir import getDATADIR, pyspec_config
import logging
import requests

__all__ = ["zenodo_download", "zenodo_upload", "create_deposition", "cmd"]

ZENODO_ID_RE = re.compile(r"^[1-9][0-9]{5,9}$")


def _raise_for_status(response, action):
    try:
        response.raise_for_status()
    except requests.HTTPError as exc:
        msg = f"{action}: {exc}"
        response_text = getattr(response, "text", "")
        if response_text != "":
            msg += f"\n{response_text}"
        logging.error(msg)
        raise requests.HTTPError(msg, response=response) from exc


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
    logging.debug("all the files are "+str(files))
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
    print(f"Downloaded from zenodo '{url}' to {dest}")
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
    _raise_for_status(r, "failed to create Zenodo deposition")
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

    headers = {"Authorization": f"Bearer {token}"}
    r = requests.get(
        "https://zenodo.org/api/deposit/depositions/" f"{deposition_id}",
        headers=headers,
    )
    _raise_for_status(r, f"failed to fetch Zenodo deposition {deposition_id}")
    bucket_url = r.json().get("links", {}).get("bucket")
    if bucket_url is None:
        raise ValueError(
            f"Zenodo deposition {deposition_id} did not include a bucket "
            "upload URL; make sure this is an editable draft deposition."
        )

    filename = os.path.basename(local_path)
    quoted_filename = urllib.parse.quote(filename, safe="")
    upload_url = bucket_url.rstrip("/") + "/" + quoted_filename
    with open(local_path, "rb") as fp:
        r = requests.put(upload_url, data=fp, headers=headers)
    _raise_for_status(r, f"failed to upload {filename!r} to Zenodo")
    info = r.json()
    print(
        "Uploaded",
        info.get("key", info.get("filename", info.get("name", filename))),
    )
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


def data_files_from_log(log_path="data_files.csv"):
    """Return local paths listed in the current directory's data file log.

    The log is the ``data_files.csv`` generated by
    :func:`pyspecdata.find_file` and :func:`pyspecdata.search_filename`.
    """
    required_fields = {"Filename", "Path", "exp_type"}
    if not os.path.exists(log_path):
        raise FileNotFoundError(
            "pyspecdata_zenodo must be run from a directory containing "
            "data_files.csv"
        )
    with open(log_path, "r", encoding="utf-8", newline="") as fp:
        reader = csv.DictReader(fp)
        if reader.fieldnames is None:
            raise ValueError(f"{log_path} is empty")
        missing = required_fields - set(reader.fieldnames)
        if len(missing) > 0:
            raise ValueError(
                f"{log_path} is missing required column(s): "
                + ", ".join(sorted(missing))
            )
        local_paths = [
            os.path.join(row["Path"], row["Filename"])
            for row in reader
            if row["Filename"] != ""
        ]
    if len(local_paths) == 0:
        raise ValueError(f"{log_path} does not list any files")
    for this_path in local_paths:
        if not os.path.exists(this_path):
            raise FileNotFoundError(
                f"{this_path} from {log_path} does not exist"
            )
    return local_paths


def cmd(argv=None):
    """Upload files from ``./data_files.csv`` to a Zenodo draft.

    The command line interface is ``pyspecdata_zenodo <draft-id-or-title>``.
    If the single argument matches ``^[1-9][0-9]{5,9}$``, it is interpreted as
    an existing Zenodo draft deposition id.  Otherwise, it is used verbatim as
    the title for a new draft deposition.  Files are read from the
    ``data_files.csv`` created in the current directory by
    :func:`pyspecdata.find_file` and :func:`pyspecdata.search_filename`.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Upload files listed in ./data_files.csv to a Zenodo draft "
            "deposition."
        )
    )
    parser.add_argument(
        "draft_id_or_title",
        help=(
            "Existing Zenodo draft id, or the title to use when creating a "
            "new draft."
        ),
    )
    args = parser.parse_args(argv)
    local_paths = data_files_from_log()
    if ZENODO_ID_RE.match(args.draft_id_or_title):
        deposition_id = args.draft_id_or_title
        title = None
    else:
        deposition_id = None
        title = args.draft_id_or_title
    for j, local_path in enumerate(local_paths):
        if j == 0 and deposition_id is None:
            deposition_id = zenodo_upload(local_path, title=title)
        else:
            deposition_id = zenodo_upload(
                local_path, deposition_id=deposition_id
            )
    print("View deposition at", f"https://zenodo.org/uploads/{deposition_id}")
    return 0
