"""Transfer files between pyspecdata and Zenodo
------------------------------------------------

This module provides :func:`zenodo_download` for downloading files from
published records or authenticated drafts, :func:`create_deposition` for
creating a draft, and :func:`zenodo_upload` for uploading a file to a draft.
"""

import os
import urllib.request
import urllib.parse
import re
import datetime
import argparse
import csv
import shutil
from ..datadir import getDATADIR, pyspec_config
import logging
import requests

__all__ = ["zenodo_download", "zenodo_upload", "create_deposition", "cmd"]

ZENODO_ID_RE = re.compile(r"^[1-9][0-9]{5,9}$")


def _raise_for_status(response, action):
    """Raise an HTTP error that includes Zenodo's response body."""
    try:
        response.raise_for_status()
    except requests.HTTPError as exc:
        msg = f"{action}: {exc}"
        response_text = getattr(response, "text", "")
        if response_text != "":
            msg += f"\n{response_text}"
        logging.error(msg)
        raise requests.HTTPError(msg, response=response) from exc


def _auth_headers():
    """Return authorization headers for authenticated Zenodo API requests.

    Zenodo accepts bearer-token authentication for the deposition API.
    Keeping the token in the Authorization header avoids putting it in
    request URLs and lets the same helper work for GET, POST, and PUT
    calls.
    """
    token_path = pyspec_config.get_setting("token_file", section="zenodo")
    with open(os.path.expanduser(token_path)) as fp:
        token = fp.read().strip()
    return {"Authorization": f"Bearer {token}"}


def zenodo_download(deposition, searchstring, exp_type=None):
    """Download the file from Zenodo ``deposition`` that matches
    ``searchstring`` and place it in the directory associated with
    ``exp_type`` using :func:`getDATADIR`.

    The public records API is tried first.  If Zenodo reports that the public
    record is not found, this function then tries the authenticated
    draft-deposition API using the configured Zenodo token.

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

    try:
        r = requests.get(f"https://zenodo.org/api/records/{deposition}")
        r.raise_for_status()
    except requests.HTTPError:
        # {{{ try for a draft deposition, since we can't find a private one
        if getattr(r, "status_code", None) != 404:
            raise
        logging.debug(
            "Zenodo public record was not found.  This may be an unpublished "
            "draft deposition; trying the authenticated draft API."
        )
        r = requests.get(
            f"https://zenodo.org/api/deposit/depositions/{deposition}/files",
            headers=_auth_headers(),
        )
        try:
            r.raise_for_status()
        except requests.HTTPError as exc:
            status_code = r.status_code
            if status_code == 401:
                detail = (
                    "draft deposition access requires a valid Zenodo token"
                )
            elif status_code == 403:
                detail = (
                    "the Zenodo token does not have permission for this "
                    "deposition"
                )
            elif status_code == 404:
                detail = (
                    "the deposition may not exist, may not belong to this "
                    "token, may not be a draft, or may be on the wrong "
                    "Zenodo host"
                )
            else:
                detail = "Unknown HTTPError"
            raise requests.HTTPError(f"{detail}\n{r.text}") from exc
        # The draft file-list endpoint returns the list of file objects
        # directly, rather than wrapping it in a top-level "files" field.
        files = r.json()
        # Draft deposition file metadata is returned directly as a list.
        # In that API, the filename is stored under "filename" and the
        # authenticated download URL is stored under links["download"].
        name_key = "filename"
        url_key = "download"
        deposition_label = "draft deposition"
        draft = True
        # }}}
    else:
        # {{{ Published-record metadata is a dictionary with a top-level
        #     "files" list.
        #     In that API, the filename is stored under "key" and the public
        #     file URL is stored under links["self"].
        files = r.json().get("files", [])
        name_key = "key"
        url_key = "self"
        deposition_label = "deposition"
        draft = False
        # }}}

    logging.debug("all the files are " + str(files))
    pattern = re.compile(searchstring)
    matches = [f for f in files if pattern.search(f.get(name_key, ""))]
    if len(matches) == 0:
        raise ValueError(
            f"no files matching {searchstring!r} in {deposition_label} "
            f"{deposition}\nfiles: {[f.get(name_key) for f in files]}"
        )
    elif len(matches) > 1:
        raise ValueError(
            f"multiple files match {searchstring!r} in {deposition_label}"
            f" {deposition}: {[f.get(name_key) for f in matches]}"
        )
    fileinfo = matches[0]
    filename = fileinfo.get(name_key)
    url = fileinfo.get("links", {}).get(url_key)
    if url is None:
        raise ValueError(
            f"no usable download link for {filename!r} in "
            f"{deposition_label} {deposition}"
        )
    dest = os.path.join(dest_dir, filename)
    # urlretrieve cannot attach the Authorization header needed by draft
    # content links, so use urlopen for both public and draft downloads.
    download_url = url
    if draft:
        download_url = urllib.request.Request(url, headers=_auth_headers())
    with (
        urllib.request.urlopen(download_url) as response,
        open(dest, "wb") as fp,
    ):
        # Stream the response so large Zenodo files are not held fully in
        # memory.
        shutil.copyfileobj(response, fp)
    logging.debug(f"downloading zenodo '{url}' to '{dest}'")
    print(f"Downloaded from zenodo '{url}' to {dest}")
    return dest


# TODO ☐: the following is fine -- if it doesn't pass the single-use test, then ask codex to edit the single use linter so that it respects the following keyword:
# SINGLE_USE_EXCEPTION (exported by __all__, even though it's also used above)
def create_deposition(title):
    """Create a new Zenodo deposition using ``title``.

    The deposition will pre-reserve a DOI, set the upload type to
    ``dataset`` and mark today's date as both the publication date and the
    availability date.
    """

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
        # Older examples often pass "access_token" as a query parameter.
        # Use the Authorization header instead so the token is not embedded
        # in a URL that can be copied into logs, shell history, or tracebacks.
        headers=_auth_headers(),
        json={"metadata": metadata},
    )
    _raise_for_status(r, "failed to create Zenodo deposition")
    return r.json()["id"]


def zenodo_upload(local_path, title=None, deposition_id=None):
    """Upload ``local_path`` to Zenodo.

    A new deposition record will be created automatically for the first file and
    the remaining files will be uploaded to that same deposition.  Keep token
    files local and never commit them to git.

    To share a draft (unpublished record), click ``Share`` on the right side
    and add your collaborator, who must also have a Zenodo account.
    Then click the ``Preview`` button, also on the right side.
    Either the panel on the right changes, indicating success, or you get a
    small red error message at the very top of the page, typically because
    author information or other required metadata is missing.

    To use this function, you must create a personal access token on the Zenodo
    website:

    Getting zenodo authorization set up
    ```````````````````````````````````

    Go to your Zenodo profile, then ``Applications``, create the
    token there, and select the ``deposit:write`` scope.
    Zenodo shows the token only once, so copy it before leaving the page.

    Create a local token file, for example ``~/.zenodo_token``,
    and write the token into that file.
    Then use the ``pyspecdata_dataconfig`` command-line
    tool (which opens a GUI as long as you have pyside installed)
    to add the token file to ``~/.pyspecdata``
    If you don't want to use the GUI, you can edit ``~/.pyspecdata``
    manually under ``[zenodo]``::

        [zenodo]
        token_file = ~/.zenodo_token


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

    if deposition_id is None:
        if title is None:
            raise ValueError(
                "must provide title when creating a new deposition"
            )
        deposition_id = create_deposition(title)

    headers = _auth_headers()
    r = requests.get(
        f"https://zenodo.org/api/deposit/depositions/{deposition_id}",
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


def cmd(argv=None):
    """Upload files from ``./data_files.csv`` to a Zenodo draft.

    The command line interface is ``pyspecdata_zenodo <draft-id-or-title>``.

    -  If the single argument matches ``^[1-9][0-9]{5,9}$``, it is interpreted
       as an existing Zenodo draft deposition id.
    -  Otherwise, it is used verbatim as the title for a new draft deposition.

    Files are read from the ``data_files.csv`` created in the current directory
    by :func:`pyspecdata.find_file` and :func:`pyspecdata.search_filename`.
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
    # In the future log_path could be set to a cmdline arg.
    log_path = "data_files.csv"
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
