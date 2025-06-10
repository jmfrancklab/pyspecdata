"""Upload a file from the examples directory to Zenodo
------------------------------------------------------

This example demonstrates how to locate an example data file, ensure it is
present locally, and then upload it to Zenodo.  The file we upload is the same
one used in :mod:`examples.UV.Cary_simple`.

To run this script you must create a personal access token on the Zenodo
website (with ``deposit:write`` scope).  Save the token in a file and reference
it from the ``[zenodo]`` section of ``~/.pyspecdata``::

    [zenodo]
    token_file = /path/to/zenodo.token
    deposition_id = 123456

``deposition_id`` is the numeric identifier of an existing deposition record.
When you create a new upload on Zenodo, the URL ends with this ID.  Files
uploaded through the API are added to that deposition.
"""

import os
import requests
from pyspecdata import find_file, search_filename
from pyspecdata.datadir import pyspec_config

# locate the data using the same search string as :mod:`examples.UV.Cary_simple`
find_file("T177R1a_pR_210615", exp_type="UV_Vis/proteorhodopsin")
local_path = search_filename(
    "T177R1a_pR_210615",
    exp_type="UV_Vis/proteorhodopsin",
    unique=True,
)

# retrieve authentication info from the config file
token_path = pyspec_config.get_setting("token_file", section="zenodo")
deposition_id = pyspec_config.get_setting("deposition_id", section="zenodo")

with open(os.path.expanduser(token_path)) as fp:
    token = fp.read().strip()

with open(local_path, "rb") as fp:
    r = requests.post(
        f"https://zenodo.org/api/deposit/depositions/{deposition_id}/files",
        params={"access_token": token},
        files={"file": fp},
    )

r.raise_for_status()
info = r.json()
print("Uploaded", info["filename"], "download URL:", info["links"]["download"]) 
