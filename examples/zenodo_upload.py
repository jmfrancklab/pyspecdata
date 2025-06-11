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

The script will create a new deposition record automatically and then upload
the file to that deposition.
"""

import os
from pyspecdata import search_filename, zenodo_upload

# locate the data using the same search string as
# :mod:`examples.UV.Cary_simple`
local_path = search_filename(
    "T177R1a_pR_210615",
    exp_type="UV_Vis/proteorhodopsin",
    unique=True,
)

zenodo_upload(local_path)
