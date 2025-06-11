"""Download a file from Zenodo using :func:`search_filename`
-----------------------------------------------------------

This example shows how to retrieve a file from Zenodo when it is not already
present locally.  The download occurs transparently via the ``zenodo``
keyword of :func:`pyspecdata.search_filename`.

This downloads the same file as the upload example.
"""

from pyspecdata import search_filename

path = search_filename(
    "Pure_T177R1a_pR_210615.BSW",
    exp_type="UV_Vis/proteorhodopsin",
    unique=True,
    zenodo=(
        "https://zenodo.org/api/records/15636512/draft/files/"
        "Pure_T177R1a_pR_210615.BSW/content"
    ),
)

print(f"Downloaded to {path}")
