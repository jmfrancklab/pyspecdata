"""Upload the UV example files to Zenodo
-------------------------------------

This example demonstrates how to locate each data file used in the
:mod:`examples.UV` package and upload them all to the same Zenodo deposition.
The deposition reserves a DOI, sets the resource type to ``dataset``, and
marks today's date as available.

To run this script you must create a personal access token on the Zenodo
website (with ``deposit:write`` scope).  Save the token in a file and reference
it from the ``[zenodo]`` section of ``~/.pyspecdata``::

    [zenodo]
    token_file = /path/to/zenodo.token

A new deposition record will be created automatically for the first file and
the remaining files will be uploaded to that same deposition.
"""

from pyspecdata import search_filename, zenodo_upload

# list of (search string, exp_type) pairs for all UV examples
files_to_upload = [
    ("T177R1a_pR_210615", "UV_Vis/proteorhodopsin"),
    (
        "221110_BSAexerciseWK_0p07-0percentBSAcalibration.BSW",
        "UV_Vis/BSA_Exercise",
    ),
    ("200703_Ellman_before_SL.DSW", "UV_Vis/Ellmans_Assay"),
    ("Ras_Stability4", "UV_Vis/Ras_stability/200803_RT"),
]

deposition_id = None
for search_str, exp_type in files_to_upload:
    local_path = search_filename(search_str, exp_type=exp_type, unique=True)
    if deposition_id is None:
        deposition_id = zenodo_upload(
            local_path,
            title="UV-Vis documentation examples for pySpecdata",
        )
    else:
        zenodo_upload(local_path, deposition_id=deposition_id)

print("View deposition at https://zenodo.org/uploads/" + str(deposition_id))
