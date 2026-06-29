"""Download files from Zenodo using :func:`search_filename`
----------------------------------------------------------

This example shows how to retrieve a file from Zenodo when it is not already
present locally.  The download occurs transparently via the ``zenodo``
keyword of :func:`pyspecdata.search_filename`.

This downloads the same files as the upload example.
"""

import os
import re

from pyspecdata import search_filename

# {{{ changeable parameters
REMOVE_LOCAL_COPIES = False

files_to_download = {
    "21041480": [
        (".*T177R1a_pR_210615.*", "UV_Vis/proteorhodopsin"),
        (
            re.escape("221110_BSAexerciseWK_0p07-0percentBSAcalibration.BSW"),
            "UV_Vis/BSA_Exercise",
        ),
        (re.escape("200703_Ellman_before_SL.DSW"), "UV_Vis/Ellmans_Assay"),
        (".*Ras_Stability4.*", "UV_Vis/Ras_stability/200803_RT"),
    ],
}
# }}}

for deposition, files in files_to_download.items():
    for search_str, exp_type in files:
        if REMOVE_LOCAL_COPIES:
            try:
                path = search_filename(
                    search_str, exp_type=exp_type, unique=True
                )
            except RuntimeError:
                pass
            else:
                os.unlink(path)
                print(f"Removed local copy {path}")
        path = search_filename(
            search_str,
            exp_type=exp_type,
            unique=True,
            zenodo=deposition,
        )
        print("Found at", path)
