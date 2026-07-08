"""Upload non-UV documentation example files to Zenodo
-----------------------------------------------------

This example demonstrates how to locate the non-UV data files used in the
documented :mod:`examples` package and upload them all to the same Zenodo
deposition.  UV-Vis example data are already covered by the Zenodo download
example.

To run this script you must create a personal access token on the Zenodo
website (with ``deposit:write`` scope).  Save the token in a file and reference
it from ``~/.pyspecdata``::

    [zenodo]
    token_file = /path/to/zenodo.token

A new deposition record will be created automatically for the first file and
the remaining files will be uploaded to that same deposition.  Keep token
files local and never commit them to git.  See
https://developers.zenodo.org/#rest-api for API and authentication details.
"""

from pyspecdata import search_filename, zenodo_upload
import re

# {{{ changeable parameters
# list of (search string, exp_type) pairs for example data files
# in older versions, we auto-added the .*, but that's no longer true
files_to_upload = [
    # examples/ESR/QESR.py
    (
        r"^" + re.escape("QESR_Test_WaterCap_Background_210923.DSC") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_Test_WaterCap_Background_210923.DTA") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_150uM_TEMPOL_1_noglyc_210923.DSC") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_150uM_TEMPOL_1_noglyc_210923.DTA") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_150uM_TEMPOL_2_noglyc_210923.DSC") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_150uM_TEMPOL_2_noglyc_210923.DTA") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_150uM_TEMPOL_3_noglyc_210923.DSC") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_150uM_TEMPOL_3_noglyc_210923.DTA") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_150uM_TEMPOL_4_wglyc_210923.DSC") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_150uM_TEMPOL_4_wglyc_210923.DTA") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_150uM_TEMPOL_5_wglyc_210923.DSC") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_150uM_TEMPOL_5_wglyc_210923.DTA") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_150uM_TEMPOL_6_wglyc_210923.DSC") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("QESR_150uM_TEMPOL_6_wglyc_210923.DTA") + r"$",
        "francklab_esr/Sam",
    ),
    # examples/ESR/esr_example.py
    (
        r"^" + re.escape("15N_S175R1a_pR_DHPC_today_200304.DSC") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^" + re.escape("15N_S175R1a_pR_DHPC_today_200304.DTA") + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^"
        + re.escape("15N_S175R1a_pR_DHPC_Power_Sat_6G_200303.DSC")
        + r"$",
        "francklab_esr/Sam",
    ),
    (
        r"^"
        + re.escape("15N_S175R1a_pR_DHPC_Power_Sat_6G_200303.DTA")
        + r"$",
        "francklab_esr/Sam",
    ),
    # examples/ESR/epr_u_domain.py
    (
        r"^" + re.escape("220307_S175_KCl.DSC") + r"$",
        "francklab_esr/Farhana",
    ),
    (
        r"^" + re.escape("220307_S175_KCl.DTA") + r"$",
        "francklab_esr/Farhana",
    ),
    # examples/calculate_covariance.py
    (
        r"^" + re.escape("230504_3p8mM_TEMPOL_stb_wt_4x.DSC") + r"$",
        "francklab_esr/alex",
    ),
    (
        r"^" + re.escape("230504_3p8mM_TEMPOL_stb_wt_4x.DTA") + r"$",
        "francklab_esr/alex",
    ),
    # examples/temp.py
    (
        r"^" + re.escape("250123_TEMPOL_100uM_AG_Covariance_2D.DSC") + r"$",
        "francklab_esr/romana",
    ),
    (
        r"^" + re.escape("250123_TEMPOL_100uM_AG_Covariance_2D.DTA") + r"$",
        "francklab_esr/romana",
    ),
    (
        r"^"
        + re.escape("250123_TEMPOL_100uM_AG_Covariance_2D_cc12.DSC")
        + r"$",
        "francklab_esr/romana",
    ),
    (
        r"^"
        + re.escape("250123_TEMPOL_100uM_AG_Covariance_2D_cc12.DTA")
        + r"$",
        "francklab_esr/romana",
    ),
]
# }}}

deposition_id = None
for search_str, exp_type in files_to_upload:
    local_path = search_filename(search_str, exp_type=exp_type, unique=True)
    if deposition_id is None:
        deposition_id = zenodo_upload(
            local_path,
            title="Documentation example data files for pySpecData",
        )
    else:
        zenodo_upload(local_path, deposition_id=deposition_id)

print("View deposition at https://zenodo.org/uploads/" + str(deposition_id))
