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
import re

# {{{ changeable parameters
# list of (search string, exp_type) pairs for example data files
# in older versions, we auto-added the .*, but that's no longer true
files_to_upload = [
    # Raw ODNP inputs for the 230410/230412 70 mM TEMPOL k_rho workflow in
    # ag_lab_notebook/tempol_series.
    (re.escape("230410_70mM_TEMPOL_ODNP_1.h5"), "ODNP_NMR_comp/ODNP"),
    (re.escape("230410_water_ODNP_2.h5"), "ODNP_NMR_comp/ODNP"),
    (re.escape("230412_70mM_TEMPOL_ODNP_1.h5"), "ODNP_NMR_comp/ODNP"),
    (re.escape("230412_water_ODNP_1.h5"), "ODNP_NMR_comp/ODNP"),
    # proc_scripts master example: FIDtoEcho_Actual.py
    (re.escape("211103_TEMPOL_269uM_HeatExch.h5"), "ODNP_NMR_comp/ODNP"),
    # proc_scripts master example: FIDtoEcho_Actual_Challenging.py
    (
        re.escape("210604_50mM_4AT_AOT_w11_cap_probe_echo.h5"),
        "ODNP_NMR_comp/Echoes",
    ),
    # proc_scripts master example: append_conc_to_h5file.py
    (re.escape("220126_Ras_M67R1a_capProbe.h5"), "ODNP_NMR_comp/ODNP"),
    # proc_scripts master examples: calibrate_tsqrtP.py and
    # verify_tsqrtP_calibration.py
    (
        re.escape("240819_test_amp0p05_calib_pulse_calib.h5"),
        "ODNP_NMR_comp/test_equipment",
    ),
    (
        re.escape("240819_amp0p1_calib_pulse_calib.h5"),
        "ODNP_NMR_comp/test_equipment",
    ),
    (
        re.escape("240819_amp0p2_calib_repeat_pulse_calib.h5"),
        "ODNP_NMR_comp/test_equipment",
    ),
    # proc_scripts master example: compare_CPMG_echo.py
    (re.escape("240704_13p5mM_TEMPOL_CPMG.h5"), "ODNP_NMR_comp/CPMG"),
    (re.escape("240704_13p5mM_TEMPOL_echo.h5"), "ODNP_NMR_comp/Echoes"),
    (re.escape("240702_13p5mM_TEMPOL_CPMG.h5"), "ODNP_NMR_comp/CPMG"),
    (re.escape("240702_13p5mM_TEMPOL_echo.h5"), "ODNP_NMR_comp/Echoes"),
    (re.escape("240702_13p5mM_TEMPOL_pm_CPMG.h5"), "ODNP_NMR_comp/CPMG"),
    (
        re.escape("240702_13p5mM_TEMPOL_pm_echo.h5"),
        "ODNP_NMR_comp/Echoes",
    ),
    (
        re.escape("240702_13p5mM_TEMPOL_pm_30dB_CPMG.h5"),
        "ODNP_NMR_comp/CPMG",
    ),
    (
        re.escape("240702_13p5mM_TEMPOL_pm_30dB_echo.h5"),
        "ODNP_NMR_comp/Echoes",
    ),
    (
        re.escape("240702_13p5mM_TEMPOL_pm_34dB_CPMG.h5"),
        "ODNP_NMR_comp/CPMG",
    ),
    (
        re.escape("240702_13p5mM_TEMPOL_pm_34dB_echo.h5"),
        "ODNP_NMR_comp/Echoes",
    ),
    # proc_scripts master example: generate_SC_PSD.py
    (re.escape("230822_BNC_RX_magon_200kHz.h5"), "ODNP_NMR_comp/Echoes"),
    # proc_scripts master example: generate_oscilloscope_PSD.py
    (
        re.escape("240328_RX_GDS_2mV_analytic.h5"),
        "ODNP_NMR_comp/noise_tests",
    ),
    # proc_scripts master example: integrate_limits_realData.py
    (
        re.escape(
            "210604_50mM_4AT_AOT_w11_cap_probe_echo_tau_11_after_noApod_ex.h5"
        ),
        "ODNP_NMR_comp/processed",
    ),
    # proc_scripts master example: phaseCycNoise_example.py
    (
        re.escape("210607_TEMPOL_100mM_cap_probe_DNP.h5"),
        "ODNP_NMR_comp/ODNP",
    ),
    # proc_scripts master example: proc_CPMG.py
    (
        re.escape("240620_200uM_TEMPOL_pm_generic_CPMG.h5"),
        "ODNP_NMR_comp/Echoes",
    ),
    # proc_scripts master example: proc_FID_nutation.py
    (
        re.escape("240805_amp0p1_27mM_TEMPOL_FID_nutation.h5"),
        "ODNP_NMR_comp/nutation",
    ),
    # proc_scripts master example: proc_GDS_capture_int.py
    (
        re.escape("240819_amp0p05_beta_max_pulse_capture.h5"),
        "ODNP_NMR_comp/test_equipment",
    ),
    (
        re.escape("240819_amp0p1_beta_max_pulse_capture.h5"),
        "ODNP_NMR_comp/test_equipment",
    ),
    (
        re.escape("240819_amp0p2_beta_max_pulse_capture.h5"),
        "ODNP_NMR_comp/test_equipment",
    ),
    # proc_scripts master example: proc_capture_nutation.py
    (re.escape("210204_gds_p90_vary_3.h5"), "nutation"),
    # proc_scripts master example: proc_fieldSweep.py
    (
        re.escape("240924_13p5mM_TEMPOL_field.h5"),
        "ODNP_NMR_comp/field_dependent",
    ),
    # proc_scripts master example: proc_gain.py
    (
        re.escape("240123_power_in_analytic.h5"),
        "ODNP_NMR_comp/noise_tests",
    ),
    (
        re.escape("240123_power_out_analytic.h5"),
        "ODNP_NMR_comp/noise_tests",
    ),
    # proc_scripts master example: proc_nutation.py
    (
        re.escape("240805_amp0p1_27mM_TEMPOL_nutation.h5"),
        "ODNP_NMR_comp/nutation",
    ),
    # proc_scripts master example: proc_rec_response.py
    (
        re.escape("240123_10mV_AFG_GDS_5mV_100MSPS_analytic.h5"),
        "ODNP_NMR_comp/noise_tests",
    ),
    (
        re.escape("240117_afg_sc_10mV_3p9kHz_zoom.h5"),
        "ODNP_NMR_comp/noise_tests",
    ),
    # proc_scripts master example: proc_square_refl.py
    (
        re.escape("210125_sqwv_cap_probe_1.h5"),
        "ODNP_NMR_comp/test_equipment",
    ),
    (
        re.escape("210111_sqwv_sol_probe_1.h5"),
        "ODNP_NMR_comp/test_equipment",
    ),
    # proc_scripts master example: proc_tune_capture.py
    (re.escape("220808_150uM_TEMPOL.npz"), "francklab_esr/alex"),
    (re.escape("220114_3mM_TEMPOL_3b.npz"), "francklab_esr/alex"),
    # proc_scripts master example: QESR.py
    (re.escape("220804_rasI36_MTSL.DSC"), "francklab_esr/Farhana"),
    (re.escape("220804_rasI36_MTSL.DTA"), "francklab_esr/Farhana"),
    (re.escape("220804_water.DSC"), "francklab_esr/Farhana"),
    (re.escape("220804_water.DTA"), "francklab_esr/Farhana"),
    (re.escape("220720_stock_4.DSC"), "francklab_esr/alex"),
    (re.escape("220720_stock_4.DTA"), "francklab_esr/alex"),
    # proc_scripts master example: QESR_samesample.py
    (re.escape("230511_water.DSC"), "francklab_esr/romana"),
    (re.escape("230511_water.DTA"), "francklab_esr/romana"),
    (re.escape("250321_OHT_control.DSC"), "francklab_esr/Warren"),
    (re.escape("250321_OHT_control.DTA"), "francklab_esr/Warren"),
    (re.escape("250321_OHT_gaindown.DSC"), "francklab_esr/Warren"),
    (re.escape("250321_OHT_gaindown.DTA"), "francklab_esr/Warren"),
    (re.escape("250321_OHT_moddown.DSC"), "francklab_esr/Warren"),
    (re.escape("250321_OHT_moddown.DTA"), "francklab_esr/Warren"),
    (
        re.escape("250321_OHT_followmodrule.DSC"),
        "francklab_esr/Warren",
    ),
    (
        re.escape("250321_OHT_followmodrule.DTA"),
        "francklab_esr/Warren",
    ),
    (re.escape("250321_OHT_morepoint.DSC"), "francklab_esr/Warren"),
    (re.escape("250321_OHT_morepoint.DTA"), "francklab_esr/Warren"),
    (re.escape("250321_OHT_lowpower.DSC"), "francklab_esr/Warren"),
    (re.escape("250321_OHT_lowpower.DTA"), "francklab_esr/Warren"),
    # proc_scripts master example: epr_overlay.py
    (re.escape("220307_S175_KCl.DSC"), "francklab_esr/Farhana"),
    (re.escape("220307_S175_KCl.DTA"), "francklab_esr/Farhana"),
    (re.escape("220729_prS175.DSC"), "francklab_esr/Farhana"),
    (re.escape("220729_prS175.DTA"), "francklab_esr/Farhana"),
    (re.escape("220307_S175_KI.DSC"), "francklab_esr/Farhana"),
    (re.escape("220307_S175_KI.DTA"), "francklab_esr/Farhana"),
    (re.escape("220307_prS175_KH2PO4.DSC"), "francklab_esr/Farhana"),
    (re.escape("220307_prS175_KH2PO4.DTA"), "francklab_esr/Farhana"),
    (
        re.escape("240404_L56_MTSL_Rasbatch240320_fraction3.DSC"),
        "francklab_esr/warren",
    ),
    (
        re.escape("240404_L56_MTSL_Rasbatch240320_fraction3.DTA"),
        "francklab_esr/warren",
    ),
    (
        re.escape("240404_L56_MTSL_Rasbatch240320_fraction4.DSC"),
        "francklab_esr/warren",
    ),
    (
        re.escape("240404_L56_MTSL_Rasbatch240320_fraction4.DTA"),
        "francklab_esr/warren",
    ),
    (
        re.escape("240404_L56_MTSL_Rasbatch240320_fraction5.DSC"),
        "francklab_esr/warren",
    ),
    (
        re.escape("240404_L56_MTSL_Rasbatch240320_fraction5.DTA"),
        "francklab_esr/warren",
    ),
    # proc_scripts master example: epr_overlay2.py
    (
        re.escape("250130_20mM_TSO4_W20_Isooct_run1.DSC"),
        "francklab_esr/romana",
    ),
    (
        re.escape("250130_20mM_TSO4_W20_Isooct_run1.DTA"),
        "francklab_esr/romana",
    ),
    (
        re.escape("240220_TempoSO4_W20_run1.DSC"),
        "francklab_esr/romana",
    ),
    (
        re.escape("240220_TempoSO4_W20_run1.DTA"),
        "francklab_esr/romana",
    ),
    (
        re.escape("250206_20mM_TSO4_W20_C8.DSC"),
        "francklab_esr/romana",
    ),
    (
        re.escape("250206_20mM_TSO4_W20_C8.DTA"),
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
            title="Documentation example data files for pySpecdata",
        )
    else:
        zenodo_upload(local_path, deposition_id=deposition_id)

print("View deposition at https://zenodo.org/uploads/" + str(deposition_id))
