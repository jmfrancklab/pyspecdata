import pyspecdata as psd
import pyspecProcScripts as psdpr

data = psd.find_file(
    "241003_27mM_TEMPOL_amp0p1_var_tau_pm_echo",
    exp_type="ODNP_NMR_comp/Echoes",
    expno="echo_1",
    lookup=psdpr.lookup_table,
)
data.ift("t2")
# FID slice
data = data["t2":(0, None)]
data *= 2
data["t2":0] *= 0.5
# }}}
# relabel nScans axis with indices
data.setaxis("nScans", "#")
with psd.figlist_var() as fl:
    write
    covariance_matrix = psdpr.select_pathway(
        data, data.get_prop("coherence_pathway")
    ).cov_mat("nScans")
    fl.next("2D Covariance Matrix")
    fl.image(abs(covariance_matrix))
