""" Calculate covariance matrix for NMR signal
==============================================
An example that shows the covariance matrix using both numpy and pyspecdata
operations. The covariance of the FID tail shows the expected diagonal.
"""

import pyspecdata as psd
import numpy as np

signal_pathway = {"ph1": 1}
thisdata = psd.find_file(
    "241003_27mM_TEMPOL_amp0p1_var_tau_pm_echo",
    exp_type="ODNP_NMR_comp/Echoes",
    expno="echo_1",
    postproc="none",
).squeeze()
with psd.figlist_var() as fl:
    thisdata = thisdata["ph1", 0]  # pull only one phcyc step to simplify
    # I explictly reorder so I know which dimension is which when I do raw
    # numpy and matplotlib
    thisdata.reorder(["t2", "nScans"])
    mycov = abs(np.cov(thisdata.data))
    fl.next("covariance")
    fl.image(mycov)
    fl.next("covariance - pyspecdata")
    fl.image(abs(thisdata.C.cov_mat("nScans")))
    # now I zoom in to show that when I'm in the part of the FID that's
    # just noise (no shifting frequencies, I get the diagonal covariance
    # that I expect)
    fl.next("covariance - FID tail")
    fl.image(abs(thisdata["t2":(0.6,0.7)].cov_mat("nScans")))
