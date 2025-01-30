""" Calculate covariance matrix for NMR signal
==============================================
First, we show (both with pyspec and then with raw nddata) that the early
points of the FID change *together* due to the fact that the peak is moving
around in frequency space. 
This means that the phase of the time domain points is not fixed -- as the
phase varies, the real and imaginary components of the signal *change together*
- covariance.
However, we show that with the application of a correlation alignment run prior
to calculating the covariance, this effect and therefore covariance is
decreased.
"""

import pyspecdata as psd
import numpy as np
import pyspecProcScripts as psdpr
import matplotlib.pyplot as plt

signal_pathway = {"ph1": 1}
d = psd.find_file(
    "241003_27mM_TEMPOL_amp0p1_var_tau_pm_echo",
    exp_type="ODNP_NMR_comp/Echoes",
    expno="echo_1",
    postproc="none",
).squeeze()
# {{{ Make copy of data and align to compare with unaligned
#     data
aligned = d.C
aligned.ft("t2", shift=True)
opt_shift, sigma, mask_fn = psdpr.correl_align(
    aligned,
    indirect_dim="nScans",
    signal_pathway=signal_pathway,
)
aligned.ift("t2")
aligned *= np.exp(-1j * 2 * np.pi * opt_shift * aligned.fromaxis("t2"))
# }}}
with psd.figlist_var() as fl:
    for thisdata, label in [(d, "Unaligned"), (aligned, "With Alignment")]:
        thisdata = thisdata["ph1", 0]  # pull only one phcyc step to simplify
        # I explictly reorder so I know which dimension is which when I do raw
        # numpy and matplotlib
        thisdata.reorder(["t2", "nScans"])
        fl.next("pyspec sub - %s" % label)
        fl.image(thisdata - thisdata.C.mean("nScans"))
        # I do still use raw numpy/matplotlib so that I can check what numpy
        # does
        fl.next("raw data, sub, numpy - %s" % label)
        fl.image(
            (
                thisdata.data - np.mean(thisdata.data, axis=1)[:, np.newaxis]
            ).real
        )
        fl.next("pyspec ft - %s" % label)
        if label == "aligned":
            fl.image(thisdata.C.ft("t2")["t2":(-200, 0)].run(abs))
        else:
            fl.image(
                thisdata.C.ft("t2", shift=True)["t2":(-100, 100)].run(abs)
            )
        fig, ax_list = plt.subplots(2, 2)
        fig.suptitle(label)
        fl.next("%s" % label, fig=fig)
        mycov = abs(np.cov(thisdata.data))
        fl.image(mycov, ax=ax_list[0][0])
        ax_list[0][0].set_title("covariance")
        fl.image(abs(thisdata.cov_mat("nScans")), ax=ax_list[0][1])
        ax_list[0][1].set_title("covariance, pyspec")
        # now I zoom in to show that when I'm in the part of the FID that's
        # just noise (no shifting frequencies, I get the diagonal covariance
        # that I expect)
        fl.image(
            abs(thisdata["$t2_i$":(0.6, 0.7)]["$t2_j$":(0.6, 0.7)]),
            ax=ax_list[1][0],
        )
        ax_list[1][0].set_title("covariance, zoomed late times")
        fl.image(
            abs(thisdata["$t2_i$":(0.9, 1.0)]["$t2_j$":(0.9, 1.0)]),
            ax=ax_list[1][1],
        )
        ax_list[1][1].set_title("covariance, zoomed even later")
