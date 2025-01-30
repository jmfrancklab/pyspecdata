import pyspecdata as psd
import numpy as np

d = psd.find_file(
    "241003_27mM_TEMPOL_amp0p1_var_tau_pm_echo",
    exp_type="ODNP_NMR_comp/Echoes",
    expno="echo_1",
    postproc="none",
)
d = d["ph1", 0].squeeze()  # pull only one phcyc step to simplify
# I explictly reorder so I know which dimension is which when I do raw numpy
# and matplotlib
d.reorder(["t2", "nScans"])
# I do still use raw numpy/matplotlib so that I can check what numpy does
with psd.figlist_var() as fl:
    fl.next("raw data")
    fl.image(d.data.real)
    # next, first with pyspec and then with raw nddata, I show that the early
    # points of the fid change *together*, which means that we expect the
    # non-diagonal covariance
    fl.next("pyspec sub")
    fl.image(d - d.C.mean("nScans"))
    fl.next("raw data, sub, numpy")
    fl.image((d.data - np.mean(d.data, axis=1)[:, np.newaxis]).real)
    # Why are the points of the FID changing together? Is the peak changing in
    # frequency, or changing in amplitude?
    fl.next("pyspec ft")
    fl.image(d.C.ft("t2", shift=True)["t2":(-100, 100)].run(abs))
    fl.next("covariance")
    mycov = abs(np.cov(d.data))
    fl.image(mycov)
    # this covariance is not diagonal! → but the previous lines show that
    # it's from the fact that the peak is moving around in frequency
    # space, so the phase of the time domain points is not fixed -- as the
    # phase varies, the real and the imaginary components of the signal
    # *change together*, which is covariance.  In order to see the *true*
    # covariance of the signal, we're going to need to run a correlation
    # alignment first!
    # confirm that the following looks the same
    fl.next("covariance, pyspec")
    fl.image(abs(d.cov_mat("nScans")))
    # now I zoom in to show that when I'm in the part of the FID that's
    # just noise (no shifting frequencies, I get the diagonal covariance
    # that I expect)
    aligned = d.C
    fl.next("covariance, zoomed late times")
    fl.image(abs(d["$t2_i$":(0.6, 0.7)]["$t2_j$":(0.6, 0.7)]))
    fl.next("covariance, zoomed even later")
    fl.image(abs(d["$t2_i$":(0.9, 1.0)]["$t2_j$":(0.9, 1.0)]))
