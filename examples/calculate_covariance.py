"""
Calculation of the Covariance Matrix
====================================
After rescaling plots, the covariance matrix is calculated
and then plotted for a 2D Field experiment
"""
from pyspecProcScripts import *
from pyspecdata import *
from pylab import *

fieldaxis = "$B_0$"
exp_type = "francklab_esr/alex"
with figlist_var() as fl:
    for filenum, (thisfile, calibration, diameter) in enumerate(
        [("230504_3p8mM_TEMPOL_stb_wt_4x.DSC", "230202", "QESR caps")]
    ):
        d = find_file(thisfile, exp_type=exp_type)["harmonic", 0]
        d.rename("'Time'", "observations")
        d.reorder(["observations", fieldaxis])
        d /= QESR_scalefactor(d, calibration_name=calibration, diameter_name=diameter)
        d.ift(fieldaxis, shift=True)
        fl.next("Covariance in U domain")
        fl.image(d.C.cov_mat("observations"))
        d.ft(fieldaxis)
        fl.next("covariance in B domain")
        fl.image(d.C.cov_mat("observations"))
