from pylab import *
from pyspecdata import *
from pyspecdata import DCCT
from pyspecProcScripts import fake_data
from numpy.random import normal, seed
from numpy.linalg import norm
import sympy as s
from collections import OrderedDict

seed(2021)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 4
init_logging(level="debug")

with figlist_var() as fl:
    # {{{ generate the fake data
    # this generates fake clean_data w/ a T1 of 0.2s
    # amplitude of 21, just to pick a random amplitude
    # offset of 300 Hz, FWHM 10 Hz

    # provide the symbols that we use for the fake data:
    t2, td, vd, ph1, ph2 = s.symbols("t2 td vd ph1 ph2")
    echo_time = 5e-3
    data = fake_data(
            # Give the functional form of the fake data.
            # This is an inversion recovery with
            # T₁ of 0.2 
            # T₂* broadening of 50 Hz
            # (echo maximum at echo_time)
            # amplitude of 21
            # resonance offset 100 Hz
            21 * (1 - 2 * s.exp(-vd / 0.2)) * s.exp(+1j * 2 * s.pi * 100 * (t2) - abs(t2) * 50 * s.pi),
            # next we give our dimensions with outer loops first, as they
            # would be acquired on the spectrometer
            # (ordering does matter, because fake_data applies a
            # time-dependent resonance variation -- see fake_data doc.)
            OrderedDict(
                [
                    ("vd", nddata(r_[0:1:40j], "vd")),
                    ("ph1", nddata(r_[0, 2] / 4.0, "ph1")),
                    ("ph2", nddata(r_[0:4] / 4.0, "ph2")),
                    ("t2", nddata(r_[0:0.2:256j] - echo_time, "t2")),
                    ]
                ),
            {"ph1": 0, "ph2": 1},
            )
    data.reorder(["ph1", "ph2", "vd", "t2"])
    DCCT(
        data,
        fl.next("DCCT - time domain"),
        total_spacing=0.15,
        label_spacing_multiplier=55,
        LHS_pad=0.05,
    )
    data.ft("t2")
    DCCT(
        data,
        fl.next("DCCT - frequency domain"),
        total_spacing=0.15,
        label_spacing_multiplier=55,
        LHS_pad=0.05,
    )
    data.ift(["ph1", "ph2"])
    DCCT(
        data,
        fl.next("phase cycling domain"),
        total_spacing=0.15,
        label_spacing_multiplier=55,
        LHS_pad=0.05,
    )
