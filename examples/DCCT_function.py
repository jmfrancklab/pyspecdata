"""Using the DCCT function
=======================

Visualize a simulated inversion recovery dataset
utilizing the benefits of the DCCT plotting style.

We can image data in both the phase cycling domain,
as well as the coherence transfer domain. Artifacts can
clearly be discerned from signal in the coherence transfer
domain as well as visualizing the inversion of phase 
using the domain colored plotting style.
"""

from pylab import *
from pyspecdata import *
import pyspecdata as psd
from pyspecdata.DCCT_function import DCCT
from numpy.random import normal, seed
import sympy as s
from collections import OrderedDict
import sys

seed(2021)
print(sys.path)
print(psd.__file__)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 2
init_logging(level="debug")

with figlist_var() as fl:
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
        21
        * (1 - 2 * s.exp(-vd / 0.2))
        * s.exp(+1j * 2 * s.pi * 100 * (t2) - abs(t2) * 50 * s.pi),
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
    # reorder into a format more suitable for plotting
    data.reorder(["ph1", "ph2", "vd", "t2"])
    # fake_data gives us data already in the coherence domain, so:
    data.ift(["ph1", "ph2"])
    LHS_pad_1 = 0.1  # left hand side padding for first DCCT
    LHS_pad_2 = 0.55  # left hand side padding for second DCCT
    R_to_L = 0.42  # width from left of decorations to right of
    #               plots should be the same for both DCCT
    bbox_bottom = 0.1  # bottom padding should be equal for both
    #                   DCCT
    # keyword arguments to use for the first DCCT
    dcct_kwargs_1 = dict(
        DCCT_bbox=[LHS_pad_1, bbox_bottom, R_to_L], total_spacing=0.2, top_pad=0.1
    )
    fig = fl.next("Data")  # Make figure object to place the DCCT
    data.ft(["ph1", "ph2"])
    DCCT(data, fig, plot_title="Time Domain", **dcct_kwargs_1)
    data.ft("t2")
    # keyword arguments to use for second DCCT
    dcct_kwargs_2 = dict(
        DCCT_bbox=[LHS_pad_2, bbox_bottom, R_to_L], total_spacing=0.2, top_pad=0.1
    )
    DCCT(data, fig, plot_title="Frequency Domain", **dcct_kwargs_2)
    # {{{ add lines indicating kwargs
    #plt.plot(
    #    [0.06, 0.06],
    #    [0.0, 0.09],
    #    "b",
    #    marker="_",
    #    linewidth=1,
    #    clip_on=False,
    #    transform=fig.transFigure,
    #)
    #plt.text(
    #    0.01,
    #    0.03,
    #    "bbox[1]",
    #    color="b",
    #    clip_on=False,
    #    transform=fig.transFigure,
    #)
    #plt.plot(
    #    [0, 0.1],
    #    [0.9, 0.9],
    #    "b",
    #    marker="|",
    #    linewidth=1,
    #    clip_on=False,
    #    transform=fig.transFigure,
    #)
    #plt.text(
    #    0.025,
    #    0.92,
    #    "bbox[0]",
    #    color="b",
    #    clip_on=False,
    #    transform=fig.transFigure,
    #)
    #plt.plot(
    #    [0.1, 0.52],
    #    [0.95, 0.95],
    #    "b",
    #    marker="|",
    #    linewidth=1,
    #    clip_on=False,
    #    transform=fig.transFigure,
    #)
    #plt.text(
    #    0.35,
    #    0.97,
    #    "bbox[2]",
    #    color="b",
    #    clip_on=False,
    #    transform=fig.transFigure,
    #)
    #plt.text(
    #    0.35,
    #    0.97,
    #    "bbox[2]",
    #    color="b",
    #    clip_on=False,
    #    transform=fig.transFigure,
    #)
    #plt.plot(
    #    [0.12, 0.15],
    #    [0.47, 0.47],
    #    "r",
    #    marker="|",
    #    linewidth=1,
    #    clip_on=False,
    #    transform=fig.transFigure,
    #)
    #plt.text(
    #    0.06,
    #    0.49,
    #    "kwarg(vert_label_space)",
    #    color="r",
    #    clip_on=False,
    #    transform=fig.transFigure,
    #)
    ## }}}
