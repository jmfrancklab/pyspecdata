"""Using the DCCT function
=======================

Visualize a simulated inversion recovery dataset utilizing the benefits of the
DCCT plotting style.

We can image data in both the phase cycling domain, or the coherence transfer
domain. Artifacts can clearly be discerned from signal in the coherence
transfer domain as well as visualizing the inversion of phase using the domain
colored plotting style. 

Here, kwargs plotted in red (e.g. vert_label_space) illustrate the kwargs are
in display coordinates while kwargs that are in blue (e.g. bbox and LHS_pad)
illustrate the kwargs are in figure coordinates """

from pylab import rcParams
import matplotlib.pyplot as plt
import pyspecdata as psd
from numpy.random import seed
import sympy as s
from collections import OrderedDict

seed(2021)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 2
psd.init_logging(level="debug")

with psd.figlist_var() as fl:
    # provide the symbols that we use for the fake data:
    t2, td, vd, ph1, ph2 = s.symbols("t2 td vd ph1 ph2")
    echo_time = 5e-3
    data = psd.fake_data(
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
                ("vd", psd.nddata(psd.r_[0:1:40j], "vd")),
                ("ph1", psd.nddata(psd.r_[0, 2] / 4.0, "ph1")),
                ("ph2", psd.nddata(psd.r_[0:4] / 4.0, "ph2")),
                ("t2", psd.nddata(psd.r_[0:0.2:256j] - echo_time, "t2")),
            ]
        ),
        {"ph1": 0, "ph2": 1},
    )
    # reorder into a format more suitable for plotting
    data.reorder(["ph1", "ph2", "vd", "t2"])
    # fake_data gives us data already in the coherence domain, so:
    fig = fl.next("Data")  # Make figure object to place the DCCT
    psd.DCCT(data, fig, plot_title="")
    # {{{ add lines indicating kwargs
    # {{{ bbox kwargs
    plt.plot(
        [0, 0.05],
        [0.1, 0.1],
        "b",
        marker="|",
        linewidth=1,
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.text(
        0.008,
        0.12,
        "bbox[0]",
        color="b",
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.plot(
        [0.16, 0.16],
        [0.0, 0.09],
        "b",
        marker="_",
        linewidth=1,
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.text(
        0.17,
        0.03,
        "bbox[1]",
        color="b",
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.plot(
        [0.05, 0.97],
        [0.97, 0.97],
        "b",
        marker="|",
        linewidth=1,
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.text(
        0.45,
        0.98,
        "bbox[2]",
        color="b",
        clip_on=False,
        transform=fig.transFigure,
    )
    # }}}
    # {{{ vert_label_space
    plt.plot(
        [0.105, 0.125],
        [0.5, 0.5],
        "r",
        marker="|",
        linewidth=1,
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.plot(
        [0.07, 0.09],
        [0.5, 0.5],
        "r",
        marker="|",
        linewidth=1,
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.text(
        0.06,
        0.52,
        "kwarg(vert_label_space)",
        color="r",
        clip_on=False,
        transform=fig.transFigure,
    )
    # }}}
    # {{{ gap
    plt.plot(
        [0.5, 0.5],
        [0.5, 0.55],
        "b",
        marker="_",
        linewidth=1,
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.text(
        0.51,
        0.52,
        "kwarg(2*gap)",
        color="b",
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.plot(
        [0.5, 0.5],
        [0.63, 0.655],
        "b",
        marker="_",
        linewidth=1,
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.text(
        0.51,
        0.635,
        "kwarg(gap)",
        color="b",
        clip_on=False,
        transform=fig.transFigure,
    )
    # }}}
