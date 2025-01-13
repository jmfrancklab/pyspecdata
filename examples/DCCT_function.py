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
# {{{ kwargs for DCCT plot
bbox = [0.05,0.1,0.85,0.75]
horiz_label_spacer = 50
gap = 0.1
# }}}

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
    fig = fl.next("Data")  # Make figure object to place the DCCT
    ax_list, allow_for_labels, total_scale_transform, ax0_origin = psd.DCCT(
            data, fig, 
            horiz_label_spacer = horiz_label_spacer, 
            gap = gap, 
            bbox=bbox, 
            plot_title="")
    # {{{ add lines indicating kwargs
    # {{{ bbox kwargs
    plt.plot(
        [0, bbox[0]],
        [0.1, 0.1],
        "b",
        marker="|",
        linewidth=1,
        clip_on=False,
        transform=fig.transFigure
    )
    plt.text(
        bbox[0]/7,
        0.12,
        "bbox[0]",
        color="b",
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.plot(
        [0.16, 0.16],
        [0.0, bbox[1]],
        "b",
        marker="_",
        linewidth=1,
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.text(
        0.17,
        bbox[1]/3,
        "bbox[1]",
        color="b",
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.plot(
        [bbox[0], bbox[2]+bbox[0]],
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
    plt.plot(
        [0.93, 0.93],
        [bbox[1], bbox[1]+bbox[3]+gap],
        "b",
        marker="_",
        linewidth=1,
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.text(
        0.95,
        0.5,
        "bbox[3]",
        color="b",
        clip_on=False,
        transform=fig.transFigure,
    )
    # }}}
    # {{{ horiz_label_space
    plt.plot(
        [-horiz_label_spacer, -2*horiz_label_spacer],
        [0.5, 0.5],
        "r",
        marker="|",
        linewidth=1,
        clip_on=False,
        transform=ax0_origin,
    )
    plt.plot(
        [0.0, -horiz_label_spacer],
        [0.55, 0.55],
        "r",
        marker="|",
        linewidth=1,
        clip_on=False,
        transform=ax0_origin,
    )
    plt.text(
        -3*horiz_label_spacer,
        0.52,
        "kwarg(horiz_label_space)",
        color="r",
        clip_on=False,
        transform=ax0_origin,
    )
    # }}}
    # {{{ gap
    ax4_x,ax4_y = (ax_list[4].transAxes+fig.transFigure.inverted()).transform(psd.r_[0.5,1])
    ax3_x,ax3_y = (ax_list[3].transAxes+fig.transFigure.inverted()).transform(psd.r_[0.5,1])
    plt.plot(
        [ax3_x, ax3_x],
        [ax3_y, ax3_y + gap/2],
        "b",
        marker="_",
        linewidth=1,
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.text(
        ax3_x+0.01,
        ax3_y+0.01,
        r"kwarg(gap) / $\text{nPh}_{\text{outer}}$",
        color="b",
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.plot(
        [ax4_x, ax4_x],
        [ax4_y, ax4_y + gap/4],
        "b",
        marker="_",
        linewidth=1,
        clip_on=False,
        transform=fig.transFigure,
    )
    plt.text(
        ax4_x+0.01,
        ax4_y + 0.007,
        r"kwarg(gap) / $\text{nPh}_{\text{inner}}$",
        color="b",
        clip_on=False,
        transform=fig.transFigure,
    )
    # }}}
    # }}}
