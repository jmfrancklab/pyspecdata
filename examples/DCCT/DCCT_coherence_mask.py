"""
Showing a Coherence Mask on a DCCT Plot
=======================================

When we process phase cycled data,
we frequently want to keep only some of the coherence pathways
(typically the pathway of the desired signal, and sometimes
a pathway that contains an artifact that we want to keep track of),
and to mask out all the others.

Here, we show a simulated inversion recovery dataset,
and use `coherence_unmask_fn` to cover the coherence pathways that are
masked out with a hatch.
The function that we pass is the same type of function that's used by
the correlation alignment in pyspecProcScripts:
it's handed an nddata that has only the coherence transfer dimensions and
is filled with zeros,
and it sets the pathways that we do *not* want to mask to 1.
"""

from pylab import rcParams
import pyspecdata as psd
from numpy.random import seed
import sympy as s
from collections import OrderedDict

seed(2021)
rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 1
signal_pathway = {"ph1": 0, "ph2": 1}


def coherence_unmask_fn(coh_array):
    """Unmask the pathway of the signal, as well as the pathway where
    neither pulse changes the coherence order (which contains, *e.g.*,
    receiver offsets)."""
    for thispathway in [
        coh_array.get_prop("coherence_pathway"),
        {"ph1": 0, "ph2": 0},
    ]:
        coh_array["ph1", thispathway["ph1"]]["ph2", thispathway["ph2"]] = 1
    return coh_array


with psd.figlist_var() as fl:
    # provide the symbols that we use for the fake data:
    t2, td, vd, ph1, ph2 = s.symbols("t2 td vd ph1 ph2")
    data = psd.fake_data(
        # inversion recovery with
        # T₁ of 0.2
        # T₂* broadening of 50 Hz
        # (echo maximum 5 ms after the start of acquisition)
        # amplitude of 21
        # resonance offset 100 Hz
        21
        * (1 - 2 * s.exp(-vd / 0.2))
        * s.exp(+1j * 2 * s.pi * 100 * (t2) - abs(t2) * 50 * s.pi),
        OrderedDict(
            [
                ("vd", psd.nddata(psd.r_[0:1:40j], "vd")),
                ("ph1", psd.nddata(psd.r_[0, 2] / 4.0, "ph1")),
                ("ph2", psd.nddata(psd.r_[0:4] / 4.0, "ph2")),
                ("t2", psd.nddata(psd.r_[0:0.2:256j] - 5e-3, "t2")),
            ]
        ),
        signal_pathway,
    )
    data.set_prop("coherence_pathway", signal_pathway)
    # reorder into a format more suitable for plotting
    data.reorder(["ph1", "ph2", "vd", "t2"])
    data.ft("t2")
    psd.DCCT(
        data,
        fl.next("Coherence mask"),
        coherence_unmask_fn=coherence_unmask_fn,
    )
