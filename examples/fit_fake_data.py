"""Fitting Data (Nonlinear + Symbolic)
======================================

This example creates fake data with noise
then fits the exponential with the fitdata
function."""
from pyspecdata import *
import sympy as sp

# {{{ this is the contents of pylab.py -- works
# need to go through and figure out which lines
# are actually needed and which are not
# -- I have already stripped out some
from lmfit import Parameters, minimize
from matplotlib.pyplot import figure, subplot, show, xlim, ylim, plot, gca
from numpy import *  # I think it wasn't importing from numpy b/c it seems we're inside sphinx

# }}}
fl = figlist_var()
# {{{creating a fake data recovery curve
tau = nddata(r_[0:2:100j], "tau")
fake_data = 102 * (1 - 2 * exp(-tau * 6.0))
fake_data.add_noise(5.0)
# }}}
# {{{ define the expression of the functional form once, and then use it
#     for both types of classes
M0, Mi, R1, vd = sp.symbols("M_0 M_inf R_1 tau", real=True)
functional_form = Mi + (M0 - Mi) * sp.exp(-vd * R1)
# }}}
with figlist_var() as fl:
    # {{{ fitting data
    f = fitdata(fake_data)
    f.functional_form = functional_form
    logger.info(strm("Functional Form", f.functional_form))
    logger.info(strm("Functional Form", f.functional_form))
    f.set_guess({M0: -500, Mi: 500, R1: 2})
    f.settoguess()
    fl.next("fit with guess")
    fl.plot(fake_data, "o", alpha=0.5, label="fake data")
    fl.plot(f.eval(100), alpha=0.5, label="fitdata guess")
    f.fit()
    print("output:", f.output())
    print("latex:", f.latex())
    # }}}
    T1 = 1.0 / f.output("R_1")
    # }}}
    # }}}
    # {{{lmfitdata method
    newfit = lmfitdata(fake_data)
    newfit.functional_form = functional_form
    newfit.set_guess(
        M_0=dict(value=-500, max=0, min=-501),
        M_inf=dict(value=500, max=501, min=0),
        R_1=dict(value=5, max=6, min=1),
    )
    newfit.settoguess()
    fl.plot(newfit.eval(100), alpha=0.5, label="lmfitdata guess")
    newfit.fit()
    # }}}
    thisline = fl.plot(f.eval(100), alpha=0.5, label="fit data fit")
    thatline = fl.plot(
        newfit.eval(100), ":", alpha=0.5, linewidth=3, label="lmfitdata fit"
    )
    # {{{ just put the text
    ax = gca()
    text(
        0.6,
        0.5,
        "LMFIT RESULT: %s" % newfit.latex(),
        ha="center",
        va="center",
        color=thatline[0].get_color(),
        transform=ax.transAxes,
    )
    text(
        0.6,
        0.25,
        "FITDATA RESULT:%s" % f.latex(),
        ha="center",
        va="center",
        color=thisline[0].get_color(),
        transform=ax.transAxes,
    )
    # }}}
