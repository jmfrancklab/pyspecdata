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
def list_symbs(f):
    # {{{ this is just to show all the parameters
    list_symbs = []
    for j,k in f.output().items():
        s_repr = sp.latex(sp.Symbol(j))
        list_symbs.append(f'${s_repr} = {k:0.5g}$')
    list_symbs = '\n'.join(list_symbs)
    # }}}
    return list_symbs
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
    # {{{ use fitdata
    fitdata_instance = fitdata(fake_data)
    fitdata_instance.functional_form = functional_form
    logger.info(strm("Functional Form", fitdata_instance.functional_form))
    logger.info(strm("Functional Form", fitdata_instance.functional_form))
    fitdata_instance.set_guess({M0: -500, Mi: 500, R1: 2})
    fitdata_instance.settoguess()
    fl.next("fit with guess")
    fl.plot(fake_data, "o", alpha=0.5, label="fake data")
    fl.plot(fitdata_instance.eval(100), alpha=0.5, label="fitdata guess")
    fitdata_instance.fit()
    print('-'*5,"Results for fitdata:",'-'*5)
    print("output:", fitdata_instance.output())
    print("latex:", fitdata_instance.latex())
    T1 = 1.0 / fitdata_instance.output("R_1")
    # }}}
    # {{{ lmfitdata method
    lmfitdata_instance = lmfitdata(fake_data)
    lmfitdata_instance.functional_form = functional_form
    lmfitdata_instance.set_guess(
        M_0=dict(value=-500, max=0, min=-501),
        M_inf=dict(value=500, max=501, min=0),
        R_1=dict(value=5, max=6, min=1),
    )
    lmfitdata_instance.settoguess()
    fl.plot(lmfitdata_instance.eval(100), alpha=0.5, label="lmfitdata guess")
    lmfitdata_instance.fit()
    print('-'*5,"Results for lmfitdata:",'-'*5)
    print("output:", lmfitdata_instance.output())
    print("latex:", lmfitdata_instance.latex())
    T1 = 1.0 / lmfitdata_instance.output("R_1")
    # }}}
    fitdata_line = fl.plot(fitdata_instance.eval(100), alpha=0.5, label="fit data fit")
    lmfitdata_line = fl.plot(
        lmfitdata_instance.eval(100), ":", alpha=0.5, linewidth=3, label="lmfitdata fit"
    )
    # {{{ just put the text
    ax = gca()
    # {{{ lmfitdata
    text(
        0.6,
        0.5,
        "lmfitdata RESULT: %s" % lmfitdata_instance.latex(),
        ha="center",
        va="center",
        color=lmfitdata_line[0].get_color(),
        transform=ax.transAxes,
    )
    text(0.6,0.5,(3*'\n')+list_symbs(lmfitdata_instance),
            ha='center',va='top',
            size=10,
            color=fitdata_line[0].get_color(),
            transform = ax.transAxes)
    # }}}
    # {{{ fitdata
    text(
        0.6,
        0.25,
        "fitdata RESULT: %s" % fitdata_instance.latex(),
        ha="center",
        va="center",
        color=fitdata_line[0].get_color(),
        transform=ax.transAxes,
    )
    text(0.6,0.25,(3*'\n')+list_symbs(fitdata_instance),
            ha='center',va='top',
            size=10,
            color=fitdata_line[0].get_color(),
            transform = ax.transAxes)
    # }}}
    # }}}
