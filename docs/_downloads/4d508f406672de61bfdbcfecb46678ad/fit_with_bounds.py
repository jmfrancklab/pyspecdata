"""
Fit Using Bounds
================

Adapt one of the examples from lmfit to use sympy and pyspecdata.
Eventually, we want to use parts of this in the base library, so we don't have
to retype them every time.

"""

import matplotlib.pyplot as plt
from numpy import exp, linspace, pi, random, sign, sin
import sympy as sp
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
import numpy as np
from pyspecdata import *

init_logging(level="debug")
np.random.seed(15816)
fl = figlist_var()
A, shift, period, decay, x = sp.symbols("A shift period decay x", real=True)
# create an empty dataset that we will drop the fake data into
thisfit = lmfitdata(nddata(r_[0:250:1500j], "x"))
thisfit.functional_form = (
    A * sp.sin(shift + x / period) * sp.exp(-((x * decay) ** 2))
)
logger.info(strm("Functional Form:", thisfit.functional_form))
# {{{ create the "true" parameters for the fake data by pretending like
#     the true values are a guess, and then outputting the guess data
true_values = {"A": 14.0, "period": 5.4321, "shift": 0.12345, "decay": 0.01000}
thisfit.set_guess(true_values)
thisfit.settoguess()
mydata = thisfit.eval()
mydata.add_noise(2.8)
# }}}
# {{{Making guess data
newfit = lmfitdata(mydata)
newfit.functional_form = thisfit.functional_form
newfit.set_guess(
    A=dict(value=13.0, max=20, min=0.0),
    shift=dict(value=0.0, max=pi / 2.0, min=-pi / 2.0),
    period=dict(value=2, max=10),
    decay=dict(value=0.02, max=0.10, min=0.00),
)
newfit.settoguess()
guess = newfit.eval(100)
# }}}
# {{{ run the fit and generate nddata
# again, now that this is a class, why is this not handled by the fit method?
newfit.fit()
# {{{plot the data with fits and guesses
plot(mydata, "ro", label="data")
plot(newfit.eval(100), "b", alpha=0.5, label="fit")
plot(guess, "g--", label="guess")
# }}}
plt.legend()
plt.show()
