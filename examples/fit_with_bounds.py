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
# {{{ a lot of what's below depends on knowing what the shape and dimension labels of my data are, so define that here
empty_data = nddata(r_[0:250:1500j], "x")
# }}}
# {{{making sympy expression
A, shift, period, decay, x = sp.symbols("A shift period decay x", real=True)
thisfit = lmfitdata(empty_data)
thisfit.functional_form = A * sp.sin(shift + x / period) * sp.exp(-((x * decay) ** 2))
logger.info(strm("Functional Form:", thisfit.functional_form))
# }}}
# {{{ create the "true" parameters for the fake data
true_values = {"A": 14.0, "period": 5.4321, "shift": 0.12345, "decay": 0.01000}
p_true = Parameters()
for k, v in true_values.items():
    p_true.add(k, value=v)
logger.info(strm("p_true is:", p_true))
# }}}
# {{{ nddata to generate the fake data
mydata = empty_data.copy(data=False)
mydata.data = thisfit.residual(p_true)
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
