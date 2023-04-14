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
A, R, nu, t = sp.symbols("A R nu t", real=True)
# create an empty dataset that we will drop the fake data into
thisfit = lmfitdata(nddata(r_[-1:1:1001j], "t"))
thisfit.functional_form = A * sp.exp(-1j*2*pi*nu*t) * sp.exp(-R*sp.pi*abs(t))
only_real = True # if you set this to True, it works
if only_real:
    thisfit.functional_form = sp.re(thisfit.functional_form)
logger.info(strm("Functional Form:", thisfit.functional_form))
# {{{ create the "true" parameters for the fake data by pretending like
#     the true values are a guess, and then outputting the guess data
true_values = {"A": 14.0, "R": 10, "nu": 25}
thisfit.set_guess(true_values)
thisfit.settoguess()
mydata = thisfit.eval()
mydata.add_noise(0.1)
# }}}
# {{{Making guess data
newfit = lmfitdata(mydata)
newfit.functional_form = thisfit.functional_form
newfit.set_guess(
    A=dict(value=13.0, max=20, min=0.0),
    R=dict(value=3, max=1000, min=0),
    nu=dict(value=20),
)
newfit.settoguess()
print("FIRS CALL TO EVAL")
guess = newfit.eval()
# }}}
# {{{ run the fit and generate nddata
# again, now that this is a class, why is this not handled by the fit method?
newfit.fit()
# {{{plot the data with fits and guesses
plt.subplot(211)
plot(mydata, "r", label="data")
print("SECOND CALL TO EVAL")
plot(newfit.eval(), "b", alpha=0.5, label="fit")
plot(guess, "g--", label="guess")
plt.ylabel('real components')
plt.subplot(212)
plot(mydata.imag, "r", label="data")
print("THIRD CALL TO EVAL")
plot(newfit.eval().imag, "b", alpha=0.5, label="fit")
plot(guess.imag, "g--", label="guess")
plt.ylabel('imag components')
# }}}
plt.legend()
plt.show()
