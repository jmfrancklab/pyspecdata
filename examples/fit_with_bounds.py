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
from lmfit import Parameters, fit_report
from lmfit.printfuncs import report_fit
import numpy as np
from pyspecdata import *

init_logging(level="debug")
np.random.seed(15816)
fl = figlist_var()
use_pinvr = False
use_jacobian = True
A, shift, freq, decay, x = sp.symbols("A shift freq decay x", real=True)
# create an empty dataset that we will drop the fake data into
thisfit = lmfitdata(nddata(r_[0:250:1500j], "x"))
thisfit.functional_form = A * sp.sin(shift + x * freq) * sp.exp(-((x * decay) ** 2))
logger.info(strm("Functional Form:", thisfit.functional_form))
# {{{ create the "true" parameters for the fake data by pretending like
#     the true values are a guess, and then outputting the guess data
true_values = {"A": 14.0, "freq": 1/5.4321, "shift": 0.12345, "decay": 0.01000}
thisfit.set_guess(true_values)
thisfit.settoguess()
mydata = thisfit.eval()
mydata.add_noise(2.8)
# }}}
# {{{Making guess data
newfit = lmfitdata(mydata)
newfit.functional_form = thisfit.functional_form
logger.info(strm(newfit.fit_axis,newfit.getaxis(newfit.fit_axis)))
newfit.set_guess(
    A=dict(value=1.0, max=20, min=0.0),
    shift=dict(value=0.0, max=pi / 2.0, min=-pi / 2.0),
    freq=dict(value=0.5, min=0, max=10),
    decay=dict(value=0.02, max=0.10, min=0.00),
)
newfit.settoguess()
logger.info(str(newfit.getaxis(newfit.fit_axis)))
guess = newfit.eval(100)
if use_pinvr:
    for j in range(10):
        newfit.pinvr_step()
    newfit.settoguess()
    newguess = newfit.eval(100)
# }}}
# {{{ run the fit and generate nddata
# again, now that this is a class, why is this not handled by the fit method?
newfit.fit(use_jacobian=use_jacobian) # True is the default -- set False to test following
logger.info(strm("number of function evaluations:",newfit.fit_output.nfev,
                 "(currently gives 78 when use_jacobian is False and 18 when True)"))
# {{{plot the data with fits and guesses
plot(mydata, "ro", label="data")
plot(newfit.eval(100), "b", alpha=0.5, lw=3, label="fit")
plot(guess, "--", alpha=0.5, label="guess")
if use_pinvr:
    plot(newguess, "--", alpha=0.5, lw=2, label="new guess")
plt.gca().text(0.5, 0.75, '$'+sp.latex(newfit.functional_form)+'$',
               transform=plt.gca().transAxes)
# }}}
plt.legend()
plt.figure()
jac = newfit.jacobian(newfit.fit_parameters).view(newfit.data.dtype)
thisline = newfit.eval()
thisline /= thisline.data.max()
plot(thisline, "k", alpha=0.5, label="fit")
for j,thisparam in enumerate(newfit.fit_parameters.keys()):
    thisline = jac[j,:]
    thisline /= thisline.max()
    plot(newfit.getaxis(newfit.dimlabels[0]),thisline,label=f'{thisparam} derivative',
         alpha=0.5)
plt.legend()
logger.info(strm("output:",newfit.output()))
plt.show()
