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
from lmfitdata import lmfitdata
init_logging(level="debug")
np.random.seed(15816)
# {{{ a lot of what's below depends on knowing what the shape and dimension labels of my data are, so define that here
x_vals = linspace(0, 250, 1500)
empty_data = nddata(x_vals, "x").copy(data=False)
# }}}
thisfit = lmfitdata(empty_data)
# {{{making sympy expression
A, shift, period, decay, x = sp.symbols("A shift period decay x")
thisfit.functional_form = (A * sp.sin(shift + x / period) * sp.exp(-((x * decay) ** 2)))
thisfit.set_guess(
        A = dict(value=13.0, max=20, min=0.0),
        shift = dict(value=0.0, max=pi / 2.0, min=-pi / 2.0),
        period = dict(value=2, max=10),
        decay = dict(value=0.02, max=0.10, min=0.00),
)
thisfit.settoguess()
# }}}
# {{{ nddata to generate the fake data
# {{{ create the "true" parameters for the fake data
true_values = {"A": 14.0, "period": 5.4321, "shift": 0.12345, "decay": 0.01000}
p_true = Parameters()
for k, v in true_values.items():
    p_true.add(k, value=v)
logger.info(strm("p_true is:", p_true))
# }}}
mydata = empty_data.copy(data=False)
mydata.data = thisfit.residual(p_true)
mydata.add_noise(2.8)
# }}}
# {{{ nddata of the guess
guess = empty_data.copy(data=False)
guess.data = thisfit.residual(thisfit.pars)
# }}}
# {{{ run the fit and generate nddata
out = minimize(
    thisfit.residual, thisfit.pars, kws={"data": mydata.data}
)
fit = empty_data.copy(data=False)
fit.data = thisfit.residual(out.params)

# }}}
# {{{ report the fit and generate the plot
report_fit(out, show_correl=True, modelpars=p_true)
plot(mydata, "ro", label="data")
plot(fit, "b", alpha=0.5, label="fit")
plot(guess, "g--", label="guess")
plot(newguess, 'k',label="please work")
# }}}
plt.legend()
plt.show()
