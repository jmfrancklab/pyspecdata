#!/usr/bin/env python3
# vim: set foldmethod=marker :

from pyspecdata import *
from numpy import r_, pi, exp, logspace, sqrt
from sklearn.linear_model import lars_path

init_logging(level="info")

# {{{ load a real cw ESR spectrum

Bname = "$B_0$"

d = find_file("S175R1a.*DHPC.*200304", exp_type="francklab_esr/Sam")
d.chunk_auto("harmonic", "phase")
d = d["harmonic", 0]["phase", 0]

d[Bname] *= 1e-4
d.set_units(Bname, "T")
d.set_ft_initial(Bname, "f").set_ft_prop(Bname, "time_not_aliased")
print(d[Bname])

# A derivative Lorentzian dictionary cannot represent a DC baseline.
# Remove the DC component explicitly before fitting.
d -= d.mean(Bname)

# }}}

# {{{ construct dense Lorentzian-derivative basis using labeled broadcasting

n_center = 800
n_lambda_L = 50

# λ_L is FWHM, not HWHM.
# For Lorentzian absorption:
#     L(B; B_c, λ_L) has HWHM = λ_L/2
#
# Fourier rule:
#     ℱ{L}(u) = exp(−π λ_L |u|) exp(−i 2π u B_c)
#
# Derivative rule:
#     ℱ{∂L/∂B} = i 2π u ℱ{L}
#
# Therefore:
#     A_u(u, B_c, λ_L)
#       = i 2π u exp(−π λ_L |u|) exp(−i 2π u B_c)
#
# This uses pyspecdata broadcasting over:
#     u        conjugate to B₀
#     center   Lorentzian center field B_c
#     λ_L      Lorentzian FWHM

print(d[Bname]);quit()
center = nddata(
    r_[d.getaxis(Bname)[0]:d.getaxis(Bname)[-1]:n_center * 1j],
    "center",
).set_units("center", "T")

lambda_L = nddata(
    logspace(-5, -2.5, n_lambda_L),
    "lambda_L",
).set_units("lambda_L", "T")

u = d.C.ift(Bname).fromaxis(Bname)

A = (
    1j
    * 2
    * pi
    * u
    * exp(-pi * lambda_L * abs(u))
    * exp(-1j * 2 * pi * u * center)
).C.ft(Bname)

# The imaginary component is transform-roundoff; the solver boundary below
# uses the real part.  Keep A as nddata until that boundary.

# Normalize each column:
#     ‖A_i‖₂ = 1
#
# Otherwise the L1 constraint would prefer some linewidths simply because
# the column norm changes with λ_L.
A /= sqrt((A**2).sum(Bname))

# Collapse the physical coefficient grid only after constructing the basis.
# The coefficient vector c still indexes individual (center, λ_L) components.
A.smoosh(["center", "lambda_L"], "basis")

# }}}

# {{{ SVD-compress residual coordinates

# Original constrained problem:
#     minimize ½ ‖A c − y‖₂²
#     subject to c ≥ 0 and 1ᵀc ≤ τ
#
# SVD:
#     A ≈ Uᵣ Σᵣ Vᵣᵀ
#
# Define:
#     Ã = Σᵣ Vᵣᵀ
#     ỹ = Uᵣᵀ y
#
# Compressed problem:
#     minimize ½ ‖Ã c − ỹ‖₂²
#     subject to c ≥ 0
#
# Positive LARS computes the LASSO path:
#     minimize ½ ‖Ã c − ỹ‖₂² + α 1ᵀc
#     with c ≥ 0
#
# Reading off 1ᵀc along that path gives the desired residual-vs-L1 curve.

U, Sigma, Vh = A.C.svd(Bname, "basis")

A_tilde = Sigma * Vh
y_tilde = U @ d

# }}}

# {{{ positive LARS solver boundary

max_active = 400

# Raw arrays appear only at the external solver boundary.
# X has shape:
#     n_samples × n_features = n_SV × n_basis
X = A_tilde.C.reorder(["SV", "basis"]).data.real
y = y_tilde.C.reorder("SV").data.real

alphas, active, coefs = lars_path(
    X,
    y,
    method="lasso",
    positive=True,
    max_iter=max_active,
    return_path=True,
)

# Immediately wrap solver output back into labeled nddata.
coef_path = nddata(coefs, ["basis", "alpha"])
coef_path.setaxis("basis", A_tilde.getaxis("basis"))
coef_path.setaxis("alpha", alphas)

# }}}

# {{{ evaluate path with pyspecdata algebra

fit_path_tilde = A_tilde @ coef_path
resid_path_tilde = fit_path_tilde - y_tilde

residual_norm = sqrt((resid_path_tilde**2).sum("SV"))
l1_norm = coef_path.sum("basis")

# Show the least-regularized point in the path.
show_idx = -1
c_show = coef_path["alpha", show_idx]
fit_show = A @ c_show

# }}}

# {{{ plots

with figlist_var() as fl:
    fl.next("positive LARS path")
    plot(l1_norm, residual_norm, "o-")
    xlabel("positive L1 mass 1ᵀc")
    ylabel("compressed residual norm ‖Ãc − ỹ‖₂")
    title("positive Lorentzian-derivative LASSO path")

    fl.next("fit at end of path")
    plot(d, label="data", alpha=0.7)
    plot(fit_show, label="fit", alpha=0.7)
    plot(d - fit_show, label="residual", alpha=0.7)
    legend()

    fl.show()

# }}}
