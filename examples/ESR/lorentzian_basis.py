#!/usr/bin/env python3
"""
fit ESR to Lorentzian basis
===========================

λ_L is FWHM, not HWHM.
For Lorentzian absorption:
    L(B; B_c, λ_L) has HWHM = λ_L/2

Fourier rule:
    ℱ{L}(u) = exp(−π λ_L |u|) exp(−i 2π u B_c)

Derivative rule:
    ℱ{∂L/∂B} = i 2π u ℱ{L}

Therefore:
    A_u(u, B_c, λ_L)
      = i 2π u exp(−π λ_L |u|) exp(−i 2π u B_c)

This uses pyspecdata broadcasting over:
    u        conjugate to B₀
    center   Lorentzian center field B_c
    λ_L      Lorentzian FWHM

The SVD compression works like this:
====================================

Original constrained problem:
    minimize ½ ‖A c − y‖₂²
    subject to c ≥ 0 and 1ᵀc ≤ τ

SVD:
    A ≈ Uᵣ Σᵣ Vᵣᵀ

Define:
    Ã = Σᵣ Vᵣᵀ
    ỹ = Uᵣᵀ y

Compressed problem:
    minimize ½ ‖Ã c − ỹ‖₂²
    subject to c ≥ 0

Positive LARS computes the LASSO path:
    minimize ½ ‖Ã c − ỹ‖₂² + α 1ᵀc
    with c ≥ 0

Reading off 1ᵀc along that path gives the desired residual-vs-L1 curve.

"""
# vim: set foldmethod=marker :

from pyspecdata import *
from numpy import r_, pi, exp, logspace, sqrt, log10
from matplotlib.pyplot import title, xlabel, ylabel, legend
from sklearn.linear_model import lars_path
import os

# {{{ changeable parameters
Bname = "$B_0$"
esr_file = "15N_S175R1a_pR_DHPC_today_200304.DSC"
# Use a tiny basis first so that we can inspect the functions and debug fast.
preview_n_center = 7
preview_n_lambda_L = 4
# The dense basis can be made larger again once we know everything is correct.
fit_n_center = 80
fit_n_lambda_L = 8
# }}}


def build_lorentzian_basis(
    d,
    Bname,
    n_center,
    n_lambda_L,
    center_limits=None,
):
    x = d.getaxis(Bname)
    # go for 5x the pixel size (decayed to 0 at end), b/c otherwise, we get weird discretization issues
    lambda_L_limits = ((x[1]-x[0])*5, (x[-1]-x[0])/2)
    if center_limits is None:
        center_limits = (x[0], x[-1])
    center = nddata(
        r_[center_limits[0] : center_limits[1] : n_center * 1j],
        "center",
    ).set_units("center", "T")
    lambda_L = nddata(
        logspace(log10(lambda_L_limits[0]), log10(lambda_L_limits[1]), n_lambda_L),
        "lambda_L",
    ).set_units("lambda_L", "T")
    u = d.C.ift(Bname).fromaxis(Bname)
    A = (
        -1j
        * 2
        * pi
        * u
        * exp(-pi * lambda_L * abs(u))
        * exp(1j * 2 * pi * u * center)
    )
    A[Bname,0] *= 0.5 # u=0 Heaviside issue
    A.ft(Bname)
    A.setaxis(Bname, x)
    A.set_units(Bname, d.get_units(Bname))
    # Normalize each column:
    #     ‖A_i‖₂ = 1
    #
    # Otherwise the L1 constraint would prefer some linewidths simply because
    # the column norm changes with λ_L.
    A /= sqrt((abs(A) ** 2).sum(Bname))
    return A


init_logging(level="info")

# {{{ load a real cw ESR spectrum

d = find_file(
    esr_file, exp_type="francklab_esr/Sam"
)
d.chunk_auto("harmonic", "phase")
d = d["harmonic", 0]["phase", 0]

d[Bname] *= 1e-4
d.set_units(Bname, "T")
d.set_ft_initial(Bname, "f").set_ft_prop(Bname, "time_not_aliased")

# A derivative Lorentzian dictionary cannot represent a DC baseline.
# Remove the DC component explicitly before fitting.
d -= d.C.mean(Bname)

# }}}

# {{{ plots

with figlist_var() as fl:
    # {{{ construct dense Lorentzian-derivative basis using labeled broadcasting
    preview_A = build_lorentzian_basis(
        d, Bname, preview_n_center, preview_n_lambda_L
    )
    print("preview basis shape", preview_A.shape)
    A = build_lorentzian_basis(d, Bname, fit_n_center, fit_n_lambda_L)

    # The imaginary component is transform-roundoff; the solver boundary below
    # uses the real part.  Keep A as nddata until that boundary.

    # Collapse the physical coefficient grid only after constructing the basis.
    # The coefficient vector c still indexes individual (center, λ_L) components.
    A.smoosh(["center", "lambda_L"], "basis")
    # }}}

    fl.next("reduced basis preview")
    for center_j in range(preview_n_center):
        for lambda_j in range(preview_n_lambda_L):
            plot(
                preview_A["center", center_j]["lambda_L", lambda_j],
                alpha=0.35,
                human_units=False,
            )
    title("reduced Lorentzian-derivative basis")

    # {{{ SVD-compress residual coordinates
    U, Sigma, Vh = A.C.svd(Bname, "basis")
    U.set_units(Bname, d.get_units(Bname))
    A_tilde = Sigma * Vh
    y_tilde = U.C.reorder(["SV", Bname]).along(Bname) @ d
    print(
        "during compression, d was reduced from",
        d.shape,
        "to",
        y_tilde.shape,
    )
    # }}}

    # {{{ positive LARS solver boundary
    print("beginning LARS path")
    alphas, active, coefs = lars_path(
        A_tilde.C.reorder(["SV", "basis"]).data.real,
        y_tilde.C.reorder("SV").data.real,
        method="lasso",
        positive=True,
        max_iter=400,
        return_path=True,
    )
    print("done with LARS path")
    # Immediately wrap solver output back into labeled nddata.
    coef_path = nddata(coefs, ["basis", "alpha"])
    coef_path.setaxis("basis", A_tilde.getaxis("basis"))
    coef_path.setaxis("alpha", alphas)
    # }}}

    # {{{ evaluate path with pyspecdata algebra
    # Show the least-regularized point in the path.
    fit_show = A.C.along("basis") @ coef_path.C["alpha", -1]
    # }}}

    fl.next("positive LARS path")
    plot(
        coef_path.C.sum("basis").data,
        sqrt(
            (
                abs((A_tilde.C.along("basis") @ coef_path) - y_tilde) ** 2
            ).sum("SV")
        ).data,
        "o-",
    )
    xlabel("positive L1 mass 1ᵀc")
    ylabel("compressed residual norm ‖Ãc − ỹ‖₂")
    title("positive Lorentzian-derivative LASSO path")

    fl.next("fit at end of path")
    plot(d, label="data", alpha=0.7)
    plot(fit_show, label="fit", alpha=0.7)
    plot(d - fit_show, label="residual", alpha=0.7)
    legend()
# }}}
