"""Verify that lmfitdata.jacobian() agrees with numerical finite differences.

The helper assert_symbolic_matches_numerical_jacobian is ported from
~/git_repos/FLinst/tests/test_plot_current_rounding_helpers.py, which has the
most complete version of this test infrastructure.

For smooth expressions (exponential, trig) all Jacobian columns are checked
with assert_allclose.  For Heaviside we check the dA column analytically and
verify the dt0 spike structure separately because the discrete DiracDelta
approximation is not in the linear finite-difference regime for sub-grid
steps.
"""
import copy

import numpy as np
import sympy as sp

import pyspecdata as psd


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _build_fit(axis, dim_name, expression, guesses):
    fit = psd.lmfitdata(
        psd.nddata(np.zeros_like(axis), [dim_name]).setaxis(dim_name, axis)
    )
    fit.functional_form = expression
    fit.set_guess(**guesses)
    fit.set_to_guess()
    return fit


def assert_symbolic_matches_numerical_jacobian(
    fit_object, params, *, rtol, atol, label
):
    """Compare symbolic Jacobian against adaptive finite-difference Jacobian.

    Ported from
    ~/git_repos/FLinst/tests/test_plot_current_rounding_helpers.py.
    Uses an adaptive step-size search to stay in the linear regime.
    """

    def adaptive_directional_derivative(
        baseline_vector, baseline_scale, base_value, shifted_vector
    ):
        step = max(1e-8, 1e-4 * max(abs(base_value), 1.0))
        for _ in range(20):
            deltas = []
            slopes = []
            for multiple in (1.0, 2.0, 3.0):
                delta = (
                    shifted_vector(base_value + multiple * step)
                    - baseline_vector
                )
                deltas.append(delta)
                slopes.append(delta / (multiple * step))
            if (
                max(np.linalg.norm(d) for d in deltas)
                < 1e3 * np.finfo(float).eps * baseline_scale
            ):
                step *= 8.0
                continue
            slope_mean = np.mean(slopes, axis=0)
            mismatch = (
                max(np.linalg.norm(s - slope_mean) for s in slopes)
                / max(np.linalg.norm(slope_mean), 1e-12)
            )
            if mismatch < 1e-2:
                return slope_mean
            step *= 0.25
        raise AssertionError("could not find a linear finite-difference step")

    sigma = fit_object.get_error()
    baseline = fit_object.residual(params, sigma)
    baseline_scale = max(np.linalg.norm(baseline), 1.0)
    symbolic_jacobian = fit_object.jacobian(params)
    numerical_columns = []
    for name in fit_object.parameter_names:
        base_value = params[name].value

        def shifted_vector(shifted_value, *, name=name):
            shifted = copy.deepcopy(params)
            shifted[name].value = shifted_value
            shifted.update_constraints()
            return fit_object.residual(shifted, sigma)

        try:
            numerical_columns.append(
                adaptive_directional_derivative(
                    baseline, baseline_scale, base_value, shifted_vector
                )
            )
        except AssertionError as exc:
            raise AssertionError(
                f"{label}: could not find a linear finite-difference step "
                f"for {name}"
            ) from exc

    np.testing.assert_allclose(
        symbolic_jacobian,
        np.array(numerical_columns),
        rtol=rtol,
        atol=atol,
        err_msg=f"Jacobian mismatch at {label}",
    )


# ---------------------------------------------------------------------------
# exponential: A * exp(-R * t)
# ---------------------------------------------------------------------------

def test_jacobian_exponential():
    t, A, R = sp.symbols("t A R", real=True)
    expr = A * sp.exp(-R * t)
    t_arr = np.linspace(0.1, 3.0, 80)
    fit = _build_fit(
        t_arr, "t", expr, dict(A=dict(value=2.0), R=dict(value=1.5))
    )
    assert_symbolic_matches_numerical_jacobian(
        fit,
        copy.deepcopy(fit.guess_parameters),
        rtol=3e-4,
        atol=3e-4,
        label="exponential A*exp(-R*t)",
    )


# ---------------------------------------------------------------------------
# Heaviside: A * Heaviside(t - t0)
#
# dA column  — Heaviside is piecewise-constant in A so the adaptive helper
#              finds a clean linear regime.
# dt0 column — uses finite_difference_heaviside_derivative (a discrete
#              DiracDelta approximation); sub-grid shifts produce no residual
#              change at all, so the adaptive helper cannot find a linear
#              regime.  We instead check structural properties: exactly one
#              nonzero entry, at the correct grid location, with the correct
#              sign.
# ---------------------------------------------------------------------------

def test_jacobian_heaviside_dA():
    """d/dA of A*Heaviside(t) matches the numerical derivative."""
    t, A = sp.symbols("t A", real=True)
    expr = A * sp.Heaviside(t)
    t_arr = np.linspace(-1.0, 1.0, 101)
    fit = _build_fit(t_arr, "t", expr, dict(A=dict(value=1.0)))
    assert_symbolic_matches_numerical_jacobian(
        fit,
        copy.deepcopy(fit.guess_parameters),
        rtol=3e-4,
        atol=3e-4,
        label="Heaviside dA",
    )


def test_jacobian_heaviside_dt0_spike_location():
    """d/dt0 of A*Heaviside(t-t0): symbolic spike must be at the grid point
    nearest t0 and must be negative (increasing t0 lowers the step)."""
    t, A, t0 = sp.symbols("t A t0", real=True)
    expr = A * sp.Heaviside(t - t0)
    t_arr = np.linspace(-1.0, 1.0, 101)
    t0_val = 0.0
    A_val = 1.0
    fit = _build_fit(
        t_arr,
        "t",
        expr,
        dict(A=dict(value=A_val), t0=dict(value=t0_val)),
    )
    fit.residual(fit.guess_parameters)  # populate nan_mask
    jac = fit.jacobian(fit.guess_parameters)

    param_names = list(fit.guess_parameters.valuesdict().keys())
    i_t0 = param_names.index("t0")
    col = jac[i_t0]

    nonzero_idx = np.flatnonzero(col)
    assert len(nonzero_idx) == 1, (
        f"expected one nonzero entry in dt0 column, got {len(nonzero_idx)}"
    )
    spike_t = t_arr[nonzero_idx[0]]
    nearest_t = t_arr[np.argmin(np.abs(t_arr - t0_val))]
    assert spike_t == nearest_t, (
        f"spike at t={spike_t}, expected nearest grid point to "
        f"t0={t0_val} which is t={nearest_t}"
    )
    assert col[nonzero_idx[0]] < 0, "dt0 spike should be negative for A > 0"
