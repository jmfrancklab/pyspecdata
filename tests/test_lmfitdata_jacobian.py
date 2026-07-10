"""Verify lmfitdata Jacobians against numerical finite differences."""

import copy
import time

import numpy as np
import pytest
import sympy as sp
from scipy.special import jv

import pyspecdata as psd
from pyspecdata import Q_


# {{{ helpers
def build_fit(axis, dim_name, expression, guesses):
    fit = psd.lmfitdata(
        psd.nddata(np.zeros_like(axis), [dim_name]).setaxis(dim_name, axis)
    )
    fit.functional_form = expression
    fit.set_guess(guesses)
    fit.set_to_guess()
    return fit


def assert_symbolic_matches_numerical_jacobian(
    fit_object, params, *, rtol, atol, label
):
    """Compare symbolic Jacobian against adaptive finite differences."""

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
                max(np.linalg.norm(delta) for delta in deltas)
                < 1e3 * np.finfo(float).eps * baseline_scale
            ):
                step *= 8.0
                continue
            slope_mean = np.mean(slopes, axis=0)
            mismatch = max(
                np.linalg.norm(slope - slope_mean) for slope in slopes
            ) / max(np.linalg.norm(slope_mean), 1e-12)
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
        if not params[name].vary:
            continue
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


# }}}


# {{{ tests
def test_functional_form_is_write_once():
    x, A = sp.symbols("x A", real=True)
    fit = build_fit(
        np.linspace(0, 1, 8),
        "x",
        A * x,
        {"A": {"value": 1.0}},
    )
    with pytest.raises(ValueError, match="fresh lmfitdata object"):
        fit.functional_form = 2 * A * x


def test_jacobian_exponential():
    t, A, R = sp.symbols("t A R", real=True)
    fit = build_fit(
        np.linspace(0.1, 3.0, 80),
        "t",
        A * sp.exp(-R * t),
        {"A": {"value": 2.0}, "R": {"value": 1.5}},
    )
    assert_symbolic_matches_numerical_jacobian(
        fit,
        copy.deepcopy(fit.guess_parameters),
        rtol=3e-4,
        atol=3e-4,
        label="exponential A*exp(-R*t)",
    )


def test_jacobian_returns_only_varying_parameter_rows():
    x, A, B = sp.symbols("x A B", real=True)
    fit = build_fit(
        np.linspace(-1.0, 1.0, 41),
        "x",
        A * x + B,
        {"A": {"value": 2.0}, "B": {"value": -0.3}},
    )
    fit.guess_parameters["B"].vary = False
    assert_symbolic_matches_numerical_jacobian(
        fit,
        copy.deepcopy(fit.guess_parameters),
        rtol=3e-4,
        atol=3e-4,
        label="linear expression with fixed offset",
    )


def test_emma_voigt_transform_jacobian_matches_numerical():
    B_axis = np.linspace(-8, 8, 2048)
    demo = psd.nddata(np.zeros_like(B_axis), "B").labels("B", B_axis)
    mod_amp = 0.0
    B = sp.symbols("B", real=True)
    A, lambda_L, Bcenter, sigma = sp.symbols(
        "A lambda_L Bcenter sigma", real=True
    )
    voigt_line = (
        A
        * (-1j * 2 * np.pi * B)
        * sp.exp(1j * 2 * np.pi * Bcenter * B)
        * sp.exp(
            -lambda_L * sp.pi * abs(B) - sp.pi**2 * abs(B) ** 2 * sigma**2
        )
    )
    A_symbols = sp.symbols("A0:1", real=True)
    A_disp_symbols = sp.symbols("A_disp0:1", real=True)
    Bcenter_symbols = sp.symbols("Bcenter0:1", real=True)
    FWHM_symbols = sp.symbols("FWHM0:1", real=True)
    L_vs_G_frac_symbols = sp.symbols("L_vs_G_frac0:1", real=True)
    voigt_fwhm_coeff = 0.5346
    voigt_fwhm_coeff_symb = sp.Float(voigt_fwhm_coeff)
    voigt_fwhm_remainder_symb = (sp.Integer(1) - voigt_fwhm_coeff_symb) ** 2

    def build_model(amplitudes):
        thefunction = 0
        for j, amplitude in enumerate(amplitudes):
            lorentzian_FWHM = FWHM_symbols[j] * L_vs_G_frac_symbols[j]
            gaussian_FWHM = sp.sqrt(
                (FWHM_symbols[j] - voigt_fwhm_coeff_symb * lorentzian_FWHM)
                ** 2
                - voigt_fwhm_remainder_symb * lorentzian_FWHM**2
            )
            thefunction += voigt_line.subs(
                {A: amplitude, Bcenter: Bcenter_symbols[j]}
            ).subs(
                {
                    lambda_L: lorentzian_FWHM,
                    sigma: gaussian_FWHM / (2 * sp.sqrt(sp.log(2))),
                }
            )
        return thefunction

    demo.ift("B", shift=True)
    demo = psd.lmfitdata(demo)
    demo.functional_form = build_model(A_symbols) - sp.I * sp.sign(
        B
    ) * build_model(A_disp_symbols)

    @demo.define_data_transform
    def my_data_transform(d_local):
        d_local["B":0] *= 0.5
        d_local.ft("B")
        return d_local.real

    @demo.define_residual_transform
    def my_residual_transform(d_local):
        h_m = Q_(mod_amp, "G").to("G").magnitude
        d_local *= d_local.fromaxis(
            "B",
            lambda B_axis: jv(0, abs(h_m * B_axis * np.pi))
            + jv(2, abs(h_m * B_axis * np.pi)),
        )
        d_local["B":0] *= 0.5
        d_local.ft("B")
        return d_local.real

    demo.set_guess(
        {
            "A0": {"value": 1.2},
            "A_disp0": {"value": -0.35},
            "Bcenter0": {"value": 0.12},
            "FWHM0": {"value": 1.1},
            "L_vs_G_frac0": {"value": 0.45, "min": 0, "max": 1},
        }
    )
    demo.set_to_guess()
    assert set(demo.parameter_names) == {
        "A0",
        "A_disp0",
        "Bcenter0",
        "FWHM0",
        "L_vs_G_frac0",
    }
    assert_symbolic_matches_numerical_jacobian(
        demo,
        copy.deepcopy(demo.guess_parameters),
        rtol=2e-3,
        atol=2e-3,
        label="Emma one-component transformed Voigt",
    )


def test_transformed_voigt_fit_is_faster_with_jacobian():
    rng = np.random.default_rng(20260710)
    B_axis = np.linspace(-8, 8, 2048)
    B = sp.symbols("B", real=True)
    A, lambda_L, Bcenter, sigma = sp.symbols(
        "A lambda_L Bcenter sigma", real=True
    )
    voigt_line = (
        A
        * (-1j * 2 * np.pi * B)
        * sp.exp(1j * 2 * np.pi * Bcenter * B)
        * sp.exp(
            -lambda_L * sp.pi * abs(B) - sp.pi**2 * abs(B) ** 2 * sigma**2
        )
    )
    A_symbols = sp.symbols("A0:3", real=True)
    Bcenter_symbols = sp.symbols("Bcenter0:3", real=True)
    FWHM_symbols = sp.symbols("FWHM0:3", real=True)
    L_vs_G_frac_symbols = sp.symbols("L_vs_G_frac0:3", real=True)
    voigt_fwhm_coeff_symb = sp.Float(0.5346)
    voigt_fwhm_remainder_symb = (sp.Integer(1) - voigt_fwhm_coeff_symb) ** 2
    thefunction = 0
    for j, amplitude in enumerate(A_symbols):
        lorentzian_FWHM = FWHM_symbols[j] * L_vs_G_frac_symbols[j]
        gaussian_FWHM = sp.sqrt(
            (FWHM_symbols[j] - voigt_fwhm_coeff_symb * lorentzian_FWHM) ** 2
            - voigt_fwhm_remainder_symb * lorentzian_FWHM**2
        )
        thefunction += voigt_line.subs(
            {A: amplitude, Bcenter: Bcenter_symbols[j]}
        ).subs(
            {
                lambda_L: lorentzian_FWHM,
                sigma: gaussian_FWHM / (2 * sp.sqrt(sp.log(2))),
            }
        )
    source = psd.nddata(np.zeros_like(B_axis), "B").labels("B", B_axis)
    source.ift("B", shift=True)
    template = psd.lmfitdata(source.C)
    template.functional_form = thefunction
    true_values = {
        "A0": 1.0,
        "A1": 0.8,
        "A2": 1.15,
        "Bcenter0": -2.4,
        "Bcenter1": 0.05,
        "Bcenter2": 2.35,
        "FWHM0": 0.9,
        "FWHM1": 1.05,
        "FWHM2": 0.95,
        "L_vs_G_frac0": 0.35,
        "L_vs_G_frac1": 0.55,
        "L_vs_G_frac2": 0.45,
    }
    raw_clean = template.fitfunc_multiarg_v2(
        *(template.getaxis(k) for k in template.variable_names),
        **true_values,
    )
    noise_level = 1e-3 * np.max(np.abs(raw_clean))
    source.data = raw_clean + noise_level * (
        rng.normal(size=raw_clean.shape)
        + 1j * rng.normal(size=raw_clean.shape)
    )
    results = {}
    for use_jacobian in (False, True):
        fit = psd.lmfitdata(source.C)
        fit.functional_form = thefunction

        @fit.define_data_transform
        def my_data_transform(d_local):
            d_local["B":0] *= 0.5
            d_local.ft("B")
            return d_local.real

        @fit.define_residual_transform
        def my_residual_transform(d_local):
            h_m = Q_(0.3, "G").to("G").magnitude
            d_local *= d_local.fromaxis(
                "B",
                lambda B_axis: jv(0, abs(h_m * B_axis * np.pi))
                + jv(2, abs(h_m * B_axis * np.pi)),
            )
            d_local["B":0] *= 0.5
            d_local.ft("B")
            return d_local.real

        guesses = {}
        for j in range(3):
            guesses[f"A{j}"] = {
                "value": true_values[f"A{j}"] * (1.0 + 0.08 * (-1) ** j)
            }
            guesses[f"Bcenter{j}"] = {
                "value": true_values[f"Bcenter{j}"] + 0.08 * (-1) ** j,
                "min": true_values[f"Bcenter{j}"] - 0.5,
                "max": true_values[f"Bcenter{j}"] + 0.5,
            }
            guesses[f"FWHM{j}"] = {
                "value": true_values[f"FWHM{j}"] * 1.12,
                "min": 0.3,
                "max": 2.0,
            }
            guesses[f"L_vs_G_frac{j}"] = {
                "value": min(true_values[f"L_vs_G_frac{j}"] + 0.08, 0.95),
                "min": 0,
                "max": 1,
            }
        fit.set_guess(guesses)
        start = time.perf_counter()
        fit.fit(use_jacobian=use_jacobian)
        elapsed = time.perf_counter() - start
        residual_norm = np.linalg.norm(fit.residual(fit.fit_parameters))
        results[use_jacobian] = (elapsed, residual_norm)
    no_jacobian_time, no_jacobian_residual = results[False]
    jacobian_time, jacobian_residual = results[True]
    assert jacobian_residual <= no_jacobian_residual * (1 + 1e-4), (
        f"jacobian residual {jacobian_residual} is worse than "
        f"non-jacobian residual {no_jacobian_residual}"
    )
    assert jacobian_time <= no_jacobian_time / 1.5, (
        f"jacobian fit took {jacobian_time:.6g} s, but needed to be at most "
        f"{no_jacobian_time / 1.5:.6g} s based on non-jacobian time "
        f"{no_jacobian_time:.6g} s"
    )


# }}}
