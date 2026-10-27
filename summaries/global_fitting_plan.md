---
title: "Global fitting architecture plan"
subtitle: "lmfitdata as a container for one or more nddata objects"
date: "October 26, 2026"
lang: en-US
documentclass: article
fontsize: 11pt
geometry: margin=1in
colorlinks: true
toc: true
toc-depth: 2
numbersections: true
header-includes:
  - \DeclareUnicodeCharacter{03B5}{\ensuremath{\varepsilon}}
  - \DeclareUnicodeCharacter{03C4}{\ensuremath{\tau}}
---

# Objective

Refactor `lmfitdata` into a container for one or more `nddata` objects,
rather than simply inheriting a single `nddata` object.
Fit across dimensions and datasets by evaluating symbolic expressions and
concatenating their residuals into one real vector.
Preserve existing single-dataset behavior and reuse the current symbolic
evaluation machinery.
The interfaces below are the implementation plan.

# Example of a target physical model

The global fit combines magnetization versus power, $M_\epsilon(p)$,
with FIR recovery curves measured at powers $p_j$.
Use $T_1$ for the relaxation time; $R_1=1/T_1$ is its inverse.
The shared relaxation model is

$$
\begin{aligned}
T_{1,w}(p) &= T_{1,w}(0)+p\,\Delta T_{1,w}, \\
k_{\rho,\mathrm{inv}}(p)
&= \frac{C_{SL}}{T_1(p)^{-1}-T_{1,w}(p)^{-1}}
= a+bp+cp^2, \\
T_1(p) &= \left[
\frac{1}{T_{1,w}(p)}+\frac{C_{SL}}{k_{\rho,\mathrm{inv}}(p)}
\right]^{-1}.
\end{aligned}
$$

Here $C_{SL}$ is the spin-label concentration.
Fix $c=0$ for the linear inverse-relaxivity model in
[fit_ODNP_data.py](../../proc_scripts/examples/broken/fit_ODNP_data.py);
allow $c$ to vary for the quadratic extension.
The water calibration supplies $T_{1,w}(0)$ and $\Delta T_{1,w}$.

$$
\begin{aligned}
s(p) &= \frac{p}{p+p_{1/2}}, \\
M_\epsilon(p) &= M_\epsilon(0)\left[
1-C_{SL}k_\sigma s(p)
\left|\frac{\omega_e}{\omega_H}\right|T_1(p)
\right].
\end{aligned}
$$

This is the source's magnetization expression
$M_0-M_0 A\,s(p)/R_1(p)$, with
$A=C_{SL}k_\sigma|\omega_e/\omega_H|$.
Use the concentration-normalized $k_\sigma$ convention and a positive,
dimensionless frequency ratio.
The saturation law has unit limiting saturation; an unknown limiting
saturation is absorbed into the fitted effective $k_\sigma$.

The FIR model from
[proc_FIR.py](../../proc_scripts/examples/proc_FIR.py) includes the finite
repetition delay:

$$
M_{\tau,j}(\tau)=M_{\tau,j}(\infty)
\left[1-\left(2-e^{-W_j/T_1(p_j)}\right)e^{-\tau/T_1(p_j)}\right].
$$

$$
W_j=10^{-6}\,\mathrm{FIR\_rep}_j
+10^{-3}\,\mathrm{acq\_time\_ms}_j.
$$

Use seconds for $W_j$, $\tau$, and $T_1$.
The acquisition metadata stores `FIR_rep` in microseconds and
`acq_time_ms` in milliseconds.
The endpoint values $M_\epsilon(0)$ and $M_{\tau,j}(\infty)$
are scalar fitting parameters, represented below by `M_p0` and `M_inf_i`.
They also equal their expressions evaluated at $p=0$ and
$\tau\to\infty$, respectively.
The FIR initial value is determined by the model:

$$
M_{\tau,j}(0)=M_{\tau,j}(\infty)
\left[e^{-W_j/T_1(p_j)}-1\right].
$$

# Datasets and functional forms

Accept a dictionary of datasets or a dictionary-like HDF5 group directly.
For example, select already loaded experimental datasets:

```python
datasets = {
    "T1_0": loaded["T1_0"],
    "T1_1": loaded["T1_1"],
    "M_ε": loaded["M_ε"],
}
T1_names = ("T1_0", "T1_1")

for name, dataset in datasets.items():
    dataset.name(name)
    for key, value in calibration.items():
        dataset.set_prop(key, value)

for name in T1_names:
    datasets[name].set_prop("p", powers[name])
    acq = datasets[name].get_prop("acq_params")
    datasets[name].set_prop(
        "W", acq["FIR_rep"] * 1e-6 + acq["acq_time_ms"] * 1e-3
    )
```

The recovery datasets already have a `tau` axis; `datasets["M_ε"]` already
has a `p` axis.
The `calibration` dictionary supplies fixed top-level properties.
Its keys are `C_SL`, `T_1w0`, `dT_1w`, and `omega_ratio`.
Set `omega_ratio` to `abs(omega_e / omega_H)`.
Likewise, each recovery dataset's power and repetition delay are top-level
properties, not nested under a `variables` property.

Supply all datasets and expressions at construction.
Each dictionary key is the string representation of the symbol on the left
of an equality, and its value is the expression on the right.
Use ordinary SymPy expressions; never construct SymPy equality objects.
With symbols named as shown, the physical model becomes

```python
functional_forms = {
    "M_τ": M_inf_i * (1 - (2 - exp(-W / T_1)) * exp(-tau / T_1)),
    "M_ε": M_p0 * (1 - C_SL * k_sigma * s * omega_ratio * T_1),
    "T_1": 1 / (1 / T_1w + C_SL / k_rho_inv),
    "T_1w": T_1w0 + p * dT_1w,
    "k_rho_inv": a + b * p + c * p**2,
    "s": p / (p + p_half),
}
fit = lmfitdata(datasets, functional_forms=functional_forms)
```

The fit shares `a`, `b`, `k_sigma`, `p_half`, and `M_p0` across the
container, with a separate `M_inf_i` for each recovery dataset.
Fix `c=0` or fit it according to the chosen inverse-relaxivity model.
Optimize the shared relaxation parameters jointly from the recovery and
magnetization data.

# Vocabulary and automatic identification of variables

Define the following vocabulary in the `lmfitdata` class docstring:

- **variables**: fixed symbols identified from dataset axis labels or
  matching top-level property keys.
  Axis variables determine the model dimensions; scalar properties fix
  values for an individual dataset.
- **parameters**: symbols altered by the optimizer.
- **globalpars** (global parameters): optimized parameters shared across
  the fitting container.
- **localpars** (local parameters): quantities local to one dataset,
  either optimized independently or computed by substitution.
- **dependent localpars**: quantities computed by substituting the supplied
  expressions, such as `T_1` evaluated at each dataset's power.
  Do not introduce numbered `T_1_i` symbols.
- **independent localpars**: per-dataset optimized parameters with trailing
  `_i`, such as `M_inf_i`.
  Rename these for each dataset, giving `M_inf_0` and `M_inf_1`.
  Raise an error if an independent localpar lacks `_i`.

The sorting routine must:

1. Normalize the existing single-expression API to
   `functional_forms = {"z_data": expression}`.
2. Substitute dependent expressions recursively.
   Use dictionary keys to identify symbols and reject circular dependencies.
   Expressions referenced by other expressions supply substitutions.
   The remaining expressions describe the measured data.
3. Identify variables from each dataset's axes and top-level properties.
   Use a measured-data expression only when its variables are available
   and it depends on all of that dataset's axes.
   Thus `M_τ` describes the recovery datasets with a `tau` axis and fixed
   `p` and `W`, while `M_ε` describes the dataset with a `p` axis.
   Require exactly one expression for each dataset; report missing or
   ambiguous choices.
4. Rename independent localpars for each dataset and classify all symbols
   using the vocabulary above.

For multi-dataset construction, print each dataset's variables and its
expression after substitution and renaming.
Also print the globalpars and both kinds of localpars.

# Container behavior and transforms

Preserve all delegated `nddata` attributes, methods, and properties for a
single dataset, including chaining such as `fit.ft("t").ift("t")`.
For multiple datasets, select a dataset explicitly:

```python
fit["M_ε"].getaxis("p")
fit["T1_0"].getaxis("tau")
```

Mapping keys must agree with the contained datasets' `nddata.name()` values.
For sequence input, derive keys from those names and require every dataset
to be named.
Both `fit.copy()` and `fit.C` deep-copy the container, including datasets,
symbolic state, parameters, and transforms.
Single-dataset `eval()` returns an `nddata`; multi-dataset `eval()` returns
a dictionary with the input dataset keys.

Keep the existing single-dataset transform decorators:

```python
@fit.define_residual_transform
def residual_domain_transform(d):
    return d.ft("t")
```

For multiple datasets, expose decorators through dataset selection:

```python
@fit["M_ε"].define_residual_transform
def magnetization_residual_transform(d):
    return d

@fit["T1_0"].define_data_transform
def recovery_data_transform(d):
    return d
```

Store transforms in container-level dictionaries keyed by dataset name.
Reuse or retool `lmfitdata` behavior to support dataset selection and
these decorators; do not introduce a separate `FitMember` class.

# Residual construction

For each dataset, evaluate its substituted expression with axis and
property variables using the existing compact multi-argument lambda
machinery, particularly `fitfunc_multiarg_v2`.
Reshape axis arrays to broadcast correctly for multiple dimensions,
including nonsquare arrays.
Apply the dataset's data and residual transforms according to the existing
convention, subtract the transformed data, and ravel the residual.
View complex residuals as real using the current convention.
Concatenate all datasets' residuals into one 1D real vector.

# Implementation stop points

## Stop point 1: container, docstring, and sorting routine

Implement the container, naming checks, deep copies, delegated chaining,
and decorator-preserving transform storage and access.
Normalize the functional-form dictionary and implement the sorting routine
and docstring vocabulary above; include architecture TODO notes.
This stage must support multiple datasets in principle, including sorting,
while deferring multi-dataset tests to stop point 3.
Run existing relevant tests and repository lint checks, then stop for review.

## Stop point 2: single-dataset multidimensional fitting

Implement and test multidimensional evaluation and axis broadcasting.
Pack residuals into a real vector using the existing evaluator.
Preserve single-dataset behavior, including complex residuals.
Run relevant tests and lint checks, then stop for review.

## Stop point 3: multi-dataset fitting

Implement and test recursive substitution, automatic identification of
variables, independent localpar renaming, dataset-specific decorators,
residual concatenation, dictionary `eval()` output, and the construction
summary showing each dataset's substituted expression.
Cover mapping input, top-level property variables, and errors for ambiguous
expressions or missing names.
Run relevant tests and lint checks, then stop for review.

## Stop point 4: synthetic ODNP test and plotted example

Fit noisy synthetic FIR recovery and magnetization-versus-power datasets
jointly to recover `k_sigma` with the physical model above.
Use known calibration values, positive relaxation times and inverse
relaxivity over the fitted power range, and the stated saturation convention.
Add a matching example plotting the data and fitted curves.
Run relevant tests and lint checks.

# Building the PDF

From the repository root:

```bash
source /home/jmfranck/base/bin/activate
export PATH="/usr/local/texlive/2024/bin/x86_64-linux:$PATH"
pandoc summaries/global_fitting_plan.md \
  --from=markdown --standalone --pdf-engine=pdflatex \
  --output=summaries/global_fitting_plan.pdf
```

Omit the `export` if `pdflatex` is already on `PATH`.
