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
with FIR curves measured at recovery times $\tau_j$ and powers $p_k$.
(Currently, an FIR at each power is stored in a different dataset).
Use $T_1$ for the relaxation time; $R_1=1/T_1$ is its inverse.
The shared relaxation model is

$$
\begin{aligned}
T_{1,w}(p) &= T_{1,w}(0)+p\,\Delta T_{1,w}, \\
k_{\rho,\mathrm{inv}}(p)
&= \frac{C_{SL}}{R_1(p)-T_{1,w}(p)^{-1}} \\
&= k_{\rho,\mathrm{inv},a}
+ k_{\rho,\mathrm{inv},b}p+k_{\rho,\mathrm{inv},c}p^2, \\
R_1(p) &=
\frac{1}{T_{1,w}(p)}+\frac{C_{SL}}{k_{\rho,\mathrm{inv}}(p)}.
\end{aligned}
$$

Here $C_{SL}$ is the spin-label concentration.
It is a **constant** -- *i.e.* it is fixed for all datapoints of all datasets.
Fix $k_{\rho,\mathrm{inv},c}=0$ for the linear inverse-relaxivity model in
[fit_ODNP_data.py](../../proc_scripts/examples/broken/fit_ODNP_data.py);
allow $k_{\rho,\mathrm{inv},c}$ to vary for the quadratic extension.
$T_{1,w}(0)$ and $\Delta T_{1,w}$ may either be constants, or fit parameters,
depending on implementation.

$$
\begin{aligned}
s(p) &= \frac{s_{max} p}{p+p_{1/2}}, \\
M_\epsilon(p) &= M_\epsilon(0)\left[
1-\frac{C_{SL}k_\sigma s(p)
\left|\omega_e/\omega_H\right|}{R_1(p)}
\right].
\end{aligned}
$$

Because all the terms on the right hand side of
$A=C_{SL}k_\sigma|\omega_e/\omega_H| s_{max}$
except $k_\sigma$ are constants,
we can rewrite this as

$$
M_\epsilon(p_k) = M_\epsilon(0)\left[
1-\frac{A s_{normed}(p_k)}{R_1(p_k)}
\right].
$$

where $s_{normed}(p)=p/(p+p_{1/2})$ does not include $s_{max}$.

The FIR model from
[proc_FIR.py](../../proc_scripts/examples/proc_FIR.py) includes the finite
repetition delay:

$$
M_{\tau}(\tau_j,p_k)=M_{\tau,k}(\infty)
\left[1-\left(2-e^{-W_k R_1(p_k)}\right)e^{-\tau_j R_1(p_k)}\right].
$$

$$
W_k=10^{-6}\,\mathtt{FIR\_rep\_us}_k
+10^{-3}\,\mathtt{acq\_time\_ms}_k.
$$

Use seconds for $W_k$, $\tau$, and $T_1$.
Use inverse seconds for $R_1$.
The acquisition metadata stores `FIR_rep_us` (or legacy `FIR_rep`)
in microseconds and `acq_time_ms` in milliseconds.
The endpoint values $M_\epsilon(0)$ and $M_{\tau,k}(\infty)$
are scalar fitting parameters, represented below by `M_p0` and `M_inf_k`.
They also equal their expressions evaluated at $p=0$ and
$\tau\to\infty$, respectively.
Here $j$ indexes recovery times and $k$ indexes powers; the recovery
asymptote and $W_k$ belong to the experiment at power $p_k$.

Note that we typically write code for the equations above as given above,
to preserve legibility and ease the process of setting up global fitting.
We rely on the `lmfitdata` module to perform the substitution and identify the
types of parameters (more on this below).
After substitution, there are only two equations:
one for $M_{\tau}(\tau_j,p_k)$ that describes the FIR experiments,
and one for $M_\epsilon(p_k)$
that describes progressive saturation.

# Datasets and functional forms

Accept a dictionary of datasets or a dictionary-like HDF5 group directly.
For example, select already loaded experimental datasets:

```python
datasets = {
    "T1_0": loaded["T1_0"],
    "T1_1": loaded["T1_1"],
    "M_ε": loaded["M_ε"],
}
```

Select recovery experiments by name using a regular expression.
Assume the loader has stored each experiment's power in watts as a
top-level `p` property via `.set_prop("p", power_in_watts)`.
This supplies a power coordinate linking the recovery datasets.

```python
import re

T1_names = [name for name in datasets if re.fullmatch(r"(?:T1|FIR)_.*", name)]

for name, dataset in datasets.items():
    dataset.name(name)
    for key, value in calibration.items():
        dataset.set_prop(key, value)

for name in T1_names:
    acq = datasets[name].get_prop("acq_params")
    rep_us = acq["FIR_rep_us"] if "FIR_rep_us" in acq else acq["FIR_rep"]
    datasets[name].set_prop("W", rep_us * 1e-6 + acq["acq_time_ms"] * 1e-3)
```

The recovery datasets already have a `tau` axis; `datasets["M_ε"]` already
has a `p` axis.
The `calibration` dictionary supplies fixed top-level properties.
Its keys are `C_SL`, `s_max`, and `omega_ratio`.
Include `T_1w0` and `dT_1w` when fixed; omit them when they are fit parameters.
Set `omega_ratio` to `abs(omega_e / omega_H)`.
Likewise, each recovery dataset's power and repetition delay are top-level
properties, not nested under a `variables` property.

Supply all datasets and expressions at construction.
Each dictionary key is the string representation of the symbol on the left
of an equality, and its value is the expression on the right.
Use ordinary SymPy expressions; never construct SymPy equality objects.
With symbols named as shown, the physical model becomes

<!-- in the following, W varies with each different power, so it should be called W_k-->
<!-- (in essence, when we set things up, we need to treat it as a dependent localvar.  We can later substitute it out or fix it ot the experimental parameters-->
<!-- similarly, C_SL is a globalpar, but it will be fixed or subbed out-->

```python
functional_forms = {
    "M_τ": M_inf_k * (1 - (2 - exp(-W * R_1)) * exp(-tau * R_1)),
    "M_ε": M_p0 * (1 - A * s_normed / R_1),
    "A": C_SL * k_sigma * omega_ratio * s_max,
    "R_1": 1 / T_1w + C_SL / k_rho_inv,
    "T_1w": T_1w0 + p * dT_1w,
    "k_rho_inv": k_rho_inv_a + k_rho_inv_b * p + k_rho_inv_c * p**2,
    "s_normed": p / (p + p_half),
}
fit = lmfitdata(datasets, functional_forms=functional_forms)
```

The fit shares `k_rho_inv_a`, `k_rho_inv_b`, `k_sigma`, `p_half`,
and `M_p0` across the container, with a separate `M_inf_k` for each recovery
dataset at power $p_k$.
Fix `k_rho_inv_c=0` or fit it according to the chosen inverse-relaxivity model.
Optimize the shared relaxation parameters jointly from the recovery and
magnetization data.

# Vocabulary and automatic identification of variables

Use the following vocabulary consistently in the plan, implementation,
and `lmfitdata` class docstring so we can discuss the architecture precisely:

- **variables**: fixed symbols identified from dataset axis labels or
  matching top-level property keys.
  Axis variables determine the model dimensions.
  Properties set with `.set_prop()` supply axis coordinates spanning
  multiple datasets, such as the powers linking the recovery experiments.
- **parameters**: symbols altered by the optimizer.
- **globalpars** (global parameters): optimized parameters shared across
  the fitting container.
- **localpars** (local parameters): quantities local to one dataset,
  either optimized independently or computed by substitution.
- **dependent localpars**: quantities computed by substituting the supplied
  expressions, such as `R_1` evaluated at each dataset's power.
  Do not introduce numbered `R_1_k` symbols.
- **independent localpars**: per-dataset optimized parameters with a trailing
  index letter, such as `_i`, `_j`, or `_k`.
  Here `M_inf_k` uses the dataset index $k$ from the math.
  The index $j$ spans recovery times within a dataset.
  Rename these for each dataset, giving `M_inf_0` and `M_inf_1`.
  Identify the dataset index letter automatically and raise an error if an
  independent localpar lacks the corresponding suffix.

The sorting routine must:

1. Normalize the existing single-expression API to
   `functional_forms = {"z_data": expression}`.
2. Identify variables from each dataset's axes and top-level properties.
   Among expressions not used as substitutions, choose the expression
   whose variables are available and which depends on all dataset axes.
   Inspect the supplied dependent expressions to identify variables such
   as `p` before substitution.
   Thus `M_τ` describes the recovery datasets with a `tau` axis and fixed
   `p` and `W`, while `M_ε` describes the dataset with a `p` axis.
   Require exactly one expression for each dataset; report missing or
   ambiguous choices.
3. Substitute dependent expressions recursively.
   Use dictionary keys to identify symbols and reject circular dependencies.
   Substitute `R_1`, `T_1w`, `k_rho_inv`, `A`, and `s_normed` into the
   expressions for `M_τ` and `M_ε`.
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
