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

for name in T1_names:
    acq = datasets[name].get_prop("acq_params")
    rep_us = acq["FIR_rep_us"] if "FIR_rep_us" in acq else acq["FIR_rep"]
    datasets[name].set_prop("W_k", rep_us * 1e-6 + acq["acq_time_ms"] * 1e-3)
```
<!-- above, you should be setting the power, as well-->

The recovery datasets already have a `tau` axis; `datasets["M_ε"]` already
has a `p` axis.
The `calibration` dictionary supplies fixed globalpars.
<!-- why is this called "calibration"??? that's a bad name-->
Its keys are `C_SL`, `s_max`, and `omega_ratio`.
Include `T_1w0` and `dT_1w` when fixed; omit them when they are fit parameters.
Set `omega_ratio` to `abs(omega_e / omega_H)`.
Each recovery dataset's `p` property supplies a variable,
while its `W_k`
property fixes a dependent localpar for that experiment.
<!-- WK is not a dependent localpar.  It's independent!-->
Keep both as top-level properties set with `.set_prop()`.

Supply all datasets and expressions at construction.
Each dictionary key is the string representation of the symbol on the left
of an equality, and its value is the expression on the right.
Use ordinary SymPy expressions; never construct SymPy equality objects.
With symbols named as shown, the physical model becomes

```python
functional_forms = {
    "M_τ": M_inf_k * (1 - (2 - exp(-W_k * R_1)) * exp(-tau * R_1)),
    "M_ε": M_p0 * (1 - A * s_normed / R_1),
    "A": C_SL * k_sigma * omega_ratio * s_max,
    "R_1": 1 / T_1w + C_SL / k_rho_inv,
    "T_1w": T_1w0 + p * dT_1w,
    "k_rho_inv": k_rho_inv_a + k_rho_inv_b * p + k_rho_inv_c * p**2,
    "s_normed": p / (p + p_half),
}
functional_forms.update(
    {key: sympify(value) for key, value in calibration.items()}
)
fit = lmfitdata(datasets, functional_forms=functional_forms)
```

The fit shares `k_rho_inv_a`, `k_rho_inv_b`, `k_sigma`, `p_half`,
and `M_p0` across the container, with a separate `M_inf_k` for each recovery
dataset at power $p_k$.
Calibration values are constant expressions in `functional_forms`.
<!-- why are you making up vocab still!!!  "calibration" is such a  HORRIBLE word to use here!-->
Thus `C_SL` remains a globalpar when fixed or substituted out.
Likewise, `W_k` remains a dependent localpar when fixed from metadata.
Fix `k_rho_inv_c=0` or fit it according to the chosen inverse-relaxivity model.
Optimize the shared relaxation parameters jointly from the recovery and
magnetization data.

# Vocabulary, classification, and validation

Use the following vocabulary consistently in the plan, implementation,
and `lmfitdata` class docstring so we can discuss the architecture precisely:

- **variables**: coordinates identified from dataset axis coordinate
  labels or from matching unindexed top-level property keys,
  <!-- this is unclear -- what does "unindexed" mean here?-->
  excluding fixed globalpars.
  Axis variables determine the model dimensions.
  Coordinate properties set with `.set_prop()` span multiple datasets,
  such as the powers linking the recovery experiments.
- **funcvals** (function values): symbols whose expressions calculate
  the values used in residuals, such as `M_ε` and `M_τ`.
  They are neither variables nor parameters.
- **parameters**: quantities that may be
  optimized, fixed, or computed;
  the overall goal of optimization is to
  adjust the parameters to choose the best
  physical model possible with the functional
  for that's given.
- **globalpars** (global parameters):
  parameters shared across the
  fitting container.
- **localpars** (local parameters): quantities local to one dataset,
  either optimized independently or computed by substitution.
- **dependent localpars**: quantities computed by substituting supplied
  expressions or experimental metadata,
  such as `R_1(p_k)` at each dataset's
  power `p_k`.
  They are dependent in the sense
  that they depend on other variables and/or
  parameters.
  Do not introduce numbered `R_1_k` symbols.
  <!-- what does this last sentence mean? why?-->
- **independent localpars**: parameters that vary with each dataset.
  As such, among the set of variables and
  parameters, they are independent.
  In order for a given model to pass
  validation,
  independent localpars must be named in a
  particular way:
  specifically, their name must end with an
  index letter, such as `_i`, `_j`, or `_k`.
  For example, `M_inf_k` is the value of the
  parameter we call $M_\infty$.
  There is never a sympy symbol with the name
  `M_inf_k`.
  Rather,
  for each of datasets indexed by $k$ with
  power $p_k$,
  we create a separate variable:
  for example `M_inf_0` at the first power, `M_inf_1` at the second, *etc.*

Classification uses the disjoint classes **variables**, **funcvals**,
**globalpars**, **dependent localpars**, and **independent localpars**.
`parameters` and `localpars` are umbrella terms for the corresponding classes.
Whether a parameter is fixed is separate from its classification.
The classification procedure and validation requirements are enumerated in
stop point 1.

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

For each dataset, evaluate its funcval's substituted expression with axis and
property variables using the existing compact multi-argument lambda
machinery, particularly `fitfunc_multiarg_v2`.
Reshape axis arrays to broadcast correctly for multiple dimensions,
including nonsquare arrays.
Apply the dataset's data and residual transforms according to the existing
convention, subtract the transformed data, and ravel the residual.
View complex residuals as real using the current convention.
Concatenate all datasets' residuals into one 1D real vector.

# Implementation stop points

## Stop point 1: container, docstring, classification, and validation

Implement the container, naming checks, deep copies, delegated chaining,
and decorator-preserving transform storage and access.
Add the vocabulary and architecture TODO notes to the class docstring.
Implement classification using sets of symbols and the vocabulary above:

1. Normalize the single-expression API to `{"z_data": expression}`.
   Let `D` be the set of symbols named by
   dictionary keys, `R` the union
   of the expressions' `free_symbols`, and `S = D | R`.
   Reject distinct symbols with the same string representation.
2. Classify `F = D - R` as funcvals: their expressions are not substituted
   into other supplied expressions.
   Build the dependency graph for `D & R` and reject circular dependencies.
3. Collect axis symbols and unindexed coordinate-property symbols into the
   variables set `V`, excluding symbols already defined in `D`.
   Consider only property names in `S`; ignore unrelated metadata.
   A definition cannot also supply an axis symbol.
   For indexed symbols supplied by properties, such as `W_k`,
   record their dataset bindings as dependent localpars.
   A symbol supplied both by a definition and a property is ambiguous.
4. Among `R - D - V`, classify unbound symbols with dataset-index suffixes
   as independent localpars; classify the remaining unbound symbols as
   globalpars.
   A parameter remains global even if only one funcval uses it.
   If an indexed symbol is bound for some applicable datasets,
   require bindings for all of them rather than optimizing missing values.
5. Follow the dependency graph to classify each symbol in `D & R`.
   An indexed definition, or one depending transitively on variables or
   localpars, is a dependent localpar.
   Other definitions are globalpars, including constants such as `C_SL` and
   expressions of globalpars such as `A`.
   Preserve these classifications when substituting or fixing values.
   Only unbound globalpars and independent localpars enter the optimizer,
   unless explicitly fixed.
6. For each dataset, choose exactly one funcval whose variables and dependent
   localpar bindings are available and whose expression depends on all dataset
   axes.
   Inspect dependencies before substitution, then substitute and rename
   independent localpars for each dataset.
7. Validate that every symbol in `S` has exactly one of the five classes.
   Their union must equal `S`, and pairwise intersections must be empty.
   Every funcval must be used and every required binding must exist.
   Each dataset must have exactly one funcval.
   Failure to classify any symbol unambiguously fails model validation.
   Report the symbols, competing or missing classes, and affected datasets;
   do not fit an invalid model.

Add classification and validation unit tests at this stop point.
For the example, expect funcvals `{M_τ, M_ε}`, variables `{tau, p}`,
independent localpars `{M_inf_k}`, and dependent localpars
`{W_k, R_1, T_1w, k_rho_inv, s_normed}`.
The remaining coefficient and calibration symbols, including `A` and
`C_SL`, are globalpars.
Test fixed values, globalpars used by only one funcval, alternative index
letters, missing bindings, conflicting classifications, ambiguous funcvals,
and circular dependencies.
Defer multi-dataset fitting tests to stop point 3.
Run existing relevant tests and repository lint checks, then stop for review.

## Stop point 2: single-dataset multidimensional fitting

Implement and test multidimensional evaluation and axis broadcasting.
Pack residuals into a real vector using the existing evaluator.
Preserve single-dataset behavior, including complex residuals.
Run relevant tests and lint checks, then stop for review.

## Stop point 3: multi-dataset fitting

Implement and test recursive substitution, classification and
validation, independent localpar renaming, dataset-specific decorators,
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
