---
title: "Global fitting architecture plan"
subtitle: "lmfitdata as a container for one or more nddata objects"
date: "Plan revised July 14, 2026; recovered October 7; corrected October 26, 2026"
lang: en-US
documentclass: article
fontsize: 11pt
geometry: margin=1in
colorlinks: true
toc: true
toc-depth: 2
numbersections: true
---

# Conversation source

This document recovers the final plan from the VS Code Codex conversation
**Sketch global fitting architecture**. The conversation began on
July 9, 2026 at 11:47 p.m. America/New_York and continued through July 14.
The final plan was recorded on July 14, 2026 at 1:24 p.m. America/New_York
(17:24:22 UTC).

Session identifier:

```text
019f4a22-fd3c-7e51-81be-1d47db0d8e00
```

The source is a local JSONL file in
`/home/jmfranck/.codex/sessions/2026/07/09/`, named:

```text
rollout-2026-07-09T23-47-14-
019f4a22-fd3c-7e51-81be-1d47db0d8e00.jsonl
```

The two lines above form one filename without a newline.

The final plan is lightly edited for standalone use. It incorporates the
last corrections about decorator transforms and top-level dataset properties.
Earlier requirements about method chaining, supplying datasets and equations
at construction, and fitting across both dimensions and datasets are made
explicit. Two wording slips in the final assistant summary are corrected:
`.C` means a deep copy, and dependent local parameters are computed quantities,
rather than independently optimized parameters. These corrections follow your
instructions in the conversation; the current `nddata.C` implementation also
returns `self.copy()`.

The document was rechecked against all seven proposed-plan revisions and your
successive corrections on October 7, 2026. The governing revision is the final
July 14 plan. Earlier requirements are retained where later corrections did
not supersede them, including implementing the sorter at stop point 1 and
showing each dataset's substituted expression.

The October 26 corrections below use magnetization, the ODNP and FIR models
from `proc_scripts`, dataset mappings, and string-to-expression dictionaries.
They supersede the corresponding recovered examples.

<!-- I don't care about any of this historical crap preceding.  We are looking for a concise accounting of the IDEA as a plan that can be ACTED ON-->

This is a recovered design plan, not a claim that the proposed API has been
implemented. The Python snippets illustrate the proposed interfaces.

# Summary

Refactor `lmfitdata` into a container that can hold one or more `nddata` objects while preserving current single-dataset behavior. The architecture should support fitting across dimensions and/or across datasets by evaluating symbolic functional forms, raveling every modeled/data residual into one 1D real vector, and optimizing global plus local parameters.

Current single-dataset examples must continue to behave the same. Global fitting intentionally spans both multiple dimensions of one dataset and multiple datasets; the container provides the common abstraction. Multi-dataset behavior is added in staged stop points and reviewed before being broadly tested.

Refactor conservatively from the existing architecture, drawing on the older
branch history and comments while reusing the current symbolic evaluation
machinery.

# Public API and compatibility

- With one dataset, existing usage remains valid:

  ```python
  fit.ft(...)
  fit.getaxis("t")
  fit.C
  fit.copy()
  ```

- With multiple datasets, direct delegated `nddata` behavior is ambiguous and must use dictionary-like dataset selection:

  ```python
  fit["M"].ft(...)
  fit["T1_0"].getaxis("tau")
  ```

- Delegation covers all `nddata` attributes, methods, and properties. Mutating
  methods that return the contained dataset must preserve container chaining,
  such as `fit.ft("t").ift("t")`.

- Accept a dictionary of datasets, including a dictionary-like HDF5 group,
  directly. Mapping keys identify datasets and must agree with their
  `nddata.name()` values. For sequence input, keys come from each contained
  `nddata.name()`; every dataset must already have a name.

- `.copy()` and `.C` are container operations, not delegated dataset operations.
  - `fit.copy()` deep-copies the whole fitting container.
  - `fit.C` returns a deep-copied fitting container, preserving symbolic state, parameters, transforms, and contained datasets.

- Single-dataset `eval()` returns an `nddata`, as now.
- Multi-dataset `eval()` returns a dictionary:

  ```python
  {
      "T1_0": fitted_T1_0,
      "T1_1": fitted_T1_1,
      "M": fitted_M,
  }
  ```

# Nomenclature to add to the `lmfitdata` class docstring

The `lmfitdata` class docstring should define the core vocabulary up front:

- `variables`: symbols defining model coördinates/abscissae. They are not altered by optimization. They come from dataset axis labels or top-level dataset properties whose names match symbols in the functional forms. Axis variables determine the dimensions of the model array that is raveled into the residual vector; a top-level scalar property fixes a variable for that dataset.
- `parameters`: symbols altered by the optimizer.
- `globalpars`: short for “global parameters”; optimized parameters shared across the whole fit.
- `localpars`: quantities local to one dataset, either optimized independently or generated by substitution; document that `localpars` means “local parameters.”
- `dependent localpars`: local quantities generated by symbolic substitution. Example: `T_1` is replaced by its power-dependent expression below, without creating `T_1_i` symbols.
- `independent localpars`: per-dataset optimized parameters written with trailing `_i`, e.g. `M_0_i`, `M_inf_i`, then renamed to dataset-specific symbols such as `M_0_0`, `M_0_1`, etc.

If the sorter determines that a symbol is an independent localpar but it lacks trailing `_i`, raise an error.

# Variables and dataset metadata

Variables are discovered from:

- axis labels:

  ```python
  datasets["M"].setaxis("p", p_axis)
  ```
  <!-- previous should be M_ε rather than M (throughout); also the user typically wouldn't setaxis, since that would come from the experimental data, so this is a little confusing.-->
  <!-- also, the actual math showing he the mathematical model, in real math, should come before this, so the reader knows what they are looking at!-->

- top-level properties whose names match symbols:

  ```python
  for name in T1_names:
      datasets[name].set_prop("p", powers[name])
  ```

Do not use nested props like:

```python
set_prop("variables", {"p": p0})
```

For each dataset, the sorter should treat a symbol as a variable if it matches either a dataset axis label or a top-level property key.

# Functional forms

Functional forms are supplied as a dictionary. All datasets and relevant equations must be supplied when creating a multi-dataset container.

Existing single-expression API is preserved by internally normalizing to:

```python
functional_forms = {"z_data": expression}
```

New multi-form API:

<!-- the following was not a good change. You should have rather shown an example of what datasets looked like using a dict construction, and you should have kept something like the functional forms.   EXCEPT that the actual math (as $$ math) should have been given first!!! -->
```python
fit = lmfitdata(datasets, functional_forms=functional_forms)
```

Here `datasets` is a dictionary or dictionary-like HDF5 group of datasets;
`functional_forms` is the dictionary shown in the ODNP example below.
Each functional-form key is the string representation of the symbol on the
left of an equality; its value is the expression on the right. Use ordinary
SymPy expressions throughout; never construct SymPy equality objects.
Keys do not manually assign forms to datasets. Dataset-to-form assignment is
auto-determined by matching variables in the fully substituted expression
against each dataset's axes and top-level props. A data model must depend on
all of that dataset's axes; scalar properties supply additional fixed inputs.
Report ambiguous assignments rather than choosing silently.

<!-- the following is vervbose, and it uses the "dataset-to-form assignment" which is a confusing and opaque term!-->

The sorter should:

- normalize old `.functional_form = expr` / single-expression construction to `{"z_data": expr}`;
- distinguish observable model roots (`M_tau`, `M_epsilon` below) from dependent symbols referenced by their expressions (`T_1`, `k_rho_inv`, `s`);
<!-- wtf is a model root???? again with making up terms!-->
- recursively substitute dependent expressions using the dictionary keys, rejecting cycles;
- assign each dataset exactly one data-model expression by variable matching;
- expand independent localpars ending in `_i`;
- classify symbols as variables, globalpars, dependent localpars, or independent localpars.

For multi-dataset construction only, print a stdout summary showing:

- each dataset and its variables;
- which functional form was assigned to each dataset;
- the post-substitution/post-renaming expression used for each residual;
- all variables, globalpars, dependent localpars, and independent localpars.

# Per-dataset transforms

Transforms must remain decorator-friendly.

Single-dataset behavior stays as now:

```python
@fit.define_residual_transform
def residual_domain_transform(d):
    return d.ft("t")
```

For multi-dataset fits, transforms are specified through dictionary-like dataset access:

```python
@fit["M"].define_residual_transform
def M_residual_transform(d):
    return d

@fit["T1_0"].define_data_transform
def T1_0_data_transform(d):
    return d
```

Implementation should store these in container-level dictionaries keyed by dataset name, while exposing the decorator API through the selected dataset view. This should be done by reusing/retooling `lmfitdata` behavior, not by introducing a separate `FitMember` class.

Residual construction applies each dataset’s own data/residual transform before raveling and concatenating.

# Residual construction

All fits produce one 1D real residual vector.

For each dataset:

- evaluate its assigned symbolic model using the existing compact multi-argument lambda machinery, especially the current `fitfunc_multiarg_v2` pattern;
- use axis variables and top-level prop variables as the fixed variable inputs;
- apply that dataset’s data/residual transforms;
- subtract transformed data;
- ravel;
- view complex residuals as real, following the current convention.

For multi-dataset fits, concatenate the per-dataset residual vectors.

# ODNP target behavior

Use $T_1$ for the longitudinal relaxation time; $R_1(p)=1/T_1(p)$ is
only its inverse, not an additional fitted quantity. The executable model in
`~/git_repos/proc_scripts/examples/broken/fit_ODNP_data.py` (the
`krho_inv`, `R1p_expr`, and `integral_vs_p_fit.functional_form` assignments)
gives the following relations, written here in terms of $T_1$:

$$
T_{1,w}(p)=T_{1,w}(0)+p\,\Delta T_{1,w},
\qquad
k_{\rho,\mathrm{inv}}(p)
=\frac{C_{SL}}{T_1(p)^{-1}-T_{1,w}(p)^{-1}}
=a+bp+cp^2,
$$

<!-- what you're doing here doesn't make sense.
it actually does make sense to define k_{\rho,inv}(p), just to make things look easier,
but then you should use it in the $T_1(p)$ equation.  The same with $T_{1,w}(p)$-->

$$
T_1(p)=\left[
\frac{1}{T_{1,w}(0)+p\,\Delta T_{1,w}}
+\frac{C_{SL}}{a+bp+cp^2}
\right]^{-1}.
$$

<!-- again, the following is just verbose nonsense -- as noted earlier, it also comes out of order-->
Here $C_{SL}$ is the spin-label concentration. The source fits the inverse
relaxivity to a **line**, so $c=0$ reproduces its model; the optional quadratic
coefficient follows the supplied sketch. Thus $a=k_{\rho,\mathrm{inv}}(0)$,
$b$ is its linear power coefficient, and $c$ its quadratic coefficient (with
no implicit factorial). The water calibration supplies $T_{1,w}(0)$ and
$\Delta T_{1,w}$; these and $C_{SL}$ are fixed metadata unless calibration
datasets are also included. Constrain the model to positive water $T_1$,
inverse relaxivity, and sample $T_1$ over the fitted power range.

The source magnetization model is $M_0-M_0 A\,s(p)/R_1(p)$, with
$s(p)=p/(p+p_{1/2})$. Its conversion of $A$ to concentration-normalized
$k_\sigma$ is $k_\sigma=A|\omega_H/\omega_e|/C_{SL}$. Equivalently,

$$
s(p)=\frac{p}{p+p_{1/2}},
\qquad
M_\epsilon(p)=M_\epsilon(0)\left[
1-C_{SL}k_\sigma s(p)
\left|\frac{\omega_e}{\omega_H}\right|T_1(p)
\right].
$$

The explicit $C_{SL}$ preserves the source's concentration-normalized
$k_\sigma$ convention; absorbing it into $k_\sigma$ gives the sketch's
prefactor. The frequency ratio is positive and dimensionless, with both
frequencies expressed in the same units. This uses the source's saturation
law with unit limiting saturation; an unknown limiting saturation would be
absorbed into the fitted effective $k_\sigma$.

For the recovery datasets, use the FIR expression from
`~/git_repos/proc_scripts/examples/proc_FIR.py`, including its finite
repetition-delay correction:

$$
M_{\tau,j}(\tau;p_j,W_j)=M_{\tau,j}(\infty)
\left[1-\left(2-e^{-W_j/T_1(p_j)}\right)e^{-\tau/T_1(p_j)}\right].
$$

$$
W_j=10^{-6}\,\mathrm{FIR\_rep}_j
    +10^{-3}\,\mathrm{acq\_time\_ms}_j.
$$

$W_j$ is in seconds: the source stores `FIR_rep` in microseconds and
`acq_time_ms` in milliseconds. Use seconds for $\tau$ and $T_1$ too.
Parentheses denote evaluations: $M_\epsilon(0)$ is the magnetization at
$p=0$, while $M_{\tau,j}(\infty)$ is the asymptote as $\tau\to\infty$.
These endpoint values are scalar fitting parameters, represented below by
`M_p0` and `M_inf_i`, respectively. They are also the corresponding model
values at those coordinates, not additional independent variables.
In this FIR model the initial value is constrained:

$$
M_{\tau,j}(0)=M_{\tau,j}(\infty)
\left[e^{-W_j/T_1(p_j)}-1\right].
$$

Thus it is not a separate free parameter. The general endpoint form is
$M_{\tau,j}(\infty)+[M_{\tau,j}(0)-M_{\tau,j}(\infty)]
e^{-\tau/T_1(p_j)}$; the **difference** ensures the correct value at
$\tau=0$.

<!-- again, just SUPER verbose.  Give the physical model up front, legibly and compactly.-->

<!-- so, the following setup is GOOD and should come FIRST, and could preceed a lot of the stuff up above, thus allowing you to preempt the dict-like syntax I was asking for above, for example-->
Example dataset setup (all datasets already loaded in a dictionary or a
dictionary-like HDF5 group, with `T1_names` identifying the recovery nodes):

```python
for name in T1_names:
    datasets[name].name(name)
    datasets[name].set_prop("p", powers[name])
    acq = datasets[name].get_prop("acq_params")
    datasets[name].set_prop(
        "W", acq["FIR_rep"] * 1e-6 + acq["acq_time_ms"] * 1e-3
    )

datasets["M"].name("M")
# Recovery datasets have a tau axis; datasets["M"] has a p axis.
# Each dataset also exposes calibrated C_SL, T_1w0, dT_1w,
# and omega_ratio = abs(omega_e / omega_H) as top-level properties.
```

Functional forms (symbols declared with their displayed string names):

```python
functional_forms = {
    "M_tau": M_inf_i * (1 - (2 - exp(-W / T_1)) * exp(-tau / T_1)),
    "M_epsilon": M_p0 * (1 - C_SL * k_sigma * s * omega_ratio * T_1),
    "T_1": 1 / (1 / (T_1w0 + p * dT_1w) + C_SL / k_rho_inv),
    "k_rho_inv": a + b * p + c * p**2,
    "s": p / (p + p_half),
}
fit = lmfitdata(datasets, functional_forms=functional_forms)
```

During construction:

- Recovery datasets use `tau` from their axis and scalar `p` and `W` from top-level properties.
- `T_1`, `k_rho_inv`, and `s` are recursively substituted; no `T_1_i` symbols are created.
- `M_inf_i` becomes an independent localpar for each recovery dataset.
- `M_epsilon` uses the `p` axis of `datasets["M"]` and the same substituted $T_1(p)$.
- `a`, `b`, `k_sigma`, `p_half`, and `M_p0` are globalpars; fix `c=0` for the source's linear model or fit it for the quadratic extension.
- The residual concatenates all recovery residuals plus the magnetization residual.

For example, the fully substituted expression for recovery dataset `T1_0`
is, with $p=p_0$ and $W=W_0$ fixed by its properties,

$$
\begin{aligned}
M_{\tau,0}(\tau)=M_{\infty,0}\Biggl[1
&-\left(2-\exp\left\{-W_0\left[
\frac{1}{T_{1,w}(0)+p_0\Delta T_{1,w}}
+\frac{C_{SL}}{a+bp_0+cp_0^2}\right]\right\}\right) \\
&\quad\times\exp\left\{-\tau\left[
\frac{1}{T_{1,w}(0)+p_0\Delta T_{1,w}}
+\frac{C_{SL}}{a+bp_0+cp_0^2}\right]\right\}
\Biggr].
\end{aligned}
$$

The magnetization expression after substitution is

$$
M_\epsilon(p)=M_\epsilon(0)\left[
1-\frac{C_{SL}k_\sigma\left|\omega_e/\omega_H\right|p}{p+p_{1/2}}
\left(\frac{1}{T_{1,w}(0)+p\Delta T_{1,w}}
+\frac{C_{SL}}{a+bp+cp^2}\right)^{-1}
\right].
$$

The constructor summary must show each recovery dataset's own renamed
`M_inf` parameter and substituted expression, as well as the magnetization
expression. Shared relaxation parameters are optimized jointly from both
kinds of data, rather than fixed from a preceding fit.

# Stop points

## Stop point 1: container scaffold, docstring, sorting routine

Implement:

- internal container structure;
- exact preservation of single-dataset behavior;
- container-level `.copy()` and `.C`;
- decorator-preserving per-dataset transform storage/access;
- multi-dataset naming checks;
- `functional_forms` dictionary normalization, including old API to `{"z_data": expression}`;
- the variable/globalpar/localpar sorting routine;
- class docstring nomenclature and architecture TODO notes.

This stage includes code that supports multiple datasets in principle,
particularly the symbol-sorting routine; it is not limited to an empty
placeholder. Explicitly defer testing multi-dataset behavior to stop point 3.
Run existing relevant tests/lint and stop for review.

## Stop point 2: single-dataset multi-D fitting

Implement tested 2D/multi-D symbolic fitting for one dataset.

All model/data residuals should be raveled into one 1D real vector, and complex data should still use the current real-view convention. Reuse the existing compact lambda/evaluator path as much as possible.

Run tests/lint and stop for review.

## Stop point 3: tested multi-dataset fitting

Implement and test:

- functional-form dictionary assignment;
- recursive dependent substitution for `T_1`, `k_rho_inv`, and `s`;
- independent localpar `_i` expansion;
- top-level prop variables via `.set_prop("p", p_value)`;
- per-dataset decorator transforms;
- residual concatenation;
- multi-dataset `eval()` dictionary output;
- stdout construction summary.

Run tests/lint and stop for review.

## Stop point 4: synthetic ODNP test and plotted example

Add:

- a noisy synthetic ODNP test fitting FIR recovery and magnetization versus power together to recover `k_sigma`;
- a matching example that plots synthetic data and fitted curves.

This becomes the end-to-end demonstration of the intended global fitting architecture.

# Historical context and recovery commentary

The July 10 investigation reported no working general 2D fitting implementation
in the code paths and branches it inspected. Its nonsquare 2D probe failed when
raw axis arrays of different lengths reached `fitfunc_multiarg_v2` without
reshaping for broadcasting. Stop point 2 therefore includes correct axis
broadcasting as well as residual packing. This is a historical finding from
the conversation, not a fresh audit of today's repository.

The following implementation questions are recovery commentary added on
October 7, rather than additional requirements agreed in the July conversation.
They identify details that the final plan does not specify completely:

<!-- don't call the following "automatic form matching" that's not what it is!!! It's automatic identification of variables.  That's the POINT of defining the vocab, which I said was so important.  Use the all important vocab, don't make up new terminology like "automatic form matching"-->
- **Automatic form matching.** Apply the axis-coverage rule above after
  substitution: recovery has a `tau` axis and scalar `p` and `W` properties,
  while magnetization has a `p` axis. The magnetization model cannot describe
  a `tau` axis. Other ambiguous matches must raise an error.
- **Selected dataset interface.** Dictionary access must support both delegated
  `nddata` operations and the per-dataset transform decorators while retaining
  the requested `lmfitdata`/`nddata` architecture. The exact mechanism was
  not settled in the final plan.

# Building the PDF

From the repository root, using the TeX Live path configured by interactive
Bash on this machine:

```bash
source /home/jmfranck/base/bin/activate
export PATH="/usr/local/texlive/2024/bin/x86_64-linux:$PATH"
pandoc summaries/global_fitting_plan.md \
  --from=markdown --standalone --pdf-engine=pdflatex \
  --output=summaries/global_fitting_plan.pdf
```

On a machine where `pdflatex` is already on `PATH`, the `export PATH` line
is unnecessary. The document metadata sets the page margins, table of contents,
section numbering, and font size.
