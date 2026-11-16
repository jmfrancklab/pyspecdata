# Writing pyspecdata-native analysis

Use this as shared guidance for analysis scripts, teaching examples, and
notebook code using `pyspecdata.nddata`. Keep the data, named axes, units,
errors, and transform state together throughout the calculation.

## Native operations first

- Before accessing `.data`, constructing an `nddata` from a raw NumPy array,
  or copying an `nddata`, stop and check again for a native approach.
  Look for named-axis slicing, arithmetic and broadcasting by axis name,
  reductions, `argmin`/`argmax`, `contiguous`, transforms, `fromaxis`, labels,
  units, and properties. Consult the installed implementation and examples
  rather than guessing an API.
- Avoid extracting arrays, operating on their positional dimensions, and
  wrapping the result in a new `nddata`. That needlessly reconstructs metadata
  and can lose units, errors, or Fourier state.
- Use axis names rather than positional axis numbers for processing. Let
  native arithmetic align named dimensions instead of manually reshaping or
  broadcasting arrays.
- Use `getaxis`, `setaxis`/`set_axis`, `set_units`, `name`, and properties when
  the task concerns coordinates or metadata. Do not drill into `.data` to
  recover information already available through these interfaces.
- Native reductions and peak/region selection should come before NumPy
  reductions, manual index searches, or masks built from extracted data.
  Check whether `argmin`/`argmax` returns coordinates or indices for the call
  you intend; do not assume NumPy semantics.

## Two kinds of named-axis slicing

Choose the syntax according to what the selection means. The distinction is
**comma versus colon**, not simply integer versus floating-point arguments.

| Intent | Syntax | Meaning |
| --- | --- | --- |
| One sample by index | `d["t", 5]` | Sixth sample along `t` |
| Samples by index range | `d["t", 5:20]` | Indices 5 through 19 |
| Every other sample | `d["t", ::2]` | Stride through sample indices |
| One point by coordinate | `d["t":1.1]` | Sample nearest coordinate 1.1 |
| Coordinate interval | `d["t":(0, 100e-6)]` | Inclusive coordinate range |
| Open coordinate interval | `d["t":(None, 100e-6)]` | Through that coordinate |

- Use coordinate slicing for physical regions, integration windows, and
  spectral bands. Values use the stored axis units; setting units does not
  make an arbitrary display prefix part of the slicing syntax.
- Use index slicing for sample counts, strides, and selection of a particular
  acquired record. An index stride subsamples existing points; it does not
  interpolate or apply an antialiasing filter.
- After `d.ft("t", shift=True)`, `d["t":(-30e3, 30e3)]` selects frequency
  coordinates on the transformed axis. The axis name remains `t`; the FT
  updates its coordinates, units, and transform state. Do not manually rebuild
  a frequency array or rename the dimension just to label the plot.
- For conditions on coordinates, consider a named-axis callable selection
  such as `d["t", lambda x: abs(x) < 1e-3]`. A bare callable selection acts on
  data values instead. Verify the relevant multidimensional behavior before
  using either, and consider `contiguous` for connected regions.
- A native slice preserves the relevant labeled-data structure. Do not add
  `.C` reflexively before slicing; check the slice's ownership/mutation
  behavior if later operations must keep the source intact.

## One evolving dataset, clear operation order

- Keep one variable name as a dataset moves through processing. Prefer
  `d.ft("t", shift=True)` over `spectrum = d.C.ft("t", shift=True)` when the
  time-domain object is no longer needed. Plot or consume the earlier state
  before mutating it.
- Chain related native operations when their return values support it and
  the expression remains readable. For example,
  `d.interp("t", new_grid).ft("t", shift=True)` describes one processing flow.
- Before `.C`, `.copy()`, or an equivalent copy, ask whether reordering the
  work removes the need to retain both states. Assignment (`other = d`)
  creates an alias, not an independent processing branch.
- When two independent branches are genuinely needed, use the minimum
  number of copies and explain why. `.C` is a deep copy. For two resampling
  branches, first plot the dense reference, then use
  `fast = fid.C.interp("t", fast_grid)` and `fid.interp("t", slow_grid)`.
  One branch preserves an independent state; the original supplies the other.
- Copies made internally by native methods are the library's responsibility.
  The caller should avoid unnecessary explicit copies, not bypass the native
  API to eliminate its internal implementation details.

## Generate functions from labeled coordinates

- Use `fromaxis` when a function is defined on an existing dataset's axes:
  `rect = t.fromaxis("t", lambda x: 1.0 * (abs(x) < 0.8))`.
  Use native arithmetic for expressions such as `exp(-(t / 0.5) ** 2)`.
- When appropriate, `fromaxis(..., overwrite=True)` lets the same object
  evolve from coordinates to a derived signal. Check whether those coordinates
  are needed for another signal before overwriting.
- A constructor at the acquisition/import boundary, or to seed a genuinely
  new synthetic coordinate grid, can be necessary. Document that boundary
  and attach axis names and units there. It does not justify extracting and
  rewrapping arrays later in the calculation.
- If a third-party algorithm truly requires NumPy arrays, keep that boundary
  narrow and document why the native API is inadequate and how relevant
  metadata is preserved.

## Compact plotting and teaching examples

- Pass `nddata` directly to pyspecdata's `plot`, or use the established
  `figlist_var` workflow. Let named axes, units, and Fourier metadata supply
  axis labels. Use `name` for a meaningful dataset label when appropriate.
- Omit redundant `xlabel`, `ylabel`, `xlim`, and `ylim` calls. Use automatic
  plot limits. If a scientific question concerns a specific region, select
  that region with native coordinate slicing and explain the selection.
- Prefer native `.real`, `.imag`, and `abs(d)` over manipulating `.data`.
- Keep example cells focused on data creation, processing, and plotting.
  Use explicit imports so students can see which operations come from NumPy,
  Matplotlib, and pyspecdata; do not use wildcard imports.
  Put shared imports and defaults such as `rcParams["figure.figsize"]` in
  one initial notebook cell, rather than repeating them in each activity.
  Put shared grids and minor-tick visibility in that same `rcParams` setup.
  Do not silently add plot styling in the export code when the displayed
  cells are intended to reproduce the slide figures.
  Shared figure styling, layout, export, and slide formatting belong in the
  surrounding builder when they are not part of the lesson.
- Use readable chains and a sensible formatter line length. Do not turn
  short calls into many lines merely to satisfy an unusually narrow setting.
  Remove redundant operations before reducing font size to fit a slide.
- Execute the exact displayed code, including successive cells in the same
  namespace, to verify that examples do what they claim.

## Notebook context

Notebook-tree `.qmd` files run through `pydifft qmd`/`qmdb`, never Quarto.
Follow nearby project notebooks for their custom formatting. Strip blank
lines from Python code blocks so gvim's `}` motion passes over the block.
Keep notebook-specific recording and `<obs>` rules in each notebook's own
instructions rather than imposing them on ordinary Python scripts.

## Sources and scope

This guide is shared guidance for programs and projects using pyspecdata.
Repository setup belongs in each project's `AGENTS.md`; notebook recording
rules belong in the relevant notebook's instructions.

Consolidated on October 13, 2026 from the author's notebook instructions and
explicit guidance on compact teaching examples. API behavior was checked
against the local implementation and documentation:

- [Named-axis operations and slicing](doc_src/nddata.rst).
- [Implementation](pyspecdata/core.py): slicing, `fromaxis`, plotting labels,
  `interp`, `name`, and deep copying with `.C`.
