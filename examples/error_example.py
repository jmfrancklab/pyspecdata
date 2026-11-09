"""
Error and units example
=======================

Here is a simple example of errors and unit propagation

Notice that the base nddata class supplies error and propagation similar
to uncertainty-type libraries.

For the units, `pint` is doing the heavy lifting here.
"""

import pyspecdata as psd
from numpy import ones, r_
import matplotlib.pyplot as plt

# %%
# As a simple example, say that we perform several measurements of a
# volume (not sure physically why we would have such variability, but
# let's roll with it to keep the example simple!)

vol = psd.nddata(r_[1.10, 1.11, 1.02, 1.03, 1.00, 1.05]).set_units("L")

# %%
# Similarly, let's take some measurements of the weight of a solute!

weight = psd.nddata(r_[2.10, 2.61, 2.002, 2.73, 2.33, 2.69]).set_units("μg")

# %%
# To test our error propagation below, we're going to divide the two
# arrays here -- because the variability of this number should be somewhat
# similar to the propagated error below (though of course, there is a
# statistical difference, and doing the two things does mean something
# different).
# Notice how, during string conversion, we always give the standard error
# 2 significant figures, and then base the significant figures of the
# number on the error.

conc_indiv = weight / vol
conc_indiv.mean(stderr=True)
print(conc_indiv)

# %%
# We take the mean, allowing it to accumulate the standard error.
# (See comment above about significant figures.)

vol.mean(stderr=True)
print(vol)

weight.mean(stderr=True)
print(weight)

print(weight / vol)

# %%
# Now, suppose that we know the uncertainty of each individual weighing
# (the balance was steadier for some than for others), that the volumes
# are known to ±0.02 L, and that we note when each sample was drawn.
# Dividing propagates the error into every individual concentration.

conc_each = (
    psd.nddata(r_[2.10, 2.61, 2.002, 2.73, 2.33, 2.69], "t")
    .setaxis("t", r_[0:30:6j])
    .set_units("t", "min")
    .set_units("μg")
    .set_error(r_[0.05, 0.40, 0.05, 0.60, 0.10, 0.30])
) / (
    psd.nddata(r_[1.10, 1.11, 1.02, 1.03, 1.00, 1.05], "t")
    .setaxis("t", r_[0:30:6j])
    .set_units("t", "min")
    .set_units("L")
    .set_error(0.02 * ones(6))
)

# %%
# Two of the samples were drawn at times we only know roughly, so the time
# coordinate also has an error.

conc_each.set_error("t", r_[0.5, 0.5, 3, 0.5, 4, 0.5])

# %%
# Error bars draw the eye to the *least* certain points, since they have
# the longest bars.  ``plot_spreading_box`` instead uses the same amount of
# pigment for each point, spread over a box whose height and width follow
# the error of the data and of the time coordinate: the certain points stay
# crisp, and the uncertain ones fade out.

fig, (ax_bars, ax_boxes) = plt.subplots(1, 2, sharey=True)
psd.plot(conc_each, "o", ax=ax_bars)
conc_each.plot_spreading_box(ax=ax_boxes)

# Because we want this to show up in sphinx gallery, we have
# to make some type of figure

fig = plt.figure()
fig.add_artist(
    plt.Text(
        x=0.5,
        y=0.5,
        text=(
            "Compare one-by-one division,\nfollowed by standard"
            f" error:\n{conc_indiv}\nto division with error"
            f" propagation:\n{weight/vol}"
        ),
        fontsize=40,
        ha="center",
        va="center",
    )
)
fig.tight_layout()
plt.show()
