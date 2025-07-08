"""
2D ILT
======
Here, we're going to provide a few demonstrations of the ILT functionality.
Let's start with Fig 1.10 in A. Beaton's thesis, which is based off the
figures in Venkataramanan.
"""

from pylab import (
    figure,
    title,
    show,
    linspace,
    logspace,
    log10,
    exp,
    sqrt,
    rcParams,
)
import numpy as np
from pyspecdata import nddata, image

rcParams["image.aspect"] = "auto"  # needed for sphinx gallery
# sphinx_gallery_thumbnail_number = 2

NT1 = 300  # Number of T1 values
NT2 = 300  # Number of T2 values
LT1_name = r"$\log(T_1)$"
LT1 = nddata(linspace(-2.5, 0.5, NT1), LT1_name)
LT2_name = r"$\log(T_2)$"
LT2 = nddata(linspace(-2.5, 0.3, NT2), LT2_name)
mu = [-1.25, -1.75]
sigma = [0.1, 0.1]
exact_data = exp(
    -((LT1 - mu[0]) ** 2) / 2 / sigma[0] ** 2
    - (LT2 - mu[1]) ** 2 / 2 / sigma[1] ** 2
)
slanted_coord1 = (LT1 + LT2) / sqrt(2)
slanted_coord2 = (LT2 - LT1) / sqrt(2)
mu = [-1.0, -0.4]
mu = [  # convert to slanted coords
    (mu[0] + mu[1]) / sqrt(2),
    (mu[1] - mu[0]) / sqrt(2),
]
sigma = [0.5, 0.05]  # in slanted
exact_data += exp(
    -((slanted_coord1 - mu[0]) ** 2) / 2 / sigma[0] ** 2
    - (slanted_coord2 - mu[1]) ** 2 / 2 / sigma[1] ** 2
)
exact_data.reorder(LT2_name)  # T₂ along y axis

figure(1)
title("exact data")
image(exact_data)

# Now add the experimental decay dimensions

tau1 = nddata(logspace(log10(5.0e-4), log10(4), 30), "tau1")
tau2 = nddata(linspace(5.0e-4, 3.8, 1000), "tau2")

# Pre-allocate the τ₁×τ₂ result via ndshape’s alloc, inline

simulated_data = (tau1.shape | tau2.shape).alloc(dtype=np.float64)

# %%
# pySpecData makes it easy to construct fake data like this. Typically this is
# very easy, but here, we must contend with the fact that we are
# memory-limited, so if we want a highly resolved fit basis, we need to chunk
# up the calculation.  Nonetheless, pySpecData still makes that easy: let's see
# how!
#
# Block sizes (tune to available RAM)

bLT1 = 20
bLT2 = 20

# %%
# Loop over LT1 and LT2 in blocks, vectorized over τ dims each time
#
# $$T_1 = 10^{\log(T_1)}$$
# $$R_1 = 10^{-\log(T_1)}$$
# $$\ln(R_1) = -\log(T_1) \ln(10)$$

print(
    "Generating the fake data can take some time.  I need to loop a"
    f" calculation in chunks over a {LT1.shape[LT1_name]} ×"
    f" {LT2.shape[LT2_name]} grid"
)
for i in range(0, LT1.shape[LT1_name], bLT1):
    LT1_blk = LT1[LT1_name, slice(i, i + bLT1)]
    B1 = 1 - 2 * exp(-tau1 / 10**LT1_blk)  # dims: (tau1, LT1_blk)
    for j in range(0, LT2.shape[LT2_name], bLT2):
        print(i, j)
        LT2_blk = LT2[LT2_name, slice(j, j + bLT2)]
        B2 = exp(-tau2 / 10**LT2_blk)  # dims: (tau2, LT2_blk)
        # Extract matching block of exact_data
        data_blk = exact_data[LT1_name, slice(i, i + bLT1)][
            LT2_name, slice(j, j + bLT2)
        ]  # dims: (tau1, tau2, LT1_blk, LT2_blk)
        # Multiply, sum out both LT axes, and accumulate
        simulated_data += (B2 * B1 * data_blk).real.sum(LT1_name).sum(LT2_name)
print("done generating")

# `simulated_data` now holds the τ₁×τ₂ synthetic data

# then add noise, and scale data so that noise has norm of 1

simulated_data.add_noise(0.1)
simulated_data /= 0.1


# %%
# Use BRD to find the value of $\lambda$ ($\alpha$).
# Note that BRD assumes that you have scaled your data so that the stdev
# of the noise is 1.0.

simulated_data.nnls(
    ("tau1", "tau2"),
    (LT1, LT2),
    (
        lambda tau1, LT1: 1 - 2 * exp(-tau1 * 10**-LT1),
        lambda tau2, LT2: exp(-tau2 * 10**-LT2),
    ),
    l="BRD",
)

figure(2)
title("BRD")
image(simulated_data)

show()
