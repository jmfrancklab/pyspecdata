"""
============
Artist tests
============

Test unit support with each of the Matplotlib primitive artist types.

The axis handles unit conversions and the artists keep a pointer to their axis
parent. You must initialize the artists with the axis instance if you want to
use them with unit data, or else they will not know how to convert the units
to scalars.

.. only:: builder_html

   This example requires :download:`basic_units.py <basic_units.py>`
"""
import random
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.text as text
import matplotlib.collections as collections

#from basic_units import cm, inch
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
# test a plain-ol-line
line = lines.Line2D([0, 0.5], [0, 1],
                    lw=2, color='r',
                    transform=ax.transAxes)
ax.add_line(line)

plt.show()
