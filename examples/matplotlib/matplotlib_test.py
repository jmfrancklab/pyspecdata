"""Matplotlib Example
==================

this is a basic example that should generate images in
sphinx, but still imports pyspecdata"""
from pyspecdata import *
from matplotlib.patches import Ellipse
from numpy import *
from matplotlib.pyplot import subplot, xlim, ylim

delta = 45.0 # degrees

angles = arange(0, 360+delta, delta)
ells = [Ellipse((1, 1), 4, 2, a) for a in angles]

a = subplot(111, aspect='equal')

for e in ells:
    e.set_clip_box(a.bbox)
    e.set_alpha(0.1)
    a.add_artist(e)

xlim(-2, 4)
ylim(-1, 3)

figure()

b = subplot(111,aspect='equal')

plot(r_[0:10])

show()
