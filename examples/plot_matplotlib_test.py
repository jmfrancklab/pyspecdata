# testing for compatability of pyspecdata + matplotlib
# plotting utils
from pyspecdata import *
from pyspecdata import fitdata 
from matplotlib.patches import Ellipse
#from pylab import figure,subplot,plot,xlim,ylim,show # doesn't work
#from pylab import * # does work
# {{{ this is the contents of pylab.py -- works
# need to go through and figure out which lines
# are actually needed and which are not
# -- I have already stripped out some
from matplotlib.pyplot import figure, subplot, show, xlim, ylim, plot
from numpy import * # I think it wasn't importing from numpy b/c it seems we're inside sphinx
# }}}


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
