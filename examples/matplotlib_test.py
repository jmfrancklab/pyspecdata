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
from matplotlib.cbook import flatten, silent_list, iterable, dedent
import matplotlib as mpl
from matplotlib.dates import (
    date2num, num2date, datestr2num, strpdate2num, drange, epoch2num,
    num2epoch, mx2num, DateFormatter, IndexDateFormatter, DateLocator,
    RRuleLocator, YearLocator, MonthLocator, WeekdayLocator, DayLocator,
    HourLocator, MinuteLocator, SecondLocator, rrule, MO, TU, WE, TH, FR,
    SA, SU, YEARLY, MONTHLY, WEEKLY, DAILY, HOURLY, MINUTELY, SECONDLY,
    relativedelta)
# bring all the symbols in so folks can import them from
# pylab in one fell swoop
## We are still importing too many things from mlab; more cleanup is needed.
from matplotlib.mlab import (
    demean, detrend, detrend_linear, detrend_mean, detrend_none,
    window_hanning, window_none)
from matplotlib import cbook, mlab, pyplot as plt
from matplotlib.pyplot import *
from numpy import *
from numpy.fft import *
from numpy.random import *
from numpy.linalg import *
import numpy as np
import numpy.ma as ma
# don't let numpy's datetime hide stdlib
import datetime
# This is needed, or bytes will be numpy.random.bytes from
# "from numpy.random import *" above
bytes = __import__("builtins").bytes
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
