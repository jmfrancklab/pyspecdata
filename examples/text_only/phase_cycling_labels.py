"""
phase cycling labels
====================

demonstrate automatic labeling of phase cycling (coherence) axes
"""
from pyspecdata import *
a = nddata(r_[0:90],[3,3,10],['ph1','ph2','t']).setaxis('ph1','#').setaxis('ph2','#').setaxis('t','#').set_units('t','s')
a.ft(['ph1','ph2','t'])
print(a.unitify_axis('ph1'))
print(a.unitify_axis('ph2'))
print(a.unitify_axis('t'))
