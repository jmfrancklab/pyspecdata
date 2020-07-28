# -*- coding: utf-8 -*-
from pyspecdata import *
#from scipy.signal import gaussian
from scipy.special import jv
from scipy.interpolate import UnivariateSpline
from scipy.signal import tukey
from matplotlib import rcParams
init_logging('debug')

def preproc_ESR_multiharmonic_centered(s):
    """load a multiharmonic spectrum, throwing out the modulation quadrature,
    and reset the field axis so the spectrum is centered at 0"""
    s.chunk_auto('harmonic','phase')
    s = s['phase',0]
    s -= s['$B_0$',:50].C.mean('$B_0$')
    # {{{ determine the center automatically
    print(s['harmonic',0].C.run_nopop(cumsum,'$B_0$'))
    s_integral = s['harmonic',0].C.run_nopop(cumsum,'$B_0$')
    x1,x2 = s_integral.getaxis('$B_0$')[r_[5,-5]]
    y1 = s_integral.data[:5].mean()
    y2 = s_integral.data[-5:].mean()
    straight_baseline = (s.fromaxis('$B_0$')-x1)*(y2-y1)/(x2-x1)
    s_integral -= straight_baseline
    s_integral /= s_integral.data.mean()
    center_field = (s_integral * s.fromaxis('$B_0$')).mean('$B_0$').item()
    # }}}
    s.setaxis('$B_0$',lambda x: x-center_field) 
    return s
s = find_file('15N_S175R1a_pR_new_9G_200309.DSC',
        exp_type='francklab_esr/Sam',
        postproc=preproc_ESR_multiharmonic_centered)
man_calc = []
for j in range(ndshape(s)['harmonic']):
    man_calc.append((s['harmonic',j].data**2).sum())
result = s.along('$B_0$') @ s
print("find the norm squared", result)
assert all(isclose(array(man_calc), result.data))
