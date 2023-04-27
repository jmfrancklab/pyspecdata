from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import QESR_scalefactor
from sympy import symbols, Symbol, latex
from sympy import exp as s_exp
import numpy as np
fieldaxis = '$B_0$'
exp_type = "francklab_esr/alex"
with figlist_var() as fl:
    for filenum, (thisfile, calibration, diameter, background) in enumerate(
            [
              ('230217_46uM_TEMPOL.DSC', '230202', 'QESR caps',
                  find_file('230217_water.DSC',exp_type=exp_type)['harmonic',0]),
             ]):
        #{{{ load data in and rescale and center about 0
        d = find_file(thisfile,exp_type = exp_type)
        if "harmonic" in d.dimlabels:
            d = d['harmonic',0]
        d /= QESR_scalefactor(d, calibration_name = calibration,
                diameter_name = diameter)
        center_field = d.getaxis(fieldaxis)[r_[0,-1]].mean()
        d -= d[fieldaxis, -100:].data.mean()
        d.setaxis(fieldaxis, lambda x: x-center_field)
        background = ndshape(d).alloc()
        background.copy_axes(d)
        background = background.interp(fieldaxis, d.getaxis(fieldaxis))
        d -= background
        #}}}
        d.ift(fieldaxis, shift=True)
        d *= 2
        d[fieldaxis:0] *= 0.5
        d.ft(fieldaxis)
        #{{{ Plot starting spectra
        d.set_units(fieldaxis,'G')
        fl.next("Starting spectrum")
        fl.plot(d, label = 'actual data')
        centers = [-18,-1,16] #make a list of the centers of the starting data for the fits
        d.ift(fieldaxis)
        fl.next('u domain')
        fl.plot(d,label = 'real data')
        fl.plot(d.imag,label = 'imag data')
        #}}}
        d.rename(fieldaxis,'u')
        u = d.fromaxis('u')
        #{{{ make symbols and prep lists for symbols
        lambda_l1, lambda_l2, lambda_l3, lambda_g1, lambda_g2, lambda_g3, u, amp1, amp2, amp3, Bc1, Bc2, Bc3 = symbols("lambda_l1 lambda_l2 lambda_l3 lambda_g1 lambda_g2 lambda_g3 u amp1 amp2 amp3 Bc1 Bc2 Bc3", real=True) 
        deriv = -1j * 2 * pi * u
        lambda_L_list = [lambda_l1, lambda_l2, lambda_l3]
        lambda_G_list = [lambda_g1, lambda_g2, lambda_g3]
        Bc_list = [Bc1, Bc2, Bc3]
        amp_list = [amp1, amp2, amp3]
        #}}}
        #{{{ index and make the functions using the symbols - we will have 3 lorentzians, 3 gaussians and 3 shifts
        Lorentz_list = [0,1,2]
        Gauss_list = [0,1,2]
        shift_list = [0,1,2]
        for j in range(len(Lorentz_list)):
            Lorentz_list[j] = s_exp(-pi*lambda_L_list[j] * abs(u))
            Gauss_list[j] = s_exp(-(pi**2*lambda_G_list[j]**2*u**2)/(4*np.log(2)))
            shift_list[j] = s_exp(+1j*pi*2*u*Bc_list[j])
        #}}}    
        #{{{ make a fitting function as a sum of the 3 voigt lines using the function lists we just made
        myfit = 0
        for j in range(len(Gauss_list)):
            myfit += (amp_list[j] * Gauss_list[j] * Lorentz_list[j] * shift_list[j] * deriv)
        #}}}
        thisdata = nddata(d.data, ['u'])
        thisdata.setaxis('u',d.getaxis('u'))
        f = lmfitdata(thisdata)
        f.functional_form = myfit
        #{{{ Make lists of guesses
        lambda_l_guess = [1.1,1.1,1.1]
        lambda_g_guess = [1.1,1.1,1.1]
        amp_guess = [30,30,30]
        #}}}
        #{{{ set your guesses
        for j in range(len(shift_list)):
            f.set_guess(
                    lambda_l1 = dict(value = lambda_l_guess[0], min = 0.1, max = 3),
                    lambda_l2 = dict(value = lambda_l_guess[1], min = 0.1, max = 3),
                    lambda_l3 = dict(value = lambda_l_guess[2], min = 0.1, max = 3),
                    lambda_g1 = dict(value = lambda_g_guess[0], min = 0.1, max = 3),
                    lambda_g2 = dict(value = lambda_g_guess[1], min = 0.1, max = 3),
                    lambda_g3 = dict(value = lambda_g_guess[2], min = 0.1, max = 3),
                    amp1 = dict(value = amp_guess[0], min = 5e-6, max = 5e6),
                    amp2 = dict(value = amp_guess[1], min = 5e-6, max = 5e6),
                    amp3 = dict(value = amp_guess[2], min = 5e-6, max = 5e6),
                    Bc1 = dict(value = centers[0], min = -30, max = 30),
                    Bc2 = dict(value = centers[1], min = -30, max = 30),
                    Bc3 = dict(value = centers[2], min = -30, max = 30)
                    )
        f.settoguess()
        guess = f.eval()
        guess.ft('u',shift = True)
        guess.set_units('u','G')
        #}}}
        #{{{fit 
        f.fit()
        fit = f.eval()
        #}}}
        fit.ft('u',shift=True)
        fit.ift('u')
        fl.next('u domain')
        fl.plot(fit,label = 'fit')
        fl.next('Starting spectrum')
        fit.ft('u')
        fit.set_units('u','G')
        fl.plot(fit, label = 'fit')

