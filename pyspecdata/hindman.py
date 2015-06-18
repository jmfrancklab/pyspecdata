from .h5nmr import *
def hindmant1(fitT,f,oxygen = 70e-3):
    'output the hindman model for temperatures fitT at nmr frequency f'
    expcor = r_[0.4554e-8 , 0.4787e4 , 0.6343e-3 , 0.1764e4]
    theird = r_[0.1074e-7,0.3063e-8,0.5954e-5,0.4224e4,0.4726e4,0.4612e4,0.6494e-3,0.3420e-3,0.4466,0.1794e4,0.1817e4,0.1650e4].reshape(4,3)

    #print lsafen('their data:',theird[:,1])
    #R1_corr = theird[0,1] * exp(theird[1,1]/fitT)
    #R1_corr += theird[2,1] * exp(theird[3,1]/fitT)
    bpp = lambda x,f: (x/(1.+(2*pi*f*x)**2)+4.*x/(1.+4*(2*pi*f*x)**2))
    lorentzian = lambda x,f: (x/(1.+(pi*f*x)**2))

    R1_corr = expcor[0] * exp(expcor[1]/fitT)
    R1_corr += expcor[2] * exp(expcor[3]/fitT)
    tau_r = 1.23e-13
    R1_corr /= bpp(tau_r,60e6)
    R1_corr *= bpp(tau_r,f)

    TrI = 1e-40*(2.9376 + 1.9187 + 1.0220)
    TrI *= (1e-3)*(1e-4) # g cm^2
    TrC2 = 4*pi**2*((-32.68e3)**2 + (-33.474e3)**2 + (-32.352e3)**2) # Hz^2
    h = 6.626068e-34
    kB = 1.3806503e-23
    tau_sr = 1.23e-13
    srRate = 2. / 9. * fitT * kB * (h**-2) * TrI * TrC2 * lorentzian(tau_sr,f)
    retval = nddata(1/(R1_corr+srRate+oxygen),[len(R1_corr)],['T'])
    retval.labels(['T'],[fitT-273])
    return retval
#thisfieldcurve = hindmant1(linspace(10.,65.100)+273,15e6)
#c,lineardata = thisfieldcurve.polyfit('T')
#propconst = c[0]/c[1] + 25.
#def magtotemp(val,ini_temp = 25.):
#    'an older form, based on the linear equation'
#    return ini_temp-val*propconst
def genTconst(ini_temp = 25.,rel_tau = 1.,Trange = r_[25.,40.]):
    'generate c for a specific initial temperature'
    rangeofT = linspace(Trange[0],Trange[1],100)
    t1p = hindmant1(rangeofT+273,15e6).data
    t1 = hindmant1(r_[ini_temp]+273,15e6).data
    m0point = 1.-exp(-rel_tau/t1)
    mtrange = 1.-exp(-rel_tau/t1p)
    mratio = nddata(rangeofT,[len(mtrange)],['Mratio']) # the x axis is the Mratio!
    mratio.labels(['Mratio'],[m0point/mtrange-1.])
    c,mratio_lin = mratio.polyfit('Mratio')
    return c
def magtoT(val,c):
    return c[0]+c[1]*val
