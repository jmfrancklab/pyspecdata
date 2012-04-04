from pylab import *
import sys
def J(omega,tau):
	st = sqrt(2)
	sot = sqrt(omega*tau)
	J = 1.0+st*sot+sot**2
	#J += st/3.*(sot**3) + 16./81.*((omega*tau)**2) +4.*st/81.*(sot**5) + (sot**3)/81.
	J += st/3.*(sot**3) + 16./81.*(sot**4) +4.*st/81.*(sot**5) + (sot**6)/81.

	J = (1.+5.*st/8.*sot+(sot**2)/4.)/J
	return J
def calcxi(omegaH,tau):
	#omega = 657.*omegaH
	omega = omegaH/1.5171e-3
	#sigma = 6.*J(omega+omegaH,tau)-J(omega-omegaH,tau)
	sigma = 6.*J(omega-omegaH,tau)-J(omega+omegaH,tau)
	#denom = 6.*J(omega+omegaH,tau)+3.*J(omegaH,tau)+J(omega-omegaH,tau)
	denom = 6.*J(omega-omegaH,tau)+3.*J(omegaH,tau)+J(omega+omegaH,tau)
	xi = sigma/denom
	return xi
def calcrho(omegaH,tau):
	omega = omegaH/1.5171e-3
	denom = 6.*J(omega-omegaH,tau)+3.*J(omegaH,tau)+J(omega+omegaH,tau)
	return denom
def interptau(xiinput,nmr_freq,simple = False,tau_short = 10**(-15),tau_long = 10**(-8)):
    tau = logspace(log10(tau_short),log10(tau_long),1000)
    B0 = double(nmr_freq)*1e6/4.258e7
    omegaH = 2*pi*B0*4.258e7
    xi = calcxi(omegaH,tau)
    order = argsort(xi)
    interpresult = interp(double(xiinput),xi[order],tau[order])
    if simple:
        return interpresult
    else:
        return interpresult,tau,xi
def interpxi(tauinput,nmr_freq,simple = False,tau_short = 10**(-15),tau_long = 10**(-8)):
    tau = logspace(log10(tau_short),log10(tau_long),1000)
    B0 = double(nmr_freq)*1e6/4.258e7
    omegaH = 2*pi*B0*4.258e7
    xi = calcxi(omegaH,tau)
    order = argsort(tau)
    interpresult = interp(double(tauinput),tau[order],xi[order])
    if simple:
        return interpresult
    else:
        return tauinput,interpresult,xi
