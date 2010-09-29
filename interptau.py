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
def calcrho(omegaH,tau):
	#omega = 657.*omegaH
	omega = omegaH/1.5171e-3
	#sigma = 6.*J(omega+omegaH,tau)-J(omega-omegaH,tau)
	sigma = 6.*J(omega-omegaH,tau)-J(omega+omegaH,tau)
	#denom = 6.*J(omega+omegaH,tau)+3.*J(omegaH,tau)+J(omega-omegaH,tau)
	denom = 6.*J(omega-omegaH,tau)+3.*J(omegaH,tau)+J(omega+omegaH,tau)
	rho = sigma/denom
	return rho
def interptau(rhoinput,nmr_freq):
    tau = logspace(-15,-9,1000)
    B0 = double(nmr_freq)*1e6/4.258e7
    omegaH = 2*pi*B0*4.258e7
    rho = calcrho(omegaH,tau)
    order = argsort(rho)
    interpresult = interp(double(rhoinput),rho[order],tau[order])
    return interpresult,tau,rho
