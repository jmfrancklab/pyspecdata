from pylab import *
from pyspecdata import *

# Here, we're going to provide a few demonstrations of the ILT functionality.  Let's start with Fig 1.10 in A. Beaton's thesis, which is based off the figures
# in Venkataramanan.
# 
# pySpecData makes it easy to construct fake data like this: let's see how!

NT1 = 200 # Number of T1 values
NT2 = 200 # Number of T2 values
LT1_name = r'$\log(T_1)$'
LT1 = nddata(linspace(-2.5,0.5,NT1), LT1_name)
LT2_name = r'$\log(T_2)$'
LT2 = nddata(linspace(-2.5,0.3,NT2), LT2_name)
mu = [-1.25,-1.75]
sigma = [0.1,0.1]
exact_data = exp(-(LT1-mu[0])**2/2/sigma[0]**2-(LT2-mu[1])**2/2/sigma[1]**2)
slanted_coord1 = (LT1+LT2)/sqrt(2)
slanted_coord2 = (LT2-LT1)/sqrt(2)
mu = [-1.0,-0.4]
mu = [ # convert to slanted coords
        (mu[0]+mu[1])/sqrt(2),
        (mu[1]-mu[0])/sqrt(2)]
sigma = [0.5,0.05] # in slanted
exact_data += exp(-(slanted_coord1-mu[0])**2/2/sigma[0]**2-(slanted_coord2-mu[1])**2/2/sigma[1]**2)
exact_data.reorder(LT2_name) # T₂ along y axis

figure(1)
title("exact data")
image(exact_data)

# Now add the experimental decay dimensions

tau1 = nddata(logspace(log10(5.0e-4),log10(4),30),'tau1')
tau2 = nddata(linspace(5.0e-4,3.8,1000),'tau2')

# $$T_1 = 10^{\log(T_1)}$$
# $$R_1 = 10^{-\log(T_1)}$$
# $$\ln(R_1) = -\log(T_1) \ln(10)$$

basis = exp(-tau2/10**LT2)*(1-2*exp(-tau1/10**LT1))

# Sum along the distribution dimensions, to create fake data
# then add noise, and scale data so that noise has norm of 1

simulated_data = basis*exact_data
simulated_data.sum(LT1_name).sum(LT2_name)
simulated_data.add_noise(0.1)
simulated_data /= 0.1


# use BRD to find the value of $\lambda$ 

simulated_data.nnls(('tau1','tau2'),(LT1,LT2),
                    (lambda tau1,LT1: 1-2*exp(-tau1*10**-LT1),
                     lambda tau2,LT2: exp(-tau2*10**-LT2)), l='BRD')

figure(2)
title("BRD")
image(simulated_data)

show()
