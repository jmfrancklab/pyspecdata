# In[1]:

get_ipython().magic(u'load_ext pyspecdata.ipy')
init_logging(level='debug')

# Here, we're going to provide a few demonstrations of the ILT functionality.  Let's start with Fig 1.10 in A. Beaton's thesis, which is based off the figures
# in Venkataramanan.
# 
# pySpecData makes it easy to construct fake data like this: let's see how!

NT1 = 50 # Number of T1 values
NT2 = 50 # Number of T2 values
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
exact_data

# In[2]:

tau1 = nddata(logspace(log10(5.0e-4),log10(4),30),'tau1')
tau2 = nddata(linspace(5.0e-4,3.8,1000),'tau2')

# $$T_1 = 10^{\log(T_1)}$$
# $$R_1 = 10^{-\log(T_1)}$$
# $$\ln(R_1) = -\log(T_1) \ln(10)$$

basis = exp(-tau2/10**LT2)*(1-2*exp(-tau1/10**LT1))
basis.shape

# In[3]:

simulated_data = basis*exact_data
simulated_data.sum(LT1_name).sum(LT2_name)

# In[4]:

simulated_data.nnls(('tau1','tau2'),(LT1,LT2),
                    (lambda tau1,LT1: 1-2*exp(-tau1*10**-LT1),
                     lambda tau2,LT2: exp(-tau2*10**-LT2)), l='BRD')

# ## 1.5D test
# Now, we use the same shape, but we pretend like one of these
# dimensions comes directly from the spectroscopy (i.e., it's a
# direct FT dimension, as in the x axis of DOSY), while the other
# requires regularization

basis = (1-2*exp(-tau1/10**LT1))
simulated_data = basis*exact_data
simulated_data.sum(LT1_name)

# now perform the regularization

simulated_data.nnls('tau1', LT1,
                    lambda tau1,LT1: 1-2*exp(-tau1*10**-LT1),
                    l='BRD')

