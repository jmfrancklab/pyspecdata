# In[1]:

get_ipython().magic(u'load_ext pyspecdata.ipy')

# Here, we're going to provide a few demonstrations of the ILT functionality.  Let's start with Fig 1.10 in A. Beaton's thesis, which is based off the figures
# in Venkataramanan.
# 
# pySpecData makes it easy to construct fake data like this: let's see how!

LT1_name = r'$\log(T_1)$'
LT1 = nddata(r_[-2.5:0.25:40j], LT1_name)
LT2_name = r'$\log(T_2)$'
LT2 = nddata(r_[-2.5:0.25:40j], LT2_name)
exact_data = exp(-(LT1+1.25)**2/2/0.1**2-(LT2+1.25)**2/2/0.1**2)
slanted_coord1 = LT1+LT2
slanted_coord2 = LT1-LT2
exact_data += exp(-(slanted_coord1+1)**2/2/0.1**2-(slanted_coord2)**2/2/0.5**2)
exact_data

# In[2]:

tau1 = nddata(logspace(-2.5,0.25+log10(4),30),'tau1')
tau2 = nddata(linspace(10**-2.5,10**(0.25+log10(4)),30),'tau2')

# $$T_1 = 10^{\log(T_1)}$$
# $$R_1 = 10^{-\log(T_1)}$$
# $$\ln(R_1) = -\log(T_1) \ln(10)$$

basis = (1-2*exp(-tau1*10**-LT1))*exp(-tau2*10**-LT2)
basis.shape

# In[3]:

simulated_data = basis*exact_data
simulated_data.sum(LT1_name).sum(LT2_name)

# In[4]:

simulated_data.nnls(('tau1','tau2'),(LT1,LT2),
                    (lambda tau1,LT1: 1-2*exp(-tau1*10**-LT1),
                     lambda tau2,LT2: exp(-tau2*10**-LT2)), l='BRD')

# In[5]:


