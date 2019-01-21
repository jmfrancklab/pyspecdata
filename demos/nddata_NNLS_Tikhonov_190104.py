#!/usr/bin/env python
# coding: utf-8

# 


try:
    get_ipython().magic(u'load_ext pyspecdata.ipy')
    in_notebook = True
except:
    from pyspecdata import *
    in_notebook = False
from pyspecdata import nnls_regularized
from numpy import random
import time
fl=figlist_var()
init_logging(level='debug')


# got the following from here:
# https://medium.com/pythonhive/python-decorator-to-measure-the-execution-time-of-methods-fa04cb6bb36d

l_line = ''
def timeit(method,n_times=1):
    def timed(*args, **kw):
        timing = zeros(n_times+1)
        timing[0] = time.time()
        for j in range(n_times):
            result = method(*args, **kw)
            timing[j+1] = time.time()
        time_diff = (timing[-1]-timing[0])/float(n_times)
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int(time_diff * 1000)
        else:
            print '%r  %2.2f ms (average of %d runs)' %                   (method.__name__, time_diff * 1000, n_times) + l_line
        return result
    return timed


# 


R = c_[1.:100:500j] # distribution of T2 relaxation rates
peaks = [(80,4,1),(20,0.5,0.5),(30,0.5,0.25)]
calcd = False
for mu,sigma,A in peaks:
    if not calcd:
        P = A*exp(-(R-mu)**2/(2*sigma**2))
        calcd = True
    else:
        P += A*exp(-(R-mu)**2/(2*sigma**2))


# Vary R as we move along the rows

# 


P = P.T
R = R.T
fl.next('distribution function')
plot(R.flatten(),P.flatten())


# 


endp = 0.2
t = c_[1e-3:endp:2048j] # column vectors give functions of time
test_signal = exp(-R*t).dot(P.T)

test_signal += random.normal(scale = 0.01,size=(2048,1))
print test_signal.shape
print t.squeeze().shape
fl.next('test data function')
test_data = nddata(test_signal.flatten(),[-1],['t']).labels('t',t)
plot(test_data)
xlim(-endp/10,endp)

# Do the basic NNLS fit

test_fit = test_data.C.nnls('t',{'R':R.ravel()},lambda x,y: exp(-y*x))
fl.next('fit an exponential distribution',legend=True)
fl.plot(test_data, alpha=0.5, label='test signal')
K = test_fit.get_prop('nnls_kernel').C
#note that order doesn't matter for the following dot (done by dimension name)
fl.plot(test_fit.C.dot(K), alpha=0.5, label='fit')
fl.next('what does the fit distribution look like?')
fl.plot(test_fit)

# Now add regularization

def L_curve(l,r_norm,x_norm, **kwargs):
    """plot L-curve using

    Parameters
    ==========
    l: double
        lambda values
    r_norm: double
        norm of the residual
    x_norm: double
        norm of solution vector"""
    plot(log10(r_norm),log10(x_norm),'o',**kwargs)
    for j,this_l in enumerate(l):
        annotate('%5g'%this_l, (log10(r_norm[j]),log10(x_norm[j])),
                 ha='left',va='bottom',rotation=45)
    ylabel('$\log_{10}(x$ norm$)$')
    xlabel('$\log_{10}($ residual $)$')

# 

l = sqrt(logspace(-8,4,10)) # I do this because it gives me a fairly even spacing of points
@timeit
def nonvec_lcurve(A,l):
    x_norm = empty_like(l)
    r_norm = empty_like(l)
    for j,this_l in enumerate(l):
        x = test_signal.C.nnls('t',
                {'R':R.ravel()},lambda x,y: exp(-y*x))
        r_norm[j] = x.get_prop('nnls_residual')
        x_norm[j] = linalg.norm(x.data)
    return x,x_norm,r_norm
x,x_norm,r_norm = nonvec_lcurve(A,l)
#x_norm = map(linalg.norm,x) # to be fair, this calculation is done outside the timing, below


# 


fl.next('L-curve', legend=True)
L_curve(l, r_norm, x_norm, markersize=10, alpha=0.5, label='manual loop')
