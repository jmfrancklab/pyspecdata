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


# 


# got the following from here:
#https://medium.com/pythonhive/python-decorator-to-measure-the-execution-time-of-methods-fa04cb6bb36d
l_line = ''
def timeit(method,n_times=5):
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
fl.next('distribution function')
plot(t.flatten(),test_signal.flatten())
xlim(-endp/10,endp)


# 


A = exp(-R*t)


# Do the basic NNLS fit

# 


print test_signal.shape
print A.shape
x,rnorm = nnls_regularized(A,test_signal.squeeze(),l=0.)
fl.next('fit an exponential',legend=True)
fl.plot(t[:],test_signal.flatten(),label='test signal')
fl.plot(t[:],A.dot(x),label='fit')


# 


print r_[c_[1:3:3j],zeros((2,1))]


# Now add regularization

# 


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
        x,r_norm[j] = nnls_regularized(A,test_signal.squeeze(),l=this_l)
        x_norm[j] = linalg.norm(x)
    return x,x_norm,r_norm
x,x_norm,r_norm = nonvec_lcurve(A,l)
#x_norm = map(linalg.norm,x) # to be fair, this calculation is done outside the timing, below


# 


fl.next('L-curve', legend=True)
L_curve(l, r_norm, x_norm, markersize=10, alpha=0.5, label='manual loop')


# 


print x_norm
print r_norm


# ## Vectorized version of lambda curve

# 


l = sqrt(logspace(-8,4,10)) # I do this because it gives me a fairly even spacing of points
@timeit
def vec_lcurve(A,l):
    return nnls_regularized(A,test_signal.squeeze(),l=l)
x,r_norm = vec_lcurve(A,l)


# 


print x.shape,r_norm.shape,l.shape,linalg.norm(x,axis=1).shape


# 


print linalg.norm(x,axis=1)
print r_norm


# 


fl.next('L-curve')
L_curve(l,r_norm,linalg.norm(x,axis=1), markersize=5, alpha=0.5, label='compiled loop')


# 


print shape(x)
print shape(x.T)


# 


P_estimated,final_rnorm = nnls_regularized(A,test_signal.squeeze(),l=0.1)
fl.next(r'show result where $\lambda$ set to knee')
fl.plot(R.flatten(),P.flatten())
fl.plot(R.flatten(),P_estimated,alpha=0.5,linewidth=2)


# Now compare to standard Tikhonov (no non-negative)

# 


m = A.shape[1]
l = sqrt(logspace(-8,4,10))
x_norm = empty_like(l)
r_norm = empty_like(l)
U,s,Vp = svd(A,full_matrices=False)
#print S
#print map(lambda x: x.shape, (U,S,Vp))
def calculate_x(this_l):
    """calculate the 'thresholded' $1/\Sigma$, which depends only on A and on
    $\lambda$, then use the SVD unitary matrices to calculate x"""
    D = diag(s/(s**2+this_l**2))
    x = Vp.T.dot(D).dot(U.T).dot(test_signal)
    return x
for j,this_l in enumerate(l):
    x = calculate_x(this_l)
    x_norm[j] = linalg.norm(x)
    r_norm[j] = linalg.norm(A.dot(x)-test_signal)
fl.next('standard (allow neg) L-curve')
plot(log10(r_norm),log10(x_norm),'.')
for j,this_l in enumerate(l):
    annotate('%5g'%this_l, (log10(r_norm[j]),log10(x_norm[j])),
             ha='left',va='bottom',rotation=45)
ylabel('$x$ norm')
xlabel('residual')


# 


#figure(figsize(8,4))
#gcf().add_axes([0.1, 0.1, 0.6, 0.5])

fl.next('compare to allow-negative and non-negative')
plot(R.flatten(), P.flatten(),label='Test distribution')
plot(R.flatten(), P_estimated, alpha=0.5, linewidth=2,label='NNLS Tikhonov (L-curve regularized)')
plot(R.flatten(), calculate_x(0.21), alpha=0.5, linewidth=2,label='standard Tikhonov (L-curve regularized)')
legend(**dict(bbox_to_anchor=(1.05,1),loc = 2,borderaxespad=0.))
xlabel('Decay rate $R$')
ylabel('Est. population of species with decay rate $P(R)$')
savefig('yerdon_170420.png',dpi=300,bbox_inches='tight')

if in_notebook:
    print "in notebook, not calling show"
else:
    print "not in notebook, calling show"
    fl.show()

