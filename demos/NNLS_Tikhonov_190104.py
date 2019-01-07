
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
fl=figlist_var()


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

m = A.shape[1]
l = sqrt(logspace(-8,4,10)) # I do this because it gives me a fairly even spacing of points
x_norm = empty_like(l)
r_norm = empty_like(l)
for j,this_l in enumerate(l):
    x,r_norm[j] = nnls_regularized(A,test_signal.squeeze(),l=this_l)
    x_norm[j] = linalg.norm(x)
fl.next('L-curve')
plot(log10(r_norm),log10(x_norm),'.')
for j,this_l in enumerate(l):
    annotate('%5g'%this_l, (log10(r_norm[j]),log10(x_norm[j])),
             ha='left',va='bottom',rotation=45)
ylabel('$x$ norm')
xlabel('residual')


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
