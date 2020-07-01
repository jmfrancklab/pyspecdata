# for 1D BRD, adapted mainly from Venkataramanan 2002
# but checked against BRD 1981
from pyspecdata import *
from scipy.optimize import nnls
# vd list
y_axis = linspace(5e-4,10,25)
# T1 axis
x_axis = r_[-4:2:100j]
def Gaussian_1d(axis,mu1,sigma1):
    this_G = exp(-(axis-mu1)**2/2/sigma1**2)
    return this_G
true_F = Gaussian_1d(10**(x_axis),6,0.3)

y_axis_2d = reshape(y_axis,(shape(y_axis)[0],1))
x_axis_2d = reshape(x_axis,(1,shape(x_axis)[0]))

print(shape(y_axis_2d))
print(shape(x_axis_2d))

K = (1.-2*exp(-y_axis_2d/10**(x_axis_2d)))
print(shape(K))
print(shape(true_F))

M = K.dot(true_F)
print(shape(M))
M = nddata(M,['vd'])
M.setaxis('vd',y_axis)
M.add_noise(0.2)
M_nd = M.C

# this is here to test the integrated 1D-BRD (for pyspecdata)
solution = M.C.nnls('vd',nddata(x_axis,'logT1'), lambda x,y: 1-2*exp(-x/10**(y)), l='BRD')

M = M.data

x,rnorm = nnls(K,M)
print(shape(x))

def nnls_reg(val):
    x_norm = empty([1])
    r_norm = empty([1])
    x,r_norm[0] = nnls(A_prime(val,dimension),b_prime)
    return x

# generate the A matrix, which should have form of the original kernel
# and then an additional length corresponding to size of the data dimension, where smothing param val is placed 
def A_prime(val,dimension):
    A_prime = r_[K,val*eye(dimension)]
    return A_prime

plot_Lcurve = False
#{{{ L-curve
l = sqrt(logspace(-10,1,25)) # adjusting the left number will adjust the right side of L-curve

M_nd.setaxis('vd',y_axis)
x_nd = nddata(x_axis,'logT1')
def vec_lcurve(l):
    return M_nd.real.C.nnls('vd',
            x_nd,lambda x,y: (1.-2*exp(-x/10**(y))), l=l)

# solution matrix for l different lambda values
x = vec_lcurve(l)
print(ndshape(x))
# norm of the residual (data - soln)
r_norm = x.get_prop('nnls_residual').data
# norm of the solution (taken along the fit axis)
x_norm = x.C.run(linalg.norm,'logT1').data

if plot_Lcurve:
    # Next plot the L-curve
    figure();title('L-Curve')
    # I do not actually know why we take the log, but this is important for the shape
    plot(log10(r_norm[:]),log10(x_norm[:]),'.')
    annotate_plot = True
    show_lambda = True
    if annotate_plot:
        if show_lambda:
            for j,this_l in enumerate(l):
                annotate('%0.4f'%this_l, (log10(r_norm[j]),log10(x_norm[j])),
                         ha='left',va='bottom',rotation=45)
        else:
            for j,this_l in enumerate(l):
                annotate('%d'%j, (log10(r_norm[j]),log10(x_norm[j])),
                         ha='left',va='bottom',rotation=45)
    show();quit()
#}}}

# From L-curve
this_L = 0.226

# generate data vector for smoothing
M = M[:,newaxis]
M_prime = r_[M,zeros((K.shape[1],1))]
M_prime = M_prime.squeeze()

print(shape(M)[0])

dimension = K.shape[1]
b_prime = M_prime

L_opt_vec = nnls_reg(this_L)

figure();title('ILT distributions')
true_F = nddata(true_F,'log(T1)')
true_F.setaxis('log(T1)',x_axis)
L_opt_vec = nddata(L_opt_vec,'log(T1)')
L_opt_vec.setaxis('log(T1)',x_axis)
plot(true_F,label='True')
plot(L_opt_vec,label='L-Curve')
plot(solution,':',label='pyspecdata-BRD')
legend()
show();quit()

