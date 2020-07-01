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

val = 1
dimension = K.shape[1]
dimension2 = M.shape[0]
b_prime = M_prime
r_norm = empty([1])
test_x,r_norm[0] = nnls(A_prime(val,dimension),b_prime)
test_x_vec = test_x[:,newaxis]
m_vec = M
c_vec = dot(K,test_x_vec) - m_vec
c_vec /= -1

#{{{ BRD code
def heaviside(product):
    if product < 0:
        return 0
    if product >= 0:
        return 1

def square_heaviside(x_vec):
    diag_heavi = []
    for q in range(shape(K.T)[0]):
        pull_val = dot(K.T[q,:],x_vec)
        diag_heavi.append(heaviside(pull_val[0]))
    diag_heavi = array(diag_heavi)
    square_heavi = diag_heavi*eye(shape(diag_heavi)[0])
    return square_heavi

def G(x_vec):
    #print shape(x_vec)
    return dot(K,dot(square_heaviside(x_vec),K.T))

def chi(x_vec,val):
    return 0.5*dot(x_vec.T,dot(dd_chi(G(x_vec),val**2),x_vec)) - dot(x_vec.T,m_vec)

def d_chi(x_vec,val):
    return dot(dd_chi(G(x_vec),val**2),x_vec) - m_vec

def dd_chi(G,val):
    return G + (val**2)*eye(shape(G)[0])

def newton_min(input_vec,val,alpha_converged):
    if alpha_converged:
        maxiter = 100
    if not alpha_converged:
        maxiter = 1
    tol = 1e-6
    vector_converged = False
    err_list = []
    for iter in range(maxiter):
        fder = dd_chi(G(input_vec),val)
        fval = d_chi(input_vec,val)
        newton_step = dot(linalg.inv(fder),fval)
        update_vec = input_vec + newton_step
        if alpha_converged:
            match = 0
            err_list = []
            for update_index,update_val in enumerate(updated_vec):
                for input_index,input_val in enumerate(input_vec):
                    error = abs(update_val[0] - input_val[0])
                    err_list.append(error)
                    if error <= otl:
                        match = match+1
            if match == len(update_vec[0]):
                print("FULLY CONVERGED!")
                print(err_list)
                vector_converged = True
    return update_vec,vector_converged,err_list

def optimize_alpha(input_vec,val):
    alpha_converged = False
    fac = sqrt(input_vec.shape[0])
    T = linalg.inv(dd_chi(G(input_vec),val**2))
    dot_prod = dot(input_vec.T,dot(T,input_vec))
    ans = dot_prod*fac
    ans = ans/linalg.norm(input_vec)
    ans = ans/(dot_prod)
    tol = 1e-6
    if abs(ans-val**2) <= tol:
        print("Alpha has converged")
        alpha_converged = True
        return ans,alpha_converged
    return ans,alpha_converged
def mod_BRD(guess,maxiter=20):
    smoothing_param = guess
    alpha_converged = False
    vector_converged = False
    for iter in range(maxiter):
        print("*** *** ITERATION NO.",iter,"*** ***")
        print("*** CURRENT LAMBDA",smoothing_param," *** ")
        x_norm = empty([1])
        r_norm = empty([1])
        soln,r_norm = nnls(A_prime(smoothing_param,dimension),b_prime)
        f_vec = soln[:,newaxis]
        alpha = smoothing_param**2
        c_vec = dot(K,f_vec) - m_vec
        c_vec /= -alpha
        new_c,vector_converged,err_list = newton_min(c_vec,smoothing_param,alpha_converged)
        new_alpha,alpha_converged = optimize_alpha(new_c,smoothing_param)
        new_alpha = new_alpha[0,0]
        new_lambda = sqrt(new_alpha)
        if alpha_converged:
            print("*** OPTIMIZED LAMBDA",new_lambda," *** ")
            break  
        if not alpha_converged:
            print("*** UPDATED LAMBDA",new_lambda," *** ")
            smoothing_param = new_lambda
    return new_lambda
#}}}

opt_val = mod_BRD(guess=2,maxiter=20)
opt_vec = nnls_reg(opt_val)
L_opt_vec = nnls_reg(this_L)

figure();title('ILT distributions')
true_F = nddata(true_F,'log(T1)')
true_F.setaxis('log(T1)',x_axis)
opt_vec = nddata(opt_vec,'log(T1)')
opt_vec.setaxis('log(T1)',x_axis)
L_opt_vec = nddata(L_opt_vec,'log(T1)')
L_opt_vec.setaxis('log(T1)',x_axis)
plot(true_F,label='True')
plot(L_opt_vec,label='L-Curve')
plot(opt_vec,label='BRD')
plot(solution,':',label='pyspecdata-BRD')
legend()
show();quit()

