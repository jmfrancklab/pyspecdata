# # This is final version of the simulation and processing code

# based on DOI 10.1109/78.995059 (Venkataramanan 2002)

get_ipython().magic(u'load_ext pyspecdata.ipy')
from scipy.optimize import nnls
import IPython.display as d

# In[1]:

get_ipython().magic(u'matplotlib inline')
from matplotlib.pylab import *
# better for plotting distributions, 2D data

# ## PART 0. Generate test dataset

# Functions for generating Gaussian and correleated Gaussian distributions

def Gaussian_2D(log_x,log_y,x_mu,y_mu,sigma_x,sigma_y):
    log_x_mu = log10(10**(x_mu))
    log_y_mu = log10(10**(y_mu))
    this_Gaussian = exp(-(log_x-log_x_mu)**2/2/sigma_x**2
        -(log_y-log_y_mu)**2/2/sigma_y**2)
    return this_Gaussian

def corr_Gaussian_2D(theta,log_x,log_y,x_mu,y_mu,sigma_x,sigma_y):
    log_x_mu = log10(10**(x_mu))
    log_y_mu = log10(10**(y_mu))
    sigma1 = sigma_x
    sigma2 = sigma_y
    log_x_diff = cos(theta)*log_x +sin(theta)*log_y
    log_x_sum = -sin(theta)*log_x + cos(theta)*log_y
    log_x_diff_mu = cos(theta)*log_x_mu +sin(theta)*log_y_mu
    log_x_sum_mu = -sin(theta)*log_x_mu + cos(theta)*log_y_mu
    this_Gaussian = exp(-(log_x_diff-log_x_diff_mu)**2/2/sigma1**2
            -(log_x_sum-log_x_sum_mu)**2/2/sigma2**2)  
    return this_Gaussian

# Generate test T1,T2 distribution
# (log axes are for consistency with Venkatarmanan)

ridge = True
multi_feature = True

Nx = 50 # Number of T1 values
Ny = 50 # Number of T2 values
log_Nx_ax = linspace(log10(3e-3),log10(3),Nx) # NOTE: log_Nx_ax == 10**log_Nx_ax
log_Ny_ax = linspace(log10(3e-3),log10(2),Ny)
log_x = nddata(log_Nx_ax.copy(), r'log$(T_{1})$')
log_y = nddata(log_Ny_ax.copy(), r'log$(T_{2})$') # log scale -- coincides with log(T1),log(T2) axes of plot

x_mu = -1.0
y_mu = -0.4
sigma_x = 0.05
sigma_y = 0.5

if ridge:
    theta = -45*pi/180
    log_rho = corr_Gaussian_2D(theta,log_x,log_y,x_mu,y_mu,sigma_x,sigma_y)/2
if not ridge:
    log_rho = Gaussian_2D(log_x,log_y,x_mu,y_mu,sigma_x,sigma_y)/2

    
if multi_feature:
    x_mu2 = -1.25
    y_mu2 = -1.75
    sigma_x = 0.1
    sigma_y = 0.1
    log_rho += Gaussian_2D(log_x,log_y,x_mu2,y_mu2,sigma_x,sigma_y)

figure('density')
title(r'F(log$(T_{1})$,log$(T_{2})$)')
log_rho.reorder(r'log$(T_{2})$')
log_rho.contour(labels=False)
log_rho.reorder(r'log$(T_{1})$')
log_rho_x = log_rho.data

# Calculate signal for range of N1,N2 dimensions, plot

figure('signal')
N1_ax = logspace(log10(5.0e-4),log10(4),30)
N2_ax = linspace(5.0e-4,3.8,1000)
N1 = nddata(N1_ax.copy(), 'N1') # Number of tau1 values (first indirect dimension)
N2 = nddata(N2_ax.copy(), 'N2') # Number of tau2 values (second indirect dimension)
s = log_rho*exp(-N2/(10**log_y))*(1.-2*exp(-N1/(10**log_x)))
s.sum(r'log$(T_{1})$').sum(r'log$(T_{2})$')
s = s/amax(s.data) # scale to 1
noise_level = 0.035 # gives SNR of approx 30 dB
s.add_noise(noise_level)
s.reorder('N1')
s.name('Simulated signal intensity')
s.rename('N1',r'$\tau_{1}$')
s.rename('N2',r'$\tau_{2}$')
s.meshplot(cmap=cm.viridis)
snr_db = 20*log10(amax(s.data)/noise_level)
print("SNR",snr_db)

# Check that data has proper dimensions

print(ndshape(s))
data = s.data
print(shape(data))

# Switching to 4-dimensional arrays to simplify matrix math with numpy

N1_4d = reshape(N1_ax,(shape(N1_ax)[0],1,1,1))
N2_4d = reshape(N2_ax,(1,shape(N2_ax)[0],1,1))
Nx_4d = reshape(10**log_Nx_ax,(1,1,shape(log_Nx_ax)[0],1))
Ny_4d = reshape(10**log_Ny_ax,(1,1,1,shape(log_Ny_ax)[0]))
print(shape(log_Nx_ax),shape(Nx_4d))
print(shape(log_Ny_ax),shape(Ny_4d))
print(shape(N1_ax),shape(N1_4d))
print(shape(N2_ax),shape(N2_4d))

# Constructing kernels then SVD

k1 = (1.-2*exp(-N1_4d/Nx_4d))
k2 = exp(-N2_4d/Ny_4d)
print(shape(k1))
print(shape(k2))
k1_sqz = squeeze(k1)
k2_sqz = squeeze(k2)
U1, S1_row, V1 = np.linalg.svd(k1_sqz, full_matrices = False)
print(list(map(lambda x: x.shape, (U1, S1_row, V1))))
U2, S2_row, V2 = np.linalg.svd(k2_sqz, full_matrices = False)
print(list(map(lambda x: x.shape, (U2, S2_row, V2))))

# Plot singular values

figure('singular values')
semilogy(S1_row,'o-',label='S1',alpha=0.2)
semilogy(S2_row,'o-',label='S2',alpha=0.2)
xlabel('Index')
ylabel('Singular values')
xlim(0,20)
axhline(10**-2, c = 'k')
legend()

# Find the number of significant singular values in each kernel

S1_max = amax(S1_row)
S2_max = amax(S2_row)
tensor_max = S1_max*S2_max
print(tensor_max)
S1_condition_n = amax(S1_row)/amin(S1_row)
print(S1_condition_n)
S2_condition_n = amax(S2_row)/amin(S1_row)
print(S2_condition_n)

condition_list = []
for S1_i,S2_i in zip(range(shape(S1_row)[0]),range(shape(S2_row)[0])):
    condition_ratio = tensor_max/(S1_row[S1_i]*S2_row[S2_i])
    condition_list.append(condition_ratio)
for list_index,value in enumerate(condition_list):
    if value > 1000:
        truncation_index = list_index+1
        print("Can truncate singular values to index position",truncation_index)
        print("This assures a condition number of at least 1000")
        break

# Define the truncation here

choose_s1 = truncation_index+2
choose_s2 = truncation_index+2

# Perform truncation

print("Uncompressed singular row vector for K2",S1_row.shape)
S1_row = S1_row[0:choose_s1]
print("Compressed singular value row vector for K2",S1_row.shape)
V1 = V1[0:choose_s1,:]
U1 = U1[:,0:choose_s1]
print("Compressed V matrix for K2",V1.shape)
print("Compressed U matrix for K2",U1.shape)


print("Uncompressed singular row vector for K2",S2_row.shape)
S2_row = S2_row[0:choose_s2]
print("Compressed singular value row vector for K2",S2_row.shape)
V2 = V2[0:choose_s2,:]
U2 = U2[:,0:choose_s2]
print("Compressed V matrix for K2",V2.shape)
print("Compressed U matrix for K2",U2.shape)


I_S1 = eye(S1_row.shape[0]) 
S1 = S1_row*I_S1
print("Singular value matrix for K1",S1.shape)
I_S2 = eye(S2_row.shape[0])
S2 = S2_row*I_S2
print("Compressed singular value matrix for K2",S2.shape)


print(shape(U1))
print(shape(data))
print(shape(U2))

# Generate projected data, plot

data_proj = U1.dot(U1.T.dot(data.dot(U2.dot(U2.T))))
for tau_1_index in range(shape(data_proj)[0]):
    plot(data_proj[tau_1_index,:])
    title('Projected data')

# This is a meshplot of projected data for easier visualization

nd_proj = nddata(data_proj,['N1','N2'])
nd_proj.name('Projected data')
nd_proj.setaxis('N1',s.getaxis(r'$\tau_{1}$')).rename('N1',r'$\tau_{1}$')
nd_proj.setaxis('N2',s.getaxis(r'$\tau_{2}$')).rename('N2',r'$\tau_{2}$')
nd_proj.meshplot(cmap=cm.viridis)

# Generate compressed data, plot

print("Projected data dimensions:",shape(data_proj))
data_compr = U1.T.dot(data.dot(U2))
print("Compressed data dimensions:",shape(data_compr))
for tau_1_compressed_index in range(shape(data_compr)[0]):
    plot(data_compr[tau_1_compressed_index,:],'.-')
    title('Compressed data')

# This generates a compressed data plot as shown in Venkatarmanan

comp = data_compr
comp = reshape(comp,(shape(data_compr))[0]*(shape(data_compr))[1])
print(comp)
plot(comp,'-.')
ylabel('Compressed data')
xlabel('Index')

# Calculating the tensor kernel (K0)

K1 = S1.dot(V1)
K2 = S2.dot(V2)
print("Compressed K1",shape(K1))
print("Compressed K2",shape(K2))

K1 = reshape(K1, (shape(K1)[0],1,shape(K1)[1],1))
K2 = reshape(K2, (1,shape(K2)[0],1,shape(K2)[1]))
K0 = K1*K2
K0 = reshape(K0, (shape(K1)[0]*shape(K2)[1],shape(K1)[2]*shape(K2)[3]))
print("Compressed tensor kernel",shape(K0))
print("* Should be (",shape(S1)[0],"*",shape(S2)[0],") x (",shape(Nx_4d)[2],"*",shape(Ny_4d)[3],")")

# END PART 0.

# ## PART I. NNLS without regularization

# Generate lexicographically ordered data

datac_lex = [] # like b_vector
for m in range(shape(data_compr)[0]):
    for l in range(shape(data_compr)[1]):
        temp = data_compr[m][l]
        datac_lex.append(temp)
print("Dimension of lexicographically ordered data:",shape(datac_lex)[0])
print("Should match first dimension of compressed tensor kernel K0 which is",shape(K0)[0])

# Find solution vector then reorder into matrix

x, rnorm = nnls(K0,datac_lex)
solution = reshape(x,(Nx,Ny))

# Plot the result using numpy, compare to original

figure()
title('True F(T1,T2)')
image(log_rho_x)
figure()
title('Estimate F(T1,T2) - no smoothing')
image(solution)

# END PART I.

# ## PART II. Implement reguliarzation

# ### PART II.A Prepare data for regularization algorithm

datac_lex = array(datac_lex) # Add second dimension of 1 to lexicographically ordered data
datac_lex = datac_lex[:,newaxis]
print(shape(datac_lex))

dimension = K0.shape[1]
def A_prime(val,dimension):
    A_prime = r_[K0,val*eye(dimension)]
    return A_prime # Expands A matrix for regularization

b_prime = r_[datac_lex,zeros((dimension,1))] # Define b vector once
print(shape(b_prime))
b_prime = b_prime.squeeze() # for NNLS input

m_vec = datac_lex # For consistency with Venkatarmanan, define data vector as m

# ### PART II.B Use modified BRD to find smoothing parameter

# Functions below are defined for the BRD algorithm in Venkataramanan

def nnls_reg(val): # Performs regularized NNLS for given smoothing parameter
    x_norm = empty([1])
    r_norm = empty([1])
    x,r_norm[0] = nnls(A_prime(val,dimension),b_prime)
    return x

def heaviside(product):
    if product < 0:
        return 0
    if product >= 0:
        return 1

def square_heaviside(x_vec):
    diag_heavi = []
    for q in range(shape(K0.T)[0]):
        pull_val = dot(K0.T[q,:],x_vec)
        diag_heavi.append(heaviside(pull_val[0]))
    diag_heavi = array(diag_heavi)
    square_heavi = diag_heavi*eye(shape(diag_heavi)[0])
    return square_heavi

def G(x_vec):
    #print shape(x_vec)
    return dot(K0,dot(square_heaviside(x_vec),K0.T))

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
    tol = 1e-2
    vector_converged=False
    err_list = []
    for iter in range(maxiter):
        fder = dd_chi(G(input_vec),val)
        fval = d_chi(input_vec,val)
        newton_step = dot(linalg.inv(fder),fval)
        update_vec = input_vec + newton_step
        if alpha_converged:
            match = 0
            err_list = []
            for update_index,update_val in enumerate(update_vec):
                for input_index,input_val in enumerate(input_vec):
                    error = abs(update_val[0] - input_val[0])
                    err_list.append(error)
                    if error <= tol:
                        match = match + 1
            if match == len(update_vec[0]):
                print("Fully converged")
                print(err_list)
                vector_converged = True
    return update_vec,vector_converged,err_list

def optimize_alpha(input_vec,val):
    alpha_converged = False
    fac = sqrt(choose_s1*choose_s2)
    T = linalg.inv(dd_chi(G(input_vec),val**2))
    dot_prod = dot(input_vec.T,dot(T,input_vec))
    ans = dot_prod*fac
    ans = ans/linalg.norm(input_vec)
    ans = ans/(dot_prod)
    tol = 1e-3
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
        c_vec = dot(K0,f_vec) - m_vec
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
    print(err_list)
    return new_lambda

# ### PART II.B1 Initial guess for smoothing parameter
# 
# 
# 
# NOTE: $\alpha = \lambda^2$

S_curve = True # This will generate the data for an S-Curve
if S_curve:
    rnorm_list = []
    smoothing_list = []
    alt_norm_list = []
    lambda_range = logspace(log10(8e-4),log10(2e4),10)
    for index,lambda_val in enumerate(lambda_range):
        print("index",index)
        soln,temp_rn = nnls(A_prime(lambda_val,dimension),b_prime)
        rnorm_list.append(temp_rn)
        f_vec = soln[:,newaxis]
        alpha = lambda_val**2
        c_vec = dot(K0,f_vec) - m_vec
        c_vec /= -alpha
        alt_temp = linalg.norm(c_vec)*alpha
        alt_norm_list.append(alt_temp)
        smoothing_list.append(lambda_val)

# Here you can generate the S-Curve using the formula given in Venkataramanan for the norm
# or you can use the norm returned from the NNLS algorithm -- plots should be the same shape

if S_curve:  
    figure('using norm')
    rnorm_axis = array(rnorm_list)
    smoothing_axis = array(smoothing_list)
    plot(log10(smoothing_axis**2),rnorm_axis)
    figure();title('using LV norm')
    altnorm_axis = array(alt_norm_list)
    smoothing_axis = array(smoothing_list)
    plot(log10(smoothing_axis**2),altnorm_axis,'-.',c='k')
    xlabel(r'log($\alpha$)')
    ylabel(r'$\chi$($\alpha$)')
    gridandtick(gca())
    # Annotations for chosen heel below
    #plot(1.96,0.69,marker='x',markersize=12,color='red') # initial value
    #plot(1.396,0.50,marker='x',markersize=12,color='orange') # optimized value

# This takes the chosen heel value (of log(alpha)) and converts to lambda for input into nnls

if S_curve:
    heel = 1.96
    heel_alpha = 10**heel
    heel_lambda = sqrt(heel_alpha)
    print("Alpha",heel_alpha)
    print("Lambda",heel_lambda)

# Also able to generate L Curve if preferred

L_curve = False # This will generate an L-Curve, adapted from JF and EGR
if L_curve:
    l = sqrt(logspace(0.5,10,10))
    x_norm = empty([shape(l)[0]])
    r_norm = empty([shape(l)[0]])
    print(shape(x_norm))
    for j, lambda_val in enumerate(l):
        x,r_norm[j] = nnls(A_prime(lambda_val,dimension),b_prime)
        x_norm[j] = linalg.norm(x)
    plot(log10(r_norm),log10(x_norm),'.')
    for j, val in enumerate(l):
        annotate('%5g'%val, (log10(r_norm[j]),log10(x_norm[j])),
        ha='left',va='bottom',rotation=45)
    ylabel('$x$ norm')
    xlabel('residual')

# Use the initial guess for the smoothing parameter to generate a distribution

figure();title(r'True F(log$(T_{1})$,log$(T_{2})$')
image(log_rho_x)
opt_vec = nnls_reg(heel_lambda)
solution = reshape(opt_vec,(Nx,Ny))
figure();title(r'Est F(log$(T_{1}$),log$(T_{2})$, $\lambda$ = %0.2f (SNR %d dB)'%(heel_lambda,snr_db))
image(solution)

# ### PART II.B2 Use the initial guess as input to BRD algorithm

# Run the BRD algorithm to find optimum alpha

opt_val = mod_BRD(guess=heel_lambda,maxiter=20)
print("OPTIMIZED LAMBDA:",opt_val)

# Reshape and visualize results

opt_vec = nnls_reg(opt_val) # Find solution vector using BRD output
solution = reshape(opt_vec,(Nx,Ny))
figure();title(r'True F(log$(T_{1})$,log$(T_{2})$')
image(log_rho_x)
figure();title(r'Est F(log$(T_{1}$),log$(T_{2})$, $\lambda$ = %0.2f (SNR %d dB)'%(opt_val,snr_db))
image(solution)

# ### PART II.B3

# Convert solution to pyspecdata object for contour plot

nd_solution = nddata(solution,[r'log$(T_{1})$',r'log$(T_{2})$'])

# Labelling axes, titles

nd_solution.setaxis(r'log$(T_{1})$',log_Nx_ax.copy())
nd_solution.setaxis(r'log$(T_{2})$',log_Ny_ax.copy())
print(ndshape(nd_solution))
print(ndshape(log_rho))

figure();title(r'True $F(log(T_{1}),log(T_{2})$')
log_rho.reorder(r'log$(T_{2})$') # For consistency with Venkataramanan
log_rho.contour(labels=False)

figure();title(r'Estiamted F(log$(T_{1})$,log($T_{2}$), $\lambda$ = %0.2f (SNR %d dB)'%(opt_val,snr_db))
nd_solution.reorder(r'log$(T_{2})$')
nd_solution.contour(labels=False)

# ### END PART II.

# ### PART III.

# Can repeat with greater size of direct dimensions (i.e., finer sampling) to produce more satisfactory results

