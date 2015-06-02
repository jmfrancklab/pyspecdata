from .propagator import *
from time import *
from .matlablike import *
class CustomError(Exception):
	def __init__(self, *value):
		if len(value)>1:
			self.value = map(str,value)
		else:
			self.value = value
	def __str__(self):
		retval = ' '.join(self.value)
		return '\n'+'\n'.join(textwrap.wrap(retval,90))
def check_for_finite(*args):
	r''' checks that args[1:] are all finite
uses the first thing in list args[0] to identify the function, and the rest to identify the arguments'''
	for j in range(1,len(args)):
		if any(~isfinite(r_[tuple(map((lambda x: x.flatten()),args[j]))])): # the weirdness compresses the list of arrays into a single array
			raise CustomError('in function',args[0][0],',',args[0][j],"isn't finite!")
		else:
			print 'verified that',args[0][j],'was finite'
def prop(t=1,p=(1),H=1,accumtime = False):
	r'''	 now takes t, p, and H in C order
	 the 'F' flag doesn't seem to operate in a super sensible way.
	 therefore, we just use .T to convert to and from the fortran arrays -- this preserves memory order, and keeps all the data in the expected place
	 otherwise, the INDECES are preserved, not the memory order, which makes things very confusing when using the array([[1,2],[3,4]]) format, or python reshaping!, information returned in U[p1...pN,t,H1,H2] format'''
	check_for_finite(['prop','t','p','H'],t,p,H)
	n_t = shape(t)[1]
	n_h = shape(H)[2]
	n_operators = shape(H)[0]
	
	# for each hamiltonian, define the multiplier p at each point on the n-d grid
	#{{{generate p_collapsed with the parameters collapsed into one dim
	p_shape = array(p[0].shape)
	for j in range(1,len(p)):
		print 'shape of parameter',j,shape(p[j])
		p_shape = array([array(p[j].shape),p_shape]).max(0) # take the max along the outer dim
	print 'p_shape: ',p_shape
	n_p = p_shape.prod()
	p_ones = ones(p_shape,dtype='double')
	p_collapsed = zeros((n_operators,n_p),dtype='double')
	print 'p_ones: ',shape(p_ones)
	print 'p_collapsed: ',shape(p_collapsed)
	for j in range(0,len(p)):
		p_collapsed[j,:] = (p_ones * p[j].copy() #expands to the correct number of dim's
				).flatten()
	print p_collapsed
	#}}}
	#{{{ do the fortran propagation
	if accumtime:
		propagator.propagators = zeros((n_p,n_t,n_h,n_h),dtype='complex128').T
	else:
		propagator.propagators = zeros((n_p,1,n_h,n_h),dtype='complex128').T
	print 'Hamiltonians:',H.T
	propagator.operators = H.T.copy() # here, i label as fortran, and transpose, to preserve the same innermost --> outermost order
	print 'time dependence:',t
	propagator.time_dependence = t.copy() # here, i leave out the transpose to make the operator dimension go on the inside
	propagator.p = p_collapsed # same here
	print 'getting ready to actually propagate'
	print 'fortran time:',
	storetime = time()
	propagator.propagate()
	propagator.nothing()
	print 'nothing next to propagate returned without segfault'
	print time() - storetime
	#}}}
	if accumtime:
		return reshape(propagator.propagators.T.copy(),r_[p_shape,n_t,n_h,n_h]) # again, use the T to convert! (note that after the T, it is also labelled as C order again!
	else:
		return reshape(propagator.propagators.T.copy(),r_[p_shape,1,n_h,n_h]) # again, use the T to convert! (note that after the T, it is also labelled as C order again!
def sandwich(U,rho):
	r'''	this does the similarity transform
	of U[p1...pN,h1,h2] on rho[] either as a matrix, or in the same form
	--> need to verify this'''
	check_for_finite(['sandwich','U','rho',],U,rho)
	# here, it's going to be best to just use a compiled subroutine
	if len(shape(rho))<4:
		rho = reshape(rho,r_[1,1,array(shape(rho))])
	# I need to collapse the parameter dimensions of U, down to 4, where the outer 2 are nominally time and parameter
	U_oldshape = array(shape(U))
	if len(U_oldshape)<4:
		U_newshape = ones(5)
		U_newshape[-len(U_oldshape):] = U_oldshape
	elif len(U_oldshape)>4:
		U_newshape = r_[U_oldshape[:-3].prod(),U_oldshape[-3:]]
	else:
		U_newshape = U_oldshape
	U = reshape(U,U_newshape)
	propagator.propagators = U.T.copy()
	propagator.rho = rho.T.copy()
	propagator.sandwich()
	U = reshape(propagator.propagators.T.copy(),U_oldshape)
	return U
def lvdot(rho,C):
	r'''	this does the liouville space dot product 
	of U[p1...pN,h1,h2] with the matrix C either as a matrix,
	or apparently as C[p,h1,h2], though it doesn't seem to
	accept multidimensional parameter arrays --> I should fix this'''
	check_for_finite(['lvdot','rho','C',],rho,C)
	# I need to collapse dimensions of rho down to 4, where the outer 2 are nominally time and parameter
	rho_oldshape = array(shape(rho))
	if len(rho_oldshape)==3:
		raise CustomError("You're trying to pass an array with 3 dimensions, and I don't know what this means! --> should be [p,t,h1,h2]\nuse -1: or 2:3 to not drop the dim")
	if len(rho_oldshape)<4:
		rho_newshape = ones(5)
		rho_newshape[-len(rho_oldshape):] = rho_oldshape
	elif len(rho_oldshape)>4:
		rho_newshape = r_[rho_oldshape[:-3].prod(),rho_oldshape[-3:]]
	else:
		rho_newshape = rho_oldshape
	rho = reshape(rho,rho_newshape)
	propagator.traceout = zeros(rho_newshape[:-2],dtype='complex128').T
	propagator.rho = rho.T.copy()
	C_shape = array(C.shape)
	if len(C_shape) == 2:
		C = reshape(C,r_[1,C_shape])
	propagator.c = C.T.copy()
	propagator.lvdot()
	lvdot = reshape(propagator.traceout.T.copy(),rho_oldshape[:-2])
	return lvdot
def grape(U,C,H_c,rho_0):
	r''' and this guy returns to gradient --> add details'''
	U_oldshape = array(shape(U))
	if len(U_oldshape)<4:
		U_newshape = ones(5)
		U_newshape[-len(U_oldshape):] = U_oldshape
	elif len(U_oldshape)>4:
		U_newshape = r_[U_oldshape[:-3].prod(),U_oldshape[-3:]]
	else:
		U_newshape = U_oldshape
	U = reshape(U,U_newshape)
	propagator.propagators = U.T
	C_shape = array(C.shape)
	if len(C_shape) != 3:
		C = reshape(C,(r_[C_shape[:-2].prod(),C_shape[-2:]]))
	rho_shape = array(shape(rho_0))
	if len(rho_shape)==2:
		rho_0 = reshape(rho_0,r_[1,1,rho_shape])
	else:
		rho_0 = reshape(rho_0,r_[rho_shape[:-2].prod,1,rho_shape[-2:]])
	print 'shape of rho_0',rho_shape
	print propagator.rho.dtype
	print propagator.rho.dtype
	propagator.rho = rho_0.T
	propagator.c = C.T.copy()
	propagator.h_c = H_c.T
	propagator.gradient = zeros(propagator.propagators.shape[2:])
	propagator.grape()
	return reshape(propagator.gradient.T.copy(),U.shape[:-2])
