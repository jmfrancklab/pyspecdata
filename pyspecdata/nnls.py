# again, this is copied liberally from scipy nnls -- see scipy licensing
from __future__ import division, print_function, absolute_import

from . import _nnls
from .general_functions import redim_F_to_C, redim_C_to_F
from numpy import asarray_chkfinite, zeros, double, isscalar
from numpy import array as np_array
import multiprocessing.dummy as mpd
from multiprocessing import cpu_count

__all__ = ['nnls_regularized']


def nnls(self, dimname, newaxis_dict, kernel_func, l=0):
    r"""Perform regularized non-negative least-squares "fit" on self.

    .. todo::
        someone can explain the math here
    
    Parameters
    ==========
    dimname: str
        Name of the "data" dimension that is to be replaced by a
        distribution (the "fit" dimension);
        *e.g.* if you are regularizing a set of functions
        :math:`\exp(-\tau*R_1)`, then this is :math:`\tau`
    newaxis_dict: dict
        a dictionary whose key is the name of the "fit" dimension
        (:math:`R_1` in the example above)
        and whose value is an array with the new axis labels
    kernel_func: function
        a function giving the kernel for the regularization.
        The first argument is the "data" variable
        and the second argument is the "fit" variable
        (in the example above, this would be something like
        ``lambda x,y: exp(-x*y)``)
    l : double (default 0)
        the regularization parameter
        :math:`lambda` -- if this is set to 0, the algorithm reverts to
        standard nnls

    Returns
    =======
    self:
        The regularized result.
        For future use, both the kernel (as an nddata, in a property called
        "nnls_kernel") and the residual (as an nddata, in a property called
        "nnls_residual") are stored as properties of the nddata.
        The regularized dimension is always last
        (innermost).
    """
    assert type(dimname) is str, "first argument is dimension name"
    assert type(newaxis_dict) is dict, "second argument is dictionary with new axis"
    assert len(newaxis_dict) == 1, "currently only set up for 1D"
    assert callable(kernel_func), "third argument is kernel function"
    # construct the kernel
    # the kernel transforms from (columns) the "fit" dimension to (rows)
    # the "data" dimension
    fitdim_name = newaxis_dict.keys().item()
    fit_axis = nddata(newaxis_dict[fitdim_name], fitdim_name)
    data_axis = self.fromaxis(dimname)
    K = fit_axis*data_axis
    self.reorder(dimname, first=False) # make the dimension we will be regularizing innermost
    data_fornnls = self.data
    if len(data_fornnls) > 2:
        data_fornnls = data_fornnls.reshape(prod(
            data_fornnls.shape[:-1]),data_fornnls.shape[-1])
    retval, residual = nnls_regularized(K, data_fornnls, l=l)
    newshape = []
    if not isscalar(l):
        newshape.append(len(l))
        self.dimlabels = ['lambda'] + self.dimlabels
        self.setaxis('lambda',l)
    newshape += list(self.data.shape)
    newshape.append(len(newaxis_dict[fitdim_name]))
    retval = retval.reshape(newshape)
    self.data = retval
    # change the dimension names and data
    self.rename(dimname, fitdim_name)
    self.setaxis(fitdim_name, newaxis_dict[fitdim_name])
    self.data = retval
    # make the residual nddata as well
    residual_nddata = ndshape(self).pop(dimname).alloc()
    residual_nddata.data[:] = residual[:]
    # store the kernel and the residual in the properties
    self.set_prop('nnls_kernel',K)
    self.set_prop('nnls_residual',residual_nddata)
    return self
def nnls_regularized(A, b, l=0, maxiter=None):
    """
    Solve math:`argmin_x || Ax - b ||_2^2 + \lambda^2 ||x||_2^2` for ``x>=0``.
    This is a wrapper for a FORTRAN non-negative least squares solver,
    with regularization (added by stacking $A$ on top an identity matrix
    times $\lambda$ and $b$ on top of a matching array of zero.

    Parameters
    ----------
    A : ndarray
        Matrix ``A`` as shown above.
    b : ndarray
        Right-hand side vector.
    l : double (default 0)
        :math:`lambda` -- if this is set to 0, the algorithm reverts to
        standard nnls (rather than stacking on top of two zero matrices
        for no reason)
    maxiter: int, optional
        Maximum number of iterations, optional.
        Default is ``3 * A.shape[1]``.

    Returns
    -------
    x : ndarray
        Solution vector.
    rnorm : float
        The residual, ``|| Ax-b ||_2``.

    Notes
    -----
    The FORTRAN code was published in the book below. The algorithm
    is an active set method. It solves the KKT (Karush-Kuhn-Tucker)
    conditions for the non-negative least squares problem.

    This was adapted from the source distributed with scipy --
    see scipy for relevant licensing.

    References
    ----------
    Lawson C., Hanson R.J., (1987) Solving Least Squares Problems, SIAM

    """

    A, b = map(asarray_chkfinite, (A, b))

    if len(A.shape) != 2:
        raise ValueError("expected matrix")
    if len(b.shape) > 2:
        raise ValueError("expected vector")

    m, n = A.shape

    if m != b.shape[-1]:
        raise ValueError("incompatible dimensions (the most quickly changing index should be the nnls dimension)")

    maxiter = -1 if maxiter is None else int(maxiter)

    if isscalar(l):
        if l == 0.0:
            w = zeros((n,), dtype=double)
            zz = zeros((m,), dtype=double)
            index = zeros((n,), dtype=int)
            x, rnorm, mode = _nnls.nnls(A, b, w, zz, index, maxiter)
        else:
            w = zeros((n,), dtype=double)
            zz = zeros((m+n,), dtype=double)
            index = zeros((n,), dtype=int)
            # choose the correct subroutine based on the dimension
            if len(b.shape) == 1:
                x, rnorm, mode = _nnls.nnls_regularized(A, b, w, zz, index, maxiter, l)
            if len(b.shape) == 2:
                x, rnorm, mode = _nnls.nnls_regularized_loop(A, redim_C_to_F(b), w, zz, index, maxiter, l)
                x = redim_F_to_C(x)
    else:
            nCPU = cpu_count() 
            #print("I found",nCPU,"CPU's")
            p = mpd.Pool(nCPU)
            if len(b.shape) == 1:
                def nnls_func(l):
                    w = zeros((n,), dtype=double)
                    zz = zeros((m+n,), dtype=double)
                    index = zeros((n,), dtype=int)
                    x, rnorm, mode = _nnls.nnls_regularized(A, b, w, zz, index, maxiter, l)
            if len(b.shape) == 2:
                def nnls_func(l):
                    w = zeros((n,), dtype=double)
                    zz = zeros((m+n,), dtype=double)
                    index = zeros((n,), dtype=int)
                    x, rnorm, mode = _nnls.nnls_regularized_loop(A, redim_C_to_F(b), w, zz, index, maxiter, l)
                    x = redim_F_to_C(x)
            retval = p.map(nnls_func,l)
            x,rnorm,mode = map(np_array,zip(*retval))
    if (isscalar(mode) and mode != 1):
        # need something for the multiple lambda
        raise RuntimeError("too many iterations")
    return x, rnorm
