import logging
import numpy as np
from numpy import newaxis
from ..general_functions import strm
from .. import nnls as this_nnls

logger = logging.getLogger("pyspecdata.matrix_math")

# {{{ local functions
def chi(x_vec, val, data_fornnls, K):
    # appears to not be used?
    return 0.5 * np.dot(
        x_vec.T, np.dot(dd_chi(G(x_vec, K), val**2), x_vec)
    ) - np.dot(x_vec.T, data_fornnls[:, newaxis])

def d_chi(x_vec, val, data_fornnls, K):
    return np.dot(dd_chi(G(x_vec, K), val**2), x_vec) - data_fornnls[:, newaxis]

def dd_chi(G, val):
    return G + (val**2) * np.eye(np.shape(G)[0])

def G(x_vec, K)
    return np.dot(K, np.dot(square_heaviside(x_vec, K), K.T))

def H(product):
    if product <= 0:
        return 0
    if product > 0:
        return 1

def square_heaviside(x_vec, K):
    K_dim = np.shape(K.T)[0]
    diag_heavi = np.empty(K_dim)
    for q in range(K_dim):
        pull_val = np.dot(K.T[q, :], x_vec)
        temp = pull_val[0]
        diag_heavi[q] = H(temp)
    return np.diag(diag_heavi)

def optimize_alpha(input_vec, val, factor, K, tol=1e-6):
    # if not twoD, then factor should be None
    alpha_converged = False
    if factor is None:
        factor = np.sqrt(input_vec.shape[0])
    T = np.linalg.inv(dd_chi(G(input_vec, K), val**2))
    dot_product = np.dot(input_vec.T, np.dot(T, input_vec))
    ans = dot_product * factor
    ans = ans / np.linalg.norm(input_vec) / dot_product
    if abs(ans - val**2) <= tol:
        logger.debug(strm("ALPHA HAS CONVERGED."))
        alpha_converged = True
        return ans, alpha_converged
    return ans, alpha_converged

def newton_min(input_vec, val, data_fornnls, K):
    fder = dd_chi(G(input_vec, K), val)
    fval = d_chi(input_vec, val, data_fornnls, K)
    return input_vec + np.dot(np.linalg.inv(fder), fval)

def mod_BRD(guess, K, factor, data_fornnls, maxiter=20):
    smoothing_param = guess
    alpha_converged = False
    for iter in range(maxiter):
        logger.debug(strm("ITERATION NO.", iter))
        logger.debug(strm("CURRENT LAMBDA", smoothing_param))
        retval, residual = this_nnls.nnls_regularized(
            K, data_fornnls, l=smoothing_param
        )
        f_vec = retval[:, newaxis]
        alpha = smoothing_param**2
        c_vec = np.dot(K, f_vec) - data_fornnls[:, newaxis]
        c_vec /= -1 * alpha
        c_update = newton_min(c_vec, smoothing_param, data_fornnls, K)
        alpha_update, alpha_converged = optimize_alpha(
            c_update, smoothing_param, factor, K
        )
        lambda_update = np.sqrt(alpha_update[0, 0])
        if alpha_converged:
            logger.debug(strm("*** OPTIMIZED LAMBDA", lambda_update, "***"))
            break
        if not alpha_converged:
            logger.debug(strm("UPDATED LAMBDA", lambda_update))
            smoothing_param = lambda_update
        if iter == maxiter - 1:
            logger.debug(strm("DID NOT CONVERGE."))
    return lambda_update

def demand_real(x, addtxt=""):
    if not x.dtype == np.float64:
        if x.dtype == np.complex128:
            raise ValueError(
                "you are not allows to pass nnls complex data:\nif it makes sense for you, try yourdata.real.nnls( np.where you now have yourdata.nnls("
                + "\n"
                + addtxt
            )
        else:
            raise ValueError(
                "I expect double-precision floating point (float64), but you passed me data of dtype "
                + str(x.dtype)
                + "\n"
                + addtxt
            )
# }}}
def nnls(self, dimname, newaxis_dict, kernel_func, l=0, default_cut=1e-3):
    r"""Perform regularized non-negative least-squares "fit" on self.

    Capable of solving for solution in 1 or 2 dimensions.

    We seek to minimize
    :math:`Q = \| Ax - b \|_2 + \|\lambda x\|_2`
    in order to obtain solution vector :math:`x` subject to non-negativity constraint
    given input matrix :math:`A`, the kernel, and input vector :math:`b`, the data.

    The first term assesses agreement between the fit :math:`Ax` and the data :math:`b`,
    and the second term accounts for noise with the regularization parameter :math:`\lambda`
    according to Tikhonov regularization.

    To perform regularized minimization in 1 dimension, provide
    :str:`dimname`, :nddata:`newaxis_dict`, :function:`kernel_func`, and
    regularization parameter `l`. One may set `l` to a :double: of the regularization
    parameter of choice (found, for instance, through L-curve analysis) or
    set `l` to :str:`BRD` to enable automatic selection of a regularization
    parameter via the BRD algorithm - namely that described in Venkataramanan et al. 2002
    but adapted for 1D case (DOI:10.1109/78.995059).

    To perform regularized minimization in 2 dimensions, set `l` to
    :str:`BRD` and provide a tuple of parameters :str:`dimname`,
    :nddata:`newaxis_dict`, and :function:`kernel_func`.  Algorithm
    described in Venkataramanan et al. 2002 is performed which determines
    optimal :math:`\lambda` for the data (DOI:10.1109/78.995059).
    Note that setting `l` to a :double: for a regularization
    parameter is supported in this 2 dimensional
    should an appropriate parameter be known.

    See: `Wikipedia page on NNLS <https://en.wikipedia.org/wiki/Non-negative_least_squares>`_,
    `Wikipedia page on Tikhonov regularization <https://en.wikipedia.org/wiki/Tikhonov_regularization>`_

    Parameters
    ==========
    dimname: str or tuple
        Name of the "data" dimension that is to be replaced by a
        distribution (the "fit" dimension);
        *e.g.* if you are regularizing a set of functions
        :math:`\exp(-\tau*R_1)`, then this is :math:`\tau`

        If you are performing 2D regularization, then this
        is a tuple (pair) of 2 names
    newaxis_dict: dict or (tuple of) nddata
        a dictionary whose key is the name of the "fit" dimension
        (:math:`R_1` in the example above)
        and whose value is an np.array with the new axis labels.

        OR
=======
        :math:`\exp(-\tau*R_1)`, then this is :math:`\tau`

        If you are performing 2D regularization, then this
        is a tuple (pair) of 2 names
    newaxis_dict: dict or (tuple of) nddata
        a dictionary whose key is the name of the "fit" dimension
        (:math:`R_1` in the example above)
        and whose value is an np.array with the new axis labels.

        OR

>>>>>>> add_mult_debug
        this can be a 1D nddata
        -- if it has an axis, the axis will be used to create the
        fit axis; if it has no axis, the data will be used

        OR

        if dimname is a tuple of 2 dimensions indicating a 2D ILT, this
        should also be a tuple of 2 nddata, representing the two axes
    kernel_func: function or tuple of functions
        a function giving the kernel for the regularization.
        The first argument is the "data" variable
        and the second argument is the "fit" variable
        (in the example above, this would be something like
        ``lambda x,y: exp(-x*y)``)

        For 2D, this must be a tuple or dictionary of functions -- the kernel is
        the product of the two.
    l : double (default 0) or str
        the regularization parameter :math:`lambda` -- if this is set to 0,
        the algorithm reverts to standard nnls.  If this is set to
        :str:`BRD`, then automatic parameter selection is executed
        according to the BRD algorithm, either in 1-dimension or
        2-dimensions depending on presence of tuple synax (i.e., specifying
        more than 1 dimension).

    Returns
    =======
    self:
        The regularized result.
        For future use, both the kernel (as an nddata, in a property called
        "nnls_kernel") and the residual (as an nddata, in a property called
        "nnls_residual") are stored as properties of the nddata.
        The regularized dimension is always last
        (innermost).

        If the tuple syntax is used to input 2 dimensions and
        :str:`BRD` is specified, then the individual,
        uncompressed kernels :math:`K_{1}` and :math:`K_{2}`
        are returned as properties of the nddata "K1" and "K2"
        respectively. The number of singular values used to
        compressed each kernel is returned in properties of the
        nddata called, respectively, "s1" and "s2".
    """
    # {{{ type checking
    demand_real(self.data)
    if type(dimname) is str and type(newaxis_dict) is type(self):
        # one dimensional with axis given as nddata
        dimname = [dimname]
        newaxis_dict = [newaxis_dict]
    elif type(dimname) is tuple and type(newaxis_dict) is tuple:
        # 2D as tuples
        assert len(dimname) == len(newaxis_dict)
        assert len(dimname) in [1, 2]
    elif type(newaxis_dict) is dict:
        # as a dictionary
        if len(newaxis_dict) == 1 and type(dimname) is str:
            newaxis_dict = [newaxis_dict[dimname]]
        else:
            # here I'm expecting 2D
            assert len(newaxis_dict) == 2
            assert len(dimname) == 2
            newaxis_dict = [newaxis_dict[j] for j in dimname]
    else:
        raise ValueError(
            strm(
                "I didn't understand what you specified for the new axes (dimension:",
                dimname,
                "and new axes",
                newaxis_dict,
            )
        )
    # make sure axes are real
    for j in dimname:
        demand_real(
            self.getaxis(j), "(this message pertains to the %s axis)" % j
        )
    for j in newaxis_dict:
        assert (
            len(j.dimlabels) == 1
        ), "must be one-dimensional, has dimensions:" + str(j.dimlabels)
        if j.getaxis(j.dimlabels[0]) is not None:
            demand_real(
                j.getaxis(j.dimlabels[0]),
                "(this message pertains to the new %s axis pulled from the second argument's axis)"
                % str(j.dimlabels[0]),
            )
        demand_real(
            j.data,
            "(this message pertains to the new %s axis pulled from the second argument's data)"
            % str(j.dimlabels[0]),
        )
    if isinstance(kernel_func, tuple):
        assert callable(kernel_func[0]) and callable(
            kernel_func[1]
        ), "third argument is tuple of kernel functions"
    elif isinstance(kernel_func, dict):
        kernel_func = [kernel_func[j] for j in dimname]
        assert all([callable(j) for j in kernel_func])
    else:
        assert callable(kernel_func), "third argument is kernel function"
        kernel_func = [kernel_func]
    assert len(kernel_func) == len(dimname)
    # }}}
    logger.debug(
        strm(
            "on first calling nnls, shape of the data is",
            self.shape,
            "is it fortran ordered?",
            np.isfortran(self.data),
        )
    )
    # at this point kernel_func and newaxis_dict are both lists with length
    # equal to dimnames (length 1 for 1D and 2 for 2D)
    twoD = len(dimname) > 1
    # construct the kernel
    # the kernel transforms from (columns) the "fit" dimension to (rows)
    # the "data" dimension
    fit_dimnames = [j.dimlabels[0] for j in newaxis_dict]
    fit_axis_coords = [
        j.getaxis(fit_dimnames[j_idx])
        for j_idx, j in enumerate(newaxis_dict)
        if j.getaxis(fit_dimnames[j_idx]) is not None
    ]
    fit_axes = [
        self.__class__(fit_axis_coords[j], fit_dimnames[j])
        for j in range(len(dimname))
    ]
    data_axes = [self.fromaxis(dimname[j]) for j in range(len(dimname))]
    # at this point, fit_axes and data_axes
    # are nddata objects that give the axis
    # coordinates in the fit and data
    # domains, respectively
    for j in range(len(dimname)):
        data_axes[j].squeeze()
        data_axes[j], fit_axes[j] = data_axes[j].aligndata(fit_axes[j])
    else:
        K = S[0] @ V[0]
        data_fornnls = U[0].T @ self.data
        logger.debug(strm(np.shape(K)))
        logger.debug(strm(np.shape(data_fornnls)))
    if len(data_fornnls.shape) > 2:
        logger.debug(strm("Reshaping data.."))
        data_fornnls = data_fornnls.reshape(
            (np.prod(data_fornnls.shape[:-1]), data_fornnls.shape[-1])
        )
    logger.debug(
        strm(
            "np.shape of the data is",
            self.shape,
            "len of axis_coords_error",
            len(self.axis_coords_error),
        )
    )
    if l == "BRD":
        if twoD:
            factor = np.sqrt(s[0]*s[1])
        else:
            factor = None
        retval, residual = this_nnls.nnls_regularized(
            K, data_fornnls, l=mod_BRD(1.0, K, factor, data_fornnls)
        )
    else:
        retval, residual = this_nnls.nnls_regularized(K, data_fornnls, l=l)
    logger.debug(
        strm("coming back from fortran, residual type is", type(residual))
        + strm(residual.dtype if isinstance(residual, np.ndarray) else "")
    )
    newshape = []
    if not np.isscalar(l):
        newshape.append(len(l))
    logger.debug(strm("test***", list(self.data.shape)[:-1]))
    newshape.append(fit_axes[0].shape[fit_dimnames[0]])
    if twoD:
        newshape.append(fit_axes[1].shape[fit_dimnames[1]])
    logger.debug(
        strm(
            "before mkd, shape of the data is",
            self.shape,
            "len of axis_coords_error",
            len(self.axis_coords_error),
        )
    )
    # {{{ store the dictionaries for later use
    axis_coords_dict = self.mkd(self.axis_coords)
    axis_units_dict = self.mkd(self.axis_coords_units)
    axis_coords_error_dict = self.mkd(self.axis_coords_error)
    # }}}
    retval = retval.reshape(newshape)
    self.data = retval
    # {{{ clear axis info
    self.axis_coords = None
    self.axis_coords_units = None
    self.axis_coords_error_dict = None
    # }}}
    # change the dimnesion names and data
    for j in range(len(dimname)):
        self.rename(dimname[j], fit_dimnames[j])
        axis_coords_dict[fit_dimnames[j]] = fit_axes[j].getaxis(
            fit_dimnames[j]
        )
        axis_units_dict[fit_dimnames[j]] = None
        axis_coords_error_dict[fit_dimnames[j]] = None
    if not np.isscalar(l):
        self.dimlabels = ["lambda"] + self.dimlabels
        axis_coords_dict["lambda"] = l
        axis_units_dict["lambda"] = None
        axis_coords_error_dict["lambda"] = None
    self.data = retval
    if not np.isscalar(residual):
        # make the residual nddata as well
        residual_nddata = self.shape.pop(fit_dimnames[0])
        if twoD:
            residual_nddata.pop(fit_dimnames[1])
        residual_nddata = residual_nddata.alloc(dtype=residual.dtype)
        residual_nddata.data[:] = residual[:]
    else:
        residual_nddata = residual
    # store the kernel and the residual as properties
    self.set_prop("nnls_kernel", K)
    self.set_prop("s1", s[0])
    self.set_prop("nnls_residual", residual_nddata)
    for j in range(len(dimname)):
        self.set_prop(f"K{j+1}", K[j])
    if twoD:
        self.set_prop("s2", s[1])
    self.axis_coords = self.fld(axis_coords_dict)
    self.axis_coords_units = self.fld(axis_units_dict)
    self.axis_coords_error = self.fld(axis_coords_error_dict)
    return self
