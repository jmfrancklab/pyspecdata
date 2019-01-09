# again, this is copied liberally from scipy nnls -- see scipy licensing
from __future__ import division, print_function, absolute_import

from . import _nnls
from numpy import asarray_chkfinite, zeros, double, isscalar

__all__ = ['nnls_regularized']


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
    if len(b.shape) != 1:
        raise ValueError("expected vector")

    m, n = A.shape

    if m != b.shape[0]:
        raise ValueError("incompatible dimensions")

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
            x, rnorm, mode = _nnls.nnls_regularized(A, b, w, zz, index, maxiter, l)
    else:
            w = zeros((n,), dtype=double)
            zz = zeros((m+n,), dtype=double)
            index = zeros((n,), dtype=int)
            x, rnorm, mode = _nnls.nnls_regularized_loop(A, b, w, zz, index, maxiter, l)
            # From the documentation, I wouldn't have thought the following is needed, but it does seem to be
            x = x.ravel('F').reshape(x.shape)
    if mode != 1:
        raise RuntimeError("too many iterations")
    return x, rnorm
