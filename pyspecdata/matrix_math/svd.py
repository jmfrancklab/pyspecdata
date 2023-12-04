import numpy as np
def svd(self, todim, fromdim):
    """Singular value decomposition.  Original matrix is unmodified.

    .. note::
        Because we are planning to upgrade with axis objects,
        FT properties, axis errors, etc, are not transferred here.
        If you are using it when this note is still around, be sure to
        `.copy_props(`

        Also, error, units, are not currently propagated, but could be relatively easily!

    If

    >>> U, Sigma, Vh = thisinstance.svd()

    then ``U``, ``Sigma``, and ``Vh`` are nddata such that ``result`` in

    >>> result = U @ Sigma @ Vh

    will be the same as ``thisinstance``.
    Note that this relies on the fact that nddata matrix multiplication doesn't care about the ordering
    of the dimensions (see :method:`~pyspecdata.core.dot`).
    The vector space that contains the singular values is called `'SV'` (see more below).

    Parameters
    ==========
    fromdim: str
        This dimension corresponds to the columns of the matrix that is
        being analyzed by SVD.
        (The matrix transforms from the vector space labeled by ``fromdim``
        and into the vector space labeled by ``todim``).
    todim: str
        This dimension corresponds to the rows of the matrix that is
        being analyzed by SVD.

    Returns
    =======
    U: nddata
        Has dimensions (all other dimensions) × 'todim' × 'SV',
        where the dimension 'SV' is the vector space of the singular
        values.
    Sigma: nddata
        Has dimensions (all other dimensions) × 'SV'.
        Only non-zero
    Vh: nddata
        Has dimensions (all other dimensions) × 'SV' × 'fromdim',
    """
    orig_order = list(self.dimlabels)
    all_but = [j for j in self.dimlabels if j not in [fromdim,todim]]
    new_order = all_but + [todim,fromdim]
    self.reorder(new_order)
    U, Sigma, Vh = np.linalg.svd(self.data, full_matrices=False)
    U = self.__class__(U,all_but + [todim,'SV'])
    Vh = self.__class__(Vh,all_but + ['SV',fromdim])
    Sigma = self.__class__(Sigma, all_but + ['SV'])
    # {{{ label the axes
    for j in all_but:
        U.setaxis(j,self.getaxis(j))
        Vh.setaxis(j,self.getaxis(j))
    U.setaxis(todim,self.getaxis(todim))
    Vh.setaxis(fromdim,self.getaxis(fromdim))
    # }}}
    self.reorder(orig_order)
    return U, Sigma, Vh
