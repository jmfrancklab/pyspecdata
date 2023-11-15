import logging
import numpy as np
from ..general_functions import strm
logger = logging.getLogger('pyspecdata.matrix_math')
#@profile
def along(self,dimname):
    """Specifies the dimension for the next matrix
    multiplication (represents the rows/columns)."""
    self._matmul_along = dimname
    return self
def dot(self,arg):
    """This will perform a dot product or a matrix multiplication.
    If one dimension in ``arg`` matches that in ``self``,
    it will dot along that dimension
    (take a matrix multiplication where that dimension represents the columns of ``self`` and the rows of ``arg``)

    Note that if you have your dimensions named "rows" and "columns", this
    will be very confusing, but if you have your dimensions named in terms
    of the vector basis they are defined/live in, this makes sense.

    If there are zero or no matching dimensions, then use
    :func:`~pyspecdata.nddata.along` to specify the dimensions for matrix
    multiplication / dot product.

    >>> a = nddata(r_[0:9],[3,3],['a','b'])
    >>> b = nddata(r_[0:3],'b')
    >>> print a.C.dot(b)
    >>> print a.data.dot(b.data)
    >>> a = nddata(r_[0:27],[3,3,3],['a','b','c'])
    >>> b = nddata(r_[0:9],[3,3],['a','b'])
    >>> print a.C.dot(b)
    >>> print tensordot(a.data,b.data,axes=((0,1),(0,1)))

    >>> a = nddata(r_[0:27],[3,3,3],['a','b','c'])
    >>> b = nddata(r_[0:9],[3,3],['a','d'])
    >>> print a.C.dot(b)
    >>> print tensordot(a.data,b.data,axes=((0),(0)))

    .. literalinclude:: ../examples/matrix_mult.py
    """
    if hasattr(self,'_matmul_along'):
        dot_dim_A = self._matmul_along
        if hasattr(arg,'_matmul_along'):
            dot_dim_B = arg._matmul_along
        else:
            dot_dim_B = dot_dim_A
    elif hasattr(arg,'_matmul_along'):
        dot_dim_B = arg._matmul_along
        dot_dim_A = dot_dim_B
    else:
        matching_dims = set(self.dimlabels) & set(arg.dimlabels)
        if len(matching_dims) == 1:
            dot_dim_A = list(matching_dims)[0]
            dot_dim_B = dot_dim_A
        else:
            raise ValueError("I can't determine the dimension for matrix"
                    " multiplication!  Either have only one matching dimension,"
                    " or use .along(")
    # {{{ unset the "along" setting
    if hasattr(self,'_matmul_along'):
        del(self._matmul_along)
    if hasattr(arg,'_matmul_along'):
        del(arg._matmul_along)
    # }}}
    if dot_dim_B != dot_dim_A:
        arg = arg.C.rename(dot_dim_B,dot_dim_A)
    A,B = self.aligndata(arg)
    # {{{ store the dictionaries for later use
    axis_coords_dict = A.mkd(A.axis_coords)
    axis_units_dict = A.mkd(A.axis_coords_units)
    axis_coords_error_dict = A.mkd(A.axis_coords_error)
    # }}}
    # manipulate "self" directly
    self.dimlabels = [j for j in A.dimlabels if j != dot_dim_A]
    if (self.get_error() is not None) or (arg.get_error() is not None):
        raise ValueError("we plan to include error propagation here, but not yet provided")
    ax_idx = A.axn(dot_dim_A) # should be the same for both, since they are aligned
    logger.debug(strm("shape of A is",ndshape(A)))
    logger.debug(strm("shape of B is",ndshape(B)))
    A.data = A.data[...,newaxis]
    B.data = B.data[...,newaxis]
    logger.debug(strm("add a new axis",A.data.shape))
    logger.debug(strm("add a new axis",B.data.shape))
    A.data = np.moveaxis(A.data,ax_idx,-1) # columns position
    B.data = np.moveaxis(B.data,ax_idx,-2) # row position
    logger.debug(strm("after movement",A.data.shape))
    logger.debug(strm("after movement",B.data.shape))
    self.data = np.matmul(A.data,B.data)
    logger.debug(strm("after mult",self.data.shape))
    self.data = self.data[...,0,0]
    logger.debug(strm("remove extras",self.data.shape))
    # {{{ use the dictionaries to reconstruct the metadata
    self.axis_coords = self.fld(axis_coords_dict)
    self.axis_coords_units = self.fld(axis_units_dict)
    self.axis_coords_error = self.fld(axis_coords_error_dict)
    # }}}
    print("self.data.shape",self.data.shape)
    print("ndshape",ndshape(self))
    return self
#@profile
def matmul(self,arg):
    assert type(arg) is nddata, "currently matrix multiplication only allowed if both are nddata"
    return self.C.dot(arg)
