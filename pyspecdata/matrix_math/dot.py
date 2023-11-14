#@profile
def dot(self,arg):
    """Tensor dot of self with arg -- dot all matching dimension labels.  This can be used to do matrix multiplication, but note that the order of doesn't matter, since the dimensions that are contracted are determined by matching the dimension names, not the order of the dimension.

    >>> a = nddata(r_[0:9],[3,3],['a','b'])
    >>> b = nddata(r_[0:3],'b')
    >>> print a.C.dot(b)
    >>> print a.data.dot(b.data)
    >>> a = nddata(r_[0:27],[3,3,3],['a','b','c'])
    >>> b = nddata(r_[0:9],[3,3],['a','b'])
    >>> print a.C.dot(b)
    >>> print np.tensordot(a.data,b.data,axes=((0,1),(0,1)))

    >>> a = nddata(r_[0:27],[3,3,3],['a','b','c'])
    >>> b = nddata(r_[0:9],[3,3],['a','d'])
    >>> print a.C.dot(b)
    >>> print np.tensordot(a.data,b.data,axes=((0),(0)))
    """
    A,B = self.aligndata(arg)
    matching_dims = list(set(self.dimlabels) & set(arg.dimlabels))
    assert len(matching_dims) > 0, "no matching dimensions!"
    # {{{ store the dictionaries for later use
    axis_coords_dict = A.mkd(A.axis_coords)
    axis_units_dict = A.mkd(A.axis_coords_units)
    axis_coords_error_dict = A.mkd(A.axis_coords_error)
    # }}}
    # manipulate "self" directly
    self.dimlabels = [j for j in A.dimlabels if j not in matching_dims]
    match_idx = [A.axn(j) for j in matching_dims]
    if (self.get_error() is not None) or (arg.get_error() is not None):
        raise ValueError("we plan to include error propagation here, but not yet provided")
    self.data = np.tensordot(A.data,B.data,axes=(match_idx,match_idx))
    logger.debug(strm("shape of A is",ndshape(A)))
    logger.debug(strm("shape of B is",ndshape(B)))
    logger.debug(strm("matching_dims are",matching_dims))
    newsize = [(A.data.shape[j] if A.data.shape[j] != 1 else B.data.shape[j])
            for j in range(len(A.data.shape)) if A.dimlabels[j] not in matching_dims]
    self.data = self.data.reshape(newsize)
    # {{{ use the dictionaries to reconstruct the metadata
    self.axis_coords = self.fld(axis_coords_dict)
    self.axis_coords_units = self.fld(axis_units_dict)
    self.axis_coords_error = self.fld(axis_coords_error_dict)
    # }}}
    return self
