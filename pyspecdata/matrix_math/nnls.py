def nnls(self, dimname, newaxis_dict, kernel_func, l=0):
    r"""Perform regularized non-negative least-squares "fit" on self.

    Capable of solving for solution in 1 or 2 dimensions.

    We seek to minimize
    :math:`Q = \| Ax - b \|_2 + \|\lambda x\|_2`
    in order to obtain solution vector :math:`x` subject to non-negativity constraint
    given input matrix :math:`A`, the kernel, and input vector :math:`b`, the data.

    The first term assesses agreement between the fit :math:`Ax` and the data :math:`b`,
    and the second term accounts for noise with the regularization parameter :math:`\lambda`
    according to Tikhonov regularization.

    To perform regularized minimization in 2 dimensions, set `l` to :str:`BRD` and provide a
    tuple of parameters :str:`dimname`, :nddata:`newaxis_dict`, and :function:`kernel_func`.
    Algorithm described in Venkataramanan et al. 2002 is performed which determines optimal :math:`\lambda`
    for the data (DOI:10.1109/78.995059).
    
    See: `Wikipedia page on NNLS <https://en.wikipedia.org/wiki/Non-negative_least_squares>`_,
    `Wikipedia page on Tikhonov regularization <https://en.wikipedia.org/wiki/Tikhonov_regularization>`_
     
    Parameters
    ==========
    dimname: str
        Name of the "data" dimension that is to be replaced by a
        distribution (the "fit" dimension);
        *e.g.* if you are regularizing a set of functions
        :math:`\np.exp(-\tau*R_1)`, then this is :math:`\tau`
    newaxis_dict: dict or nddata
        a dictionary whose key is the name of the "fit" dimension
        (:math:`R_1` in the example above)
        and whose value is an np.array with the new axis labels.
        OR
        this can be a 1D nddata
        -- if it has an axis, the axis will be used to create the
        fit axis; if it has no axis, the data will be used
    kernel_func: function
        a function giving the kernel for the regularization.
        The first argument is the "data" variable
        and the second argument is the "fit" variable
        (in the example above, this would be something like
        ``lambda x,y: np.exp(-x*y)``)
    l : double (default 0) or str
        the regularization parameter
        :math:`lambda` -- if this is set to 0, the algorithm reverts to
        standard nnls.
        If this is set to :str:`BRD`, then algorithm expects tuple of each parameter
        described above in order to perform a 2-dimensional fit.

    Returns
    =======
    self:
        The regularized result.
        For future use, both the kernel (as an nddata, in a property called
        "nnls_kernel") and the residual (as an nddata, in a property called
        "nnls_residual") are stored as properties of the nddata.
        The regularized dimension is always last
        (innermost).
        If :str:`BRD` is specified, then the individual, uncompressed kernels :math:`K_{1}` and :math:`K_{2}` are returned as properties of the nddata "K1" and "K2" respectively. The number of singular values used to compressed each kernel is returned in properties of the nddata called, respectively, "s1" and "s2". 
    """
    logger.debug(strm('on first calling nnls, shape of the data is',ndshape(self),'is it fortran ordered?',np.isfortran(self.data)))
    tuple_syntax = False
    if isinstance(dimname, tuple):
        tuple_syntax = True
        assert len(dimname) == 2, "tuple of two dimension names only"
        assert type(dimname[0]) and isinstance(dimname[1], str), "first argument is tuple of two dimension names"
    else:
        assert isinstance(dimname, str), "first argument is dimension name or tuple of two dimension names"
    if isinstance(newaxis_dict, tuple):
        assert len(newaxis_dict) == 2, "tuple of two nddatas only"
        if isinstance(newaxis_dict[0],nddata) and isinstance(newaxis_dict[1],nddata):
            assert len(newaxis_dict[0].dimlabels) and len(newaxis_dict[1].dimlabels) == 1, "currently only set up for 1D"
    elif isinstance(newaxis_dict, dict):
        assert len(newaxis_dict) == 1, "currently only set up for 1D"
    elif isinstance(newaxis_dict,nddata):
        assert len(newaxis_dict.dimlabels) == 1, "currently only set up for 1D"
    else:
        raise ValueError("second argument is dictionary or nddata with new axis, or tuple of nddatas with new axes")
    if isinstance(kernel_func, tuple):
        assert callable(kernel_func[0]) and callable(kernel_func[1]), "third argument is tuple of kernel functions"
    else:
        assert callable(kernel_func), "third argument is kernel function"
    # construct the kernel
    # the kernel transforms from (columns) the "fit" dimension to (rows)
    # the "data" dimension
    if tuple_syntax:
        if isinstance(newaxis_dict[0],nddata):
            assert len(newaxis_dict[0].dimlabels) and len(newaxis_dict[1].dimlabels) == 1, "must be 1 dimensional!!"
            fitdim_name1 = newaxis_dict[0].dimlabels[0]
            fitdim_name2 = newaxis_dict[1].dimlabels[0]
            fit_axis1 = newaxis_dict[0].getaxis(fitdim_name1)
            fit_axis2 = newaxis_dict[1].getaxis(fitdim_name2)
            if fit_axis1 is None:
                fit_axis1 = newaxis_dict[0].data
            if fit_axis2 is None:
                fit_axis2 = newaxis_dict[1].data
    elif isinstance(newaxis_dict,nddata):
        assert len(newaxis_dict.dimlabels) == 1, "must be 1 dimensional!!"
        fitdim_name = newaxis_dict.dimlabels[0]
        fit_axis = newaxis_dict.getaxis(fitdim_name)
        if fit_axis is None:
            fit_axis = newaxis_dict.data
    else:
        fitdim_name = list(newaxis_dict.keys())[0]
        logger.debug(strm('shape of fit dimension is',newaxis_dict[fitdim_name].shape))
        fit_axis = newaxis_dict[fitdim_name]
    if tuple_syntax:
        fit_axis1 = nddata(fit_axis1,fitdim_name1)
        fit_axis2 = nddata(fit_axis2,fitdim_name2)
        data_axis1 = self.fromaxis(dimname[0])
        data_axis2 = self.fromaxis(dimname[1])
        data_axis1.squeeze()
        data_axis2.squeeze()
        data_axis1,fit_axis1 = data_axis1.aligndata(fit_axis1)
        data_axis2,fit_axis2 = data_axis2.aligndata(fit_axis2)
        K1 = kernel_func[0](data_axis1,fit_axis1).squeeze()
        K1_ret = K1
        logger.debug(strm('K1 dimlabels',K1.dimlabels,'and raw shape',K1.data.shape))
        K2 = kernel_func[1](data_axis2,fit_axis2).squeeze()
        K2_ret = K2
        logger.debug(strm('K2 dimlabels',K2.dimlabels,'and raw shape',K2.data.shape))
        # SVD and truncation of kernels
        U1,S1,V1 = np.linalg.svd(K1.data,full_matrices=False)
        U2,S2,V2 = np.linalg.svd(K2.data,full_matrices=False)
        logger.debug(strm('Uncompressed SVD K1:',[x.shape for x in (U1,S1,V1)]))
        logger.debug(strm('Uncompressed SVD K2:',[x.shape for x in (U2,S2,V2)]))
        default_cut = 1e-2
        s1 = np.where(S1 > default_cut)[0][-1]
        s2 = np.where(S2 > default_cut)[0][-1]
        U1 = U1[:,0:s1]
        S1 = S1[0:s1]
        V1 = V1[0:s1,:]
        U2 = U2[:,0:s2]
        S2 = S2[0:s2]
        V2 = V2[0:s2,:]
        S1 = S1*np.eye(s1)
        S2 = S2*np.eye(s2)
        logger.debug(strm('Compressed SVD of K1:',[x.shape for x in (U1,S1,V1)]))
        logger.debug(strm('Compressed SVD K2:',[x.shape for x in (U2,S2,V2)]))
        # would generate projected data here
        # compress data here
        K1 = S1.dot(V1)
        K2 = S2.dot(V2)
        K = K1[:,newaxis,:,newaxis]*K2[newaxis,:,newaxis,:]
        K = K.reshape(K1.shape[0]*K2.shape[0],K1.shape[1]*K2.shape[1])
        logger.debug(strm('Compressed K0, K1, and K2:',[x.shape for x in (K,K1,K2)]))

        data_compressed = U1.T.dot(self.data.dot(U2))
        logger.debug(strm('Compressed data:',data_compressed.shape))
        # data_for_nnls = nddata(data_compressed,[dimname[0],dimname[1]])
        # data_for_nnls.smoosh([dimname[0],dimname[1]],dimname=dimname[0])
        data_fornnls = np.empty(s1*s2)
        for s1_index in range(s1):
            for s2_index in range(s2):
                temp = data_compressed[s1_index][s2_index]
                data_fornnls[s1_index*s2+s2_index] = temp
        logger.debug(strm('Lexicographically ordered data:',data_fornnls.shape))
        if len(data_fornnls.shape) > 2:
            logger.debug(strm('Reshpaing data..'))
            data_fornnls = data_fornnls.reshape((np.prod(data_fornnls.shape[:-1]),data_fornnls.shape[-1]))
        
        if l == 'BRD':
            def chi(x_vec,val):
                return 0.5*np.dot(x_vec.T,np.dot(dd_chi(G(x_vec),val**2),x_vec)) - np.dot(x_vec.T,data_fornnls[:,newaxis])
            def d_chi(x_vec,val):
                return np.dot(dd_chi(G(x_vec),val**2),x_vec) - data_fornnls[:,newaxis]
            def dd_chi(G,val):
                return G + (val**2)*np.eye(np.shape(G)[0])
            def G(x_vec):
                return np.dot(K,np.dot(square_heaviside(x_vec),K.T))
            def H(product):
                if product <= 0:
                    return 0
                if product > 0:
                    return 1
            def square_heaviside(x_vec):
                diag_heavi = []
                for q in range(np.shape(K.T)[0]):
                    pull_val = np.dot(K.T[q,:],x_vec)
                    temp = pull_val[0]
                    temp = H(temp)
                    diag_heavi.append(temp)
                diag_heavi = np.array(diag_heavi)
                square_heavi = diag_heavi*np.eye(np.shape(diag_heavi)[0])
                return square_heavi
            def optimize_alpha(input_vec,val):
                alpha_converged = False
                factor = sqrt(s1*s2)
                T = np.linalg.inv(dd_chi(G(input_vec),val**2))
                dot_product = np.dot(input_vec.T,np.dot(T,input_vec))
                ans = dot_product*factor
                ans = ans/np.linalg.norm(input_vec)/dot_product
                tol = 1e-3
                if abs(ans-val**2) <= tol:
                    logger.debug(strm('ALPHA HAS CONVERGED.'))
                    alpha_converged = True
                    return ans,alpha_converged
                return ans,alpha_converged
            def newton_min(input_vec,val):
                fder = dd_chi(G(input_vec),val)
                fval = d_chi(input_vec,val)
                return (input_vec + np.dot(np.linalg.inv(fder),fval))
            def mod_BRD(guess,maxiter=20):
                smoothing_param = guess
                alpha_converged = False
                for iter in range(maxiter):
                    logger.debug(strm('ITERATION NO.',iter))
                    logger.debug(strm('CURRENT LAMBDA',smoothing_param))
                    retval,residual = this_nnls.nnls_regularized(K,data_fornnls,l=smoothing_param)
                    f_vec = retval[:,newaxis]
                    alpha = smoothing_param**2
                    c_vec = np.dot(K,f_vec) - data_fornnls[:,newaxis]
                    c_vec /= -1*alpha
                    c_update = newton_min(c_vec,smoothing_param)
                    alpha_update,alpha_converged = optimize_alpha(c_update,smoothing_param)
                    lambda_update = sqrt(alpha_update[0,0])
                    if alpha_converged:
                        logger.debug(strm('*** OPTIMIZED LAMBDA',lambda_update,'***'))
                        break
                    if not alpha_converged:
                        logger.debug(strm('UPDATED LAMBDA',lambda_update))
                        smoothing_param = lambda_update
                    if iter == maxiter-1:
                        logger.debug(strm('DID NOT CONVERGE.'))
                return lambda_update
            retval, residual = this_nnls.nnls_regularized(K,data_fornnls,l=mod_BRD(guess=1.0))
        else:
            retval, residual = this_nnls.nnls_regularized(K,data_fornnls,l=l)
        logger.debug(strm('coming back from fortran, residual type is',type(residual))+ strm(residual.dtype if isinstance(residual, np.ndarray) else ''))
        newshape = []
        if not np.isscalar(l):
            newshape.append(len(l))
        logger.debug(strm('test***',list(self.data.shape)[:-1]))
        #newshape += list(self.data.shape)[:-1] # this would return parametric axis
        newshape.append(ndshape(fit_axis1)[fitdim_name1])
        newshape.append(ndshape(fit_axis2)[fitdim_name2])
        logger.debug(strm('before mkd, shape of the data is',ndshape(self),'len of axis_coords_error',len(self.axis_coords_error)))
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
        self.rename(dimname[0], fitdim_name1)
        self.rename(dimname[1], fitdim_name2)
        axis_coords_dict[fitdim_name1] = fit_axis1.getaxis(fitdim_name1)
        axis_units_dict[fitdim_name1] = None
        axis_coords_error_dict[fitdim_name1] = None
        axis_coords_dict[fitdim_name2] = fit_axis2.getaxis(fitdim_name2)
        axis_units_dict[fitdim_name2] = None
        axis_coords_error_dict[fitdim_name2] = None
        if not np.isscalar(l):
            self.dimlabels = ['lambda'] + self.dimlabels
            axis_coords_dict['lambda'] = l
            axis_units_dict['lambda'] = None
            axis_coords_error_dict['lambda'] = None
        self.data = retval
        if not np.isscalar(residual):
            # make the residual nddata as well
            residual_nddata = ndshape(self).pop(fitdim_name2).pop(fitdim_name1).alloc(dtype=residual.dtype)
            residual_nddata.data[:] = residual[:]
        else:
            residual_nddata = residual
        # store the kernel and the residual as properties
        self.set_prop('nnls_kernel',K)
        self.set_prop('s1',s1)
        self.set_prop('s2',s2)
        self.set_prop('nnls_residual',residual_nddata)
        self.set_prop('K1',K1_ret)
        self.set_prop('K2',K2_ret)
        # {{{ use the info from the dictionaries
        self.axis_coords = self.fld(axis_coords_dict)
        self.axis_coords_units = self.fld(axis_units_dict)
        self.axis_coords_error = self.fld(axis_coords_error_dict)
        return self
    else:
        fit_axis = nddata(fit_axis, fitdim_name)
        data_axis = self.fromaxis(dimname)
        data_axis.squeeze()
        data_axis, fit_axis = data_axis.aligndata(fit_axis)
        K = kernel_func(data_axis, fit_axis).squeeze()
        logger.debug(strm('K dimlabels',K.dimlabels,'and raw shape',K.data.shape))
        self.reorder(dimname, first=False) # make the dimension we will be regularizing innermost
        logger.debug(strm('shape of the data is',ndshape(self),'is it fortran ordered?',np.isfortran(self.data)))
        data_fornnls = self.data
        if len(data_fornnls.shape) > 2:
            data_fornnls = data_fornnls.reshape((np.prod(
                data_fornnls.shape[:-1]),data_fornnls.shape[-1]))
        logger.debug(strm('shape of the data is',ndshape(self),"len of axis_coords_error",len(self.axis_coords_error)))
        retval, residual = this_nnls.nnls_regularized(K.data, data_fornnls, l=l)
        logger.debug(strm("coming back from fortran, residual type is",type(residual))+ strm(residual.dtype if isinstance(residual, np.ndarray) else ''))
        newshape = []
        if not np.isscalar(l):
            newshape.append(len(l))
        newshape += list(self.data.shape)[:-1] # exclude data dimension
        newshape.append(ndshape(fit_axis)[fitdim_name])
        logger.debug(strm('before mkd, shape of the data is',ndshape(self),"len of axis_coords_error",len(self.axis_coords_error)))
        # {{{ store the dictionaries for later use
        axis_coords_dict = self.mkd(self.axis_coords)
        axis_units_dict = self.mkd(self.axis_coords_units)
        axis_coords_error_dict = self.mkd(self.axis_coords_error)
        # }}}
        retval = retval.reshape(newshape)
        self.data = retval
        # {{{ clear all the axis info
        self.axis_coords = None
        self.axis_coords_units = None
        self.axis_coords_error_dict = None
        # }}}
        # change the dimension names and data
        self.rename(dimname, fitdim_name)
        # {{{ manipulate the dictionaries, and call fld below
        axis_coords_dict[fitdim_name] = fit_axis.getaxis(fitdim_name)
        axis_units_dict[fitdim_name] = None
        axis_coords_error_dict[fitdim_name] = None
        if not np.isscalar(l):
            self.dimlabels = ['lambda'] + self.dimlabels
            axis_coords_dict['lambda'] = l
            axis_units_dict['lambda'] = None
            axis_coords_error_dict['lambda'] = None
        # }}}
        self.data = retval
        if not np.isscalar(residual):
            # make the residual nddata as well
            residual_nddata = ndshape(self).pop(fitdim_name).alloc(dtype=residual.dtype)
            residual_nddata.data[:] = residual[:]
        else:
            residual_nddata = residual
        # store the kernel and the residual in the properties
        self.set_prop('nnls_kernel',K)
        self.set_prop('nnls_residual',residual_nddata)
        # {{{ use the info from the dictionaries
        self.axis_coords = self.fld(axis_coords_dict)
        self.axis_coords_units = self.fld(axis_units_dict)
        self.axis_coords_error = self.fld(axis_coords_error_dict)
        # }}}
        return self
