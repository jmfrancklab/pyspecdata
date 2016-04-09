from numpy import *
def secsy_transform_manual(self,
        direct_dim,
        indirect_dim,
        has_indirect = True,
        truncate = False):
    r'''Shift the time-domain data backwards by the echo time.
    As opposed to :func:`secsy_transform <pyspecdata.axis_manipulation.secsy_transform>`, this directlly manipulates the phase of the function, rather than calling on :func:`skew <pyspecdata.skew>`.

    Parameters
    ----------
    has_indirect : bool
        (This option is largely specific to data loaded by :func:`acert_hdf5 <pyspecdata.acert_hdf5>`)

        Does the data actually have an indirect dimension?
        If not, assume that there is a constant echo time,
        that can be retrieved with ``.get_prop('te')``.
    truncate : bool
        If this is set, `register_axis <pyspecdata.axis_manipulation.register_axis>` to :math:`t_{direct}=0`,
        and then throw out the data for which :math:`t_{direct}<0`.
    '''
    if has_indirect:
        self *= self.fromaxis([indirect_dim,direct_dim],lambda t1,t2: exp(1j*2*pi*t1*t2))
    else:
        echo_time = self.get_prop('te')
        self *= self.fromaxis(direct_dim,lambda t2: exp(1j*2*pi*echo_time*t2)) # positive time shift corrects positive slope
    if truncate:
        raise ValueError("truncate is not yet supported")
    return self
def secsy_transform(self,
        direct_dim,
        indirect_dim,
        has_indirect = True,
        method = 'fourier',
        truncate = True):
    r'''Shift the time-domain data backwards by the echo time.

    As opposed to :func:`secsy_transform_manual <pyspecdata.axis_manipulation.secsy_transform_manual>`, this calls on on :func:`skew <pyspecdata.skew>`,
    rather than directly manipulating the phase of the function, which can lead to aliasing.

    Parameters
    ----------
    has_indirect : bool
        (This option is largely specific to data loaded by :func:`acert_hdf5 <pyspecdata.acert_hdf5>`)

        Does the data actually have an indirect dimension?
        If not, assume that there is a constant echo time,
        that can be retrieved with ``.get_prop('te')``.
    truncate : bool
        If this is set, `register_axis <pyspecdata.axis_manipulation.register_axis>` to :math:`t_{direct}=0`,
        and then throw out the data for which :math:`t_{direct}<0`.
    method : str
        The shear method (linear or fourier).
    '''
    if has_indirect:
        self.shear(direct_dim,indirect_dim,-1,method = method)
    else:
        echo_time = self.get_prop('te')
        self.extend(direct_dim,self.getaxis(direct_dim)[0]-echo_time)
        self *= self.fromaxis(direct_dim,lambda t2: exp(1j*2*pi*echo_time*t2)) # positive time shift corrects positive slope
    if truncate:
        for thisaxis in [direct_dim,indirect_dim]:
            self.register_axis({thisaxis:0.})
            #self = self[thisaxis:(0.,)] 
            keep_zeros = False # just for debugging
            if keep_zeros:
                self[thisaxis,
                        lambda x: x<0] = 0.
            else:
                newdata = self[thisaxis:(0.,)] 
                self.data = newdata.data
                self.setaxis(thisaxis,newdata.getaxis(thisaxis))
    return self
