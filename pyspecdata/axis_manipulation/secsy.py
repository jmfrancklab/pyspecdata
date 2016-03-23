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
        self *= self.fromaxis(['t1','t2'],lambda t1,t2: exp(1j*2*pi*t1*t2))
    else:
        echo_time = self.get_prop('te')
        self *= self.fromaxis('t2',lambda t2: exp(1j*2*pi*echo_time*t2)) # positive time shift corrects positive slope
    if truncate:
        raise ValueError("truncate is not yet supported")
    return self
