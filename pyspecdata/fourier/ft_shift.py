"shift-related helper functions"
from numpy import r_
import numpy as np
import logging
from ..general_functions import *
thinkaboutit_message = ("If you think about it, you"
                        " probably don't want to do this.  You either want to fill with"
                        " zeros from zero up to the start or you want to first set the"
                        " start point to zero.")
def ft_state_to_str(self,*axes):
    """Return a string that lists the FT domain for the given axes.
    
    :math:`u` refers to the original domain (typically time) and :math:`v` refers to the FT'd domain (typically frequency)
    If no axes are passed as arguments, it does this for all axes."""
    retstr = []
    if len(axes) == 0:
        axes = self.dimlabels
    for j in axes:
        if self.get_ft_prop(j):
            retstr.append(j+' --> v (typ. freq) domain')
        else:
            retstr.append(j+' --> u (typ. time) domain')
    return '\n'.join(retstr)
def get_ft_prop(self,axis,propname = None):
    ("Gets the FT property given by `propname`.  For both setting"
            " and getting, `None` is equivalent to an unset value"
            " if no `propname` is given, this just sets the"
            " `FT` property, which tells if a dimension is"
            " frequency or time domain")
    if propname is None:
        propname = []
    elif isinstance(propname, str):
        propname = [propname]
    key_name = '_'.join(['FT'] + propname)
    this_dict = self.get_prop(key_name)
    if this_dict is None:
        return None
    elif axis in list(this_dict.keys()):
        return this_dict[axis]
    else:
        return None
def set_ft_prop(self,axis,propname = None,value = True):
    ("Sets the FT property given by `propname`.  For both setting"
            " and getting, `None` is equivalent to an unset value"
            " if `propname` is a boolean, and `value` is True"
            " (the default), it's assumed that propname is"
            " actually None, and that `value` is set to the"
            " `propname` argument (this allows us to set the"
            " `FT` property more easily)")
    #{{{ if I pass no propname, but a value
    if isinstance(propname, bool) and value is True:
        value = propname
        propname = None
    #}}}
    if propname is None:
        propname = []
    elif isinstance(propname, str):
        propname = [propname]
    key_name = '_'.join(['FT'] + propname)
    this_dict = self.get_prop(key_name)
    if value is None:# unset
        if this_dict is not None and axis in list(this_dict.keys()):
            this_dict.pop(axis)
            if len(this_dict) == 0:
                self.unset_prop(key_name)
    else:
        if this_dict is None:
            self.set_prop(key_name,{axis:value})
        else:
            this_dict[axis] = value
    return self# important, so that I can chain operations
def _ft_shift(self,thisaxis,p2,shift_axis = None,verbose = False):
    ("perform a generalized fftshift along the axis indicated by the integer `thisaxis`, where `p2` gives the index that will become the first index"
    "\n this is derived from the numpy fftshift routine, but defines slices instead of index numbers"
    "\n `shift_axis` is only used after the (i)fft.  It assumes that the axis labels start at zero, and it aliases them over in the same way the data was aliased")
    if p2 == 0: # nothing to be done
        return self
    newdata = np.empty_like(self.data)
    n = self.data.shape[thisaxis]
    sourceslice = [slice(None,None,None)] * len(self.data.shape)
    targetslice = [slice(None,None,None)] * len(self.data.shape)
    # move second half first -- the following are analogous to the numpy function, but uses slices instead
    sourceslice[thisaxis] = slice(p2,n)
    targetslice[thisaxis] = slice(None,n-p2)
    newdata[tuple(targetslice)]  = self.data[tuple(sourceslice)]
    # move first half second (the negative frequencies)
    sourceslice[thisaxis] = slice(None,p2)
    targetslice[thisaxis] = slice(-p2,None)
    newdata[tuple(targetslice)]  = self.data[tuple(sourceslice)]
    self.data = newdata
    if shift_axis is not None and shift_axis:
        axisname = self.dimlabels[thisaxis]
        x = self.getaxis(axisname)
        newaxis = np.empty_like(x)
        n = len(x)
        # move second half first -- the following are analogous to the numpy function, but uses slices instead
        sourceslice = slice(p2,n)
        targetslice = slice(None,n-p2)
        assert x[0] == 0.
        newaxis[targetslice]  = x[sourceslice] - x[n-1] - x[1] # when I alias back over, I can't subtract x[n-1], which would give 0, so I subtract one more dx (note the assertion above)
        # move first half second (the negative frequencies)
        sourceslice = slice(None,p2)
        targetslice = slice(-p2,None)
        newaxis[targetslice]  = x[sourceslice]
        self.setaxis(axisname,newaxis)
    return self
def ft_clear_startpoints(self,axis,t=None,f=None,nearest=None):
    r"""Clears memory of where the origins in the time and frequency domain are.
    This is useful, *e.g.* when you want to ift and center about time=0.
    By setting shift=True you can also manually set the points.

    Parameters
    ==========
    t: float, 'current', 'reset', or None
        keyword arguments `t` and `f` can be set by (1) manually setting
        the start point (2) using the string 'current' to leave the
        current setting a lone (3) 'reset', which clears the startpoint
        and (4) None, which will be changed to 'current' when the other is set to a number or 'rest' if both are set to None.
    t: float, 'current', 'reset', or None
        see `t`
    nearest: bool
        Shifting the startpoint can only be done
        by an integral number of datapoints
        (*i.e.* an integral number of dwell
        times, dt, in the time domain or
        integral number of df in the frequency
        domain).
        While it is possible to shift by a
        non-integral number of datapoints,
        this is done by applying a
        phase-dependent shift in the inverse
        domain.
        Applying such a axis-dependent shift
        can have vary unexpected effects if the
        data in the inverse domain is aliased,
        and is therefore heavily discouraged.
        (For example, consider what happens if
        we attempt to apply a
        frequency-dependent phase shift to data
        where a peak at 110 Hz is aliased and
        appears at the 10 Hz position.)

        Setting `nearest` to **True**
        will choose a startpoint at the closest
        integral datapoint to what you have
        specified.

        Setting `nearest` to **False**
        will explicitly override the safeties --
        essentially telling the code that you
        know the data is not aliased in the
        inverse domain and/or are willing to
        deal with the consequences.
    """
    if f is None and t is None:
        t='reset'
        f='reset'
    if f is None:
        assert t != 'current'
        f='current'
    elif t is None:
        assert f != 'current'
        t='current'
    if f == 'current':
        pass
    else:
        if f == 'reset':
            f=None
        df = _get_ft_df(self,axis)
        orig_f = self.get_ft_prop(axis,['start','freq'])
        if orig_f is None and self.get_ft_prop(axis):
            orig_f = self.getaxis(axis)[0]
        if f is not None:
            n_df = (orig_f-f)/df # number of df's shifted by
            if abs((n_df - round(n_df))/n_df) > 1e-3:
                if nearest is None:
                    raise ValueError(strm("You need to explicitly"
                        " set `nearest`, since you are trying to shift"
                        " the start point from",orig_f,"to",f,
                        "which is a non-integral number ",abs(orig_f-f)/df,
                        " of df=",df,"intervals (n_df ", n_df,
                        ".  If you don't know why"
                        " you're getting this error, see the documentation"
                        " for ft_clear_startpoints!!"))
                elif nearest:
                    f = round((f-orig_f)/df)*df + orig_f
            else:
                f = orig_f - round(n_df)*df
        self.set_ft_prop(axis,['start_freq'],f)
        self.set_ft_prop(axis,['freq','not','aliased'],None)
        if nearest is False:
            self.set_ft_prop(axis,['time','not','aliased'],True)
    if t == 'current':
        pass
    else:
        if t == 'reset':
            t=None
        dt = _get_ft_dt(self,axis)
        orig_t = self.get_ft_prop(axis,['start','time'])
        if orig_t is None and not self.get_ft_prop(axis):
            orig_t = self.getaxis(axis)[0]
        if t is not None:
            n_dt = (orig_t-t)/dt # number of dt's shifted by
            logging.debug(strm("trying to shift by",n_dt))
            if n_dt != 0 and abs((n_dt - round(n_dt))/n_dt) > 1e-3:
                if nearest is None:
                    logging.debug(strm("discrepancy",abs(orig_t-t) % dt))
                    raise ValueError(strm("You need to explicitly"
                        " set `nearest`, since you are trying to shift"
                        " the start point from",orig_t,"to",t,
                        "which is a non-integral number ",abs(orig_t-t)/dt,
                        " of dt=",dt,"intervals (n_dt ", n_dt,
                        ".  If you don't know why"
                        " you're getting this error, see the documentation"
                        " for ft_clear_startpoints!!"))
                elif nearest:
                    t = round((t-orig_t)/dt)*dt + orig_t
                    logging.debug(strm("nearest t is",t))
            else:
                t = orig_t - round(n_dt)*dt
                logging.debug(strm("setting t to",t))
        self.set_ft_prop(axis,['start_time'],t)
        self.set_ft_prop(axis,['time','not','aliased'],None)
        if nearest is False:
            self.set_ft_prop(axis,['freq','not','aliased'],True)
    return self
def _get_ft_df(self,axis):
    if self.get_ft_prop(axis):
        # axis is in frequency domain
        return np.diff(self.getaxis(axis)[r_[0,1]]).item()
    else:
        # axis in time domain
        N = len(self.getaxis(axis))
        return (N)/(N+1)/np.diff(self.getaxis(axis)[r_[0,-1]]).item()
def _get_ft_dt(self,axis):
    if self.get_ft_prop(axis):
        # axis is in frequency domain
        N = len(self.getaxis(axis))
        return (N)/(N+1)/np.diff(self.getaxis(axis)[r_[0,-1]]).item()
    else:
        # axis in time domain
        return np.diff(self.getaxis(axis)[r_[0,1]]).item()
def _find_index(u,origin = 0.0,tolerance = 1e-4,verbose = False):
    ("identify the index of `u` (represents either time or frequency) where"
            " `origin` lives -- if it finds a value exactly equal"
            " to `origin`, returns `(p2,None)` -- otherwise, `None` is replaced by"
            " u-shift indicating how far the position"
            " marked by `p2` is ahead of `origin`"
            "\n\tIf `origin` is outside the range of `u`, it assumes that you want"
            " to find the appropriate index in one of the higher or lower-frequency"
            " replicates (*i.e.*, identical to within $n\\timesSW$, with $n$ integer")
    du = check_ascending_axis(u,tolerance,"In order to determine the FT shift index")
    N = len(u)
    SW = du*N # should be u[-1]+du
    alias_number = 0
    if not u[0] < origin < u[-1]:
        alias_number = (origin - u[0]) // SW # subtracting this many SW's from origin
        #                                     will land me back inside the range of u
        if verbose: logging.debug(strm("(_find_index) range of axis:",u[0],u[-1]))
        if verbose: logging.debug(strm("(_find_index) alias number is",alias_number))
        if verbose: logging.debug(strm("(_find_index) set origin from",origin, end=' '))
        origin -= alias_number * SW
        if verbose: logging.debug(strm("to",origin))
    p2 = np.argmin(abs(u-origin))
    assert np.count_nonzero(u[p2] == u) == 1, ("there seem to be"
            " "+repr(np.count_nonzero(u[p2] == u))+" values equal"
            " to "+repr(u[p2])+" but there should be only one")
    if abs(u[p2] - origin) > tolerance * max(abs(u[p2]),abs(origin)):
        p2_discrepancy = origin - u[p2]
    else:
        p2_discrepancy = None
    if verbose: logging.debug(strm("(_find_index) for origin",origin,"I am returning p2",p2,"out of",N,"and discrepancy",p2_discrepancy))
    alias_number += 1 # because the way that _ft_shift works essentially entails one aliasing
    return p2,p2_discrepancy,alias_number*SW
