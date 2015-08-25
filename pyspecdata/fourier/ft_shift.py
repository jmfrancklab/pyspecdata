"shift-related helper functions"
from numpy import zeros,r_,nonzero,isclose,empty_like,argmin,count_nonzero
thinkaboutit_message = ("If you think about it, you"
                        " probably don't want to do this.  You either want to fill with"
                        " zeros from zero up to the start or you want to first set the"
                        " start point to zero.")
def _ft_shift(self,thisaxis,p2,shift_axis = None):
    ("perform a generalized fftshift along the axis indicated by the integer `thisaxis`, where `p2` gives the index that will become the first index"
    "\n this is derived from the numpy fftshift routine, but defines slices instead of index numbers"
    "\n `shift_axis` is only used after the (i)fft.  It assumes that the axis labels start at zero, and it aliases them over in the same way the data was aliased")
    if p2 == 0: # nothing to be done
        return self
    newdata = empty_like(self.data)
    n = self.data.shape[thisaxis]
    sourceslice = [slice(None,None,None)] * len(self.data.shape)
    targetslice = [slice(None,None,None)] * len(self.data.shape)
    # move second half first -- the following are analogous to the numpy function, but uses slices instead
    sourceslice[thisaxis] = slice(p2,n)
    targetslice[thisaxis] = slice(None,n-p2)
    newdata[targetslice]  = self.data[sourceslice]
    # move first half second (the negative frequencies)
    sourceslice[thisaxis] = slice(None,p2)
    targetslice[thisaxis] = slice(-p2,None)
    newdata[targetslice]  = self.data[sourceslice]
    self.data = newdata
    if shift_axis is not None and shift_axis:
        axisname = self.dimlabels[thisaxis]
        x = self.getaxis(axisname)
        newaxis = empty_like(x)
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
def ft_clear_startpoints(self,axis,t=None,f=None):
    ("clears memory of where the origins in the time and frequency domain are"
            " this is useful, e.g. when you want to ift and center about time=0"
            " by setting shift=True you can also manually set the points:"
            "\n\tkeyword arguments `t` and `f` can be set by (1) manually"
            " setting the start point (2) using the string 'current' to leave the"
            " current setting alone or (3) None, which clears the startpoint")
    if f is 'current':
        pass
    else:
        startf_dict = self.get_prop("FT_start_freq")
        if f is None: # then unset
            if startf_dict is not None and axis in startf_dict.keys():
                startf_dict.pop(axis)
                if len(startf_dict) == 0:
                    self.unset_prop('FT_start_freq')
        elif startf_dict is None:
            self.set_prop("FT_start_freq",{axis:f})
        else:
            startf_dict[axis] = f
    if t is 'current':
        pass
    else:
        startt_dict = self.get_prop("FT_start_time")
        if t is None: # then unset
            if startt_dict is not None and axis in startt_dict.keys():
                startt_dict.pop(axis)
                if len(startt_dict) == 0:
                    self.unset_prop('FT_start_time')
        elif startt_dict is None:
            self.set_prop("FT_start_time",{axis:t})
        else:
            startt_dict[axis] = t
    return self
def _find_index(u,origin = 0.0,tolerance = 1e-5):
    ("identify the index of `u` (represents either time or frequency) where"
            " `origin` lives -- if it finds a value exactly equal"
            " to `origin`, returns `(p2,None)` -- otherwise, `None` is replaced by"
            " u-shift indicating how far the position"
            " marked by `p2` is ahead of `origin`")
    assert origin >= 0, ("The origin must be >=0.  If you want to find a"
            " negative frequency that's aliased over, you need to add the"
            " spectral width yourself")
    p2 = argmin(abs(u-origin))
    assert count_nonzero(u[p2] == u) == 1, ("there seem to be"
            " "+repr(count_nonzero(u[p2] == u))+" values equal"
            " to "+repr(u[p2])+" but there should be only one")
    du = u[1]-u[0]
    if abs(u[p2] - origin) > tolerance * max(abs(u[p2]),abs(origin)):
        p2_discrepancy = u[p2] - origin # marks where the p2 position really is vs. where we want it to be
    else:
        p2_discrepancy = None
    print "for origin",origin,"I am returning p2",p2,"and discrepancy",p2_discrepancy
    return p2,p2_discrepancy
