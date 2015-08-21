"shift-related helper functions"
from numpy import zeros,r_,nonzero,isclose,empty_like,argmin
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
def ft_clear_startpoints(self,axis):
    ("clears memory of where the origins in the time and frequency domain are"
            " this is useful, e.g. when you want to ift and center about time=0"
            " by setting shift=True")
    startf_dict = self.get_prop("FT_start_freq")
    if startf_dict is not None and axis in startf_dict.keys():
        startf_dict.pop(axis)
        if len(startf_dict) == 0:
            self.unset_prop('FT_start_freq')
        else:
            self.set_prop('FT_start_freq',startf_dict)
    startt_dict = self.get_prop("FT_start_time")
    if startt_dict is not None and axis in startt_dict.keys():
        startt_dict.pop(axis)
        if len(startt_dict) == 0:
            self.unset_prop('FT_start_time')
        else:
            self.set_prop('FT_start_time',startt_dict)
    return self
def _find_zero_index(t):
    "identify the index where zero lives"
    assert t[0] <= 0, ("You seem to be trying to FT a sequence whose axis"
                        " starts at a value greater than zero."+thinkaboutit_message)
    p2 = nonzero(isclose(0.,t,atol = 0))[0]
    try:
        assert len(p2) == 1, ("The axis should have exactly one value that is equal to zero -- otherwise we're going to need to add code to (1) incorporate a carrier frequency or (2) add a phase shift")
    except:
        print "------------extra info for error dump----------------"
        print "p2 gives",p2
        if len(p2) == 0:
            print "t is",t
            print "t values that are closest to zero are",
            temp = argmin(abs(t-0))
            print t[temp-2:temp+3]
        else:
            print "t is",t
            print "t around the first value is",t[p2[0]-3:p2[0]+3]
        raise
    p2 = p2[0]
    return p2
