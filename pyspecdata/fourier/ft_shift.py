"shift-related helper functions"
from numpy import zeros,r_,nonzero,isclose,empty_like
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
    print targetslice,sourceslice
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
def _find_zero_index(t):
    "identify the index where zero lives"
    p2 = nonzero(isclose(t,0.))[0]
    try:
        assert len(p2) == 1
    except:
        print "p2 gives",p2
        print "t around the first value is",t[p2[0]-3:p2[0]+3]
        raise ValueError("The axis should have exactly one value that is equal to zero -- otherwise we're going to need to add code to (1) incorporate a carrier frequency or (2) add a phase shift")
    p2 = p2[0]
    return p2
