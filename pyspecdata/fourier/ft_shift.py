"shift-related helper functions"
def _ft_shift(self,thisaxis,p2,shift_axis = True):
    ("perform a generalized fftshift along the axis indicated by the integer `thisaxis`, where `p2` gives the index that will become the first index"
    "\n this is derived from the numpy fftshift routine, but defines slices instead of index numbers"
    "\n `shift_axis` is an option because it is wasteful to shift the axis before the fft")
    newdata = list(shape(self.data))
    newdata[thisaxis] = padded_length
    newdata = zeros(tuple(newdata),dtype = self.data.dtype)
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
    if shift_axis:
        axisname = self.dimlabels[thisaxis]
        x = self.getaxis(axisname)
        n = len(x)
        # move second half first -- the following are analogous to the numpy function, but uses slices instead
        sourceslice = slice(p2,n)
        targetslice = slice(None,n-p2)
        newaxis[targetslice]  = x[sourceslice]
        # move first half second (the negative frequencies)
        sourceslice[thisaxis] = slice(None,p2)
        targetslice[thisaxis] = slice(-p2,None)
        newaxis[targetslice]  = x[sourceslice]
        self.setaxis(axisname,newaxis)
    return self
def _find_zero_index(t):
    "identify the index where zero lives"
    p2 = nonzero(t==0.)[0]
    try:
        assert len(p2) == 1
    except:
        raise ValueError("The axis should have exactly one value that is equal to zero -- otherwise we're going to need to add code to (1) incorporate a carrier frequency or (2) add a phase shift")
    p2 = p2[0]
    return p2
