def _ft_shift(self,thisaxis,p2):
    ("perform a generalized fftshift along the axis indicated by the integer `thisaxis`, where `p2` gives the index that will become the first index"
    "\n this is derived from the numpy fftshift routine, but defines slices instead of index numbers")
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
    return self
