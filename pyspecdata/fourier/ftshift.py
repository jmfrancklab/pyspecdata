def ftshift(self,axis):
    self.data = fftshift(self.data,axes = self.axn(axis))
    x = self.getaxis(axis)
    x[:] = fftshift(x)
    j = len(x)/2 # given the floor, this works out to be the central index
    x_subset = x[:j]
    x_subset -= x_subset[-1] + x[j+1] # zero and set to this
    return self
