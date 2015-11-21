from ..general_functions import *
from pylab import * 

def skew(self,along_axis,by_amount,propto_axis,zero_fill = True):
    kwargs = {}
    if zero_fill:
        n = self.data.shape[self.axn(along_axis)]
        n = int(2**(ceil(log2(n))))
        n *= 2
        {'pad':n}
    self.ft(along_axis,shift = None,**kwargs)
    newdata = self * self.fromaxis([along_axis,propto_axis],
            lambda x,y: exp(-2j*pi*by_amount*x*y))
    self.data = newdata.data
    self.ift(along_axis)
    return self
