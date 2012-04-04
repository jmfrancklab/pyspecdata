from fornotebook import *
import sys
import interptau
interpresult,tau,rho = interptau.interptau(sys.argv[1],sys.argv[2])
print "interp: %0.1fps",interpresult*1e12
