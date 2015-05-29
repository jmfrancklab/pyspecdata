from fornotebook import *
import sys
import interptau
interpresult,tau,rho = interptau.interptau(sys.argv[1],sys.argv[2])
print "interp: %0.1fps"%(interpresult*1e12)
if len(sys.argv) > 3:
    print "retardation: %0.1f"%(interpresult*1e12/double(sys.argv[3]))
