from pyspecdata.fornotebook import *
import sys
import interptau
if len(sys.argv) > 4 or len(sys.argv) < 3:
    print "example:"
    print "\tpython printtau.py 0.01 14.8 54"
    print "\t0.01 -- xi"
    print "\t14.8 -- NMR freq / MHz (for field)"
    print "\t54   -- tau / ps (OPTIONAL, reference for bulk water)"
    exit()
Dw_default = 2.3e-9 # Brandon's
Dsl_default = 4.1e-10 # Brandon's
tau_reference = 54
if len(sys.argv) > 3:
    tau_reference = double(sys.argv[3])
tau_reference *= 1e-12
xi = double(sys.argv[1])
nmr_freq = double(sys.argv[2])
interpresult,tau,rho = interptau.interptau(xi,nmr_freq)
print "xi input:",xi
print "interp: %0.1fps"%(interpresult*1e12)
retardation = interpresult / tau_reference
print "retardation: %0.1f"%(retardation)
print "diffusion: %0.3g m^2/s x 1e-11"%(1/retardation*(Dw_default)/1e-11)
if len(sys.argv) > 3:
    print "retardation: %0.1f"%(interpresult*1e12/double(sys.argv[3]))
