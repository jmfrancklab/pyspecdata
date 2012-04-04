from matlablike import *
import sys
import interptau
if len(sys.argv) < 2:
    raise CustomError('call as python dettau.py xi nmrfreq/MHz')
interpresult,tau,rho = interptau.interptau(sys.argv[1],sys.argv[2])
print "interp:",interpresult
semilogx(tau,rho)
semilogx(r_[interpresult],r_[double(sys.argv[1])],'o')
axis('tight')
show()
