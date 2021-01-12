"""
Error Propagation
=================

Demonstrating automatic error propagation for a set of simple functions.
"""
from pylab import *
from pyspecdata import *
fl = figlist_var()
a = nddata(r_[-1:1:100j],'x').set_units('x','m')
a.set_error(0.1)
b = a.C
b.run(lambda t: t**2+0.1)
b.set_error(0.1)
print(b.get_error('x'))
fl.next('$a$', figsize=(9,4))
fl.plot(a)
fl.next('$b$', figsize=(9,4))
fl.plot(b)
fl.next('aob', figsize=(9,4))
fl.plot(a/b)
title('$a/b$')
fl.show()

