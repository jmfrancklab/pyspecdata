from pyspecdata import *
fl = figlist_var()
a = nddata(r_[-1:1:100j],'x').set_units('x','m')
a.set_error(0.1)
b = a.C
b.run(lambda t: t**2+0.1)
b.set_error(0.1)
print(b.get_error('x'))
fl.next('$a$', figsize=(4.5,2))
fl.plot(a)
fl.next('$b$', figsize=(4.5,2))
fl.plot(b)
fl.next('aob', figsize=(4.5,2))
fl.plot(a/b)
title('$a/b$')
fl.show()

