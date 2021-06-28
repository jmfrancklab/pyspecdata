"""
Polynomial Fitting
==================

A quick demonstration of polynomial fitting.
"""
from pylab import *
from pyspecdata import *
# {{{ generate fake data
x = nddata(r_[0:10:7j], 'x')
y = (x-2.0)**2
y.add_noise(5)
# }}}
plot(y,'o')
c = y.polyfit('x', order=2)
assert len(c)==3
# generate a polynomial that's more finely spaced
x = nddata(r_[0:10:100j], 'x')
plot(x.eval_poly(c,'x'))
# {{{ not a good idea, but force the y intercept to 0
#     to show the code works
c = y.polyfit('x', order=3, force_y_intercept=0)
x = nddata(r_[0:10:100j], 'x')
plot(x.eval_poly(c,'x'))
# }}}
show()
