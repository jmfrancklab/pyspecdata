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
y.add_noise(2)
# }}}
plot(y,'o')
c = y.polyfit('x', order=2)
assert len(c)==3
# math for min:
# a(x-b)²= ax² - 2abx + ab²
# c₂ = a
# c₁ = -2ab
# c₀ = ab²
# b = -c₁/(c₂2)
print("I found the minimum here at",-c[1]/c[2]/2)
# generate a polynomial that's more finely spaced by setting the
# `npts` parameter.  This is a shortcut for:
# x = nddata(r_[0:10:300j], 'x')
# followed by
# plot(x.eval_poly(c,'x'))
plot(y.eval_poly(c,'x', npts=300))
# {{{ not a good idea, but force the y intercept to 0
#     to show the code works
c = y.polyfit('x', order=3, force_y_intercept=0)
x = nddata(r_[0:10:100j], 'x')
plot(y.eval_poly(c,'x', npts=300))
# }}}
show()
