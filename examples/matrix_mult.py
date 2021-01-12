"""
Matrix Multiplication
=====================

Various ways of implementing different matrix multiplications.
Read the documentation embedded in the code.
"""
# -*- coding: utf-8 -*-
from pylab import *
from pyspecdata import *
from numpy.random import random
import time
init_logging('debug')

# %%
# In this example, the assertions essentially tell the story of what's going on

# %%
# Note that in all these examples, the pyspecdata version *appears* more
# complicated.
# But, that's because these are toy examples, where we have no need for the
# dimension names or axes.
# Nonetheless, we wanted to give the simplest working example possible.

# %%
# First, we demonstrate matrix multiplication
# for all the below, I attach an axis to make sure the routines work with the
# axes attached

a_nd = nddata(random(10*2048),[10,2048],['x','y']).setaxis('x','#').setaxis('y','#')
a = a_nd.data
# %%
# in the next line, note how only the dimension that goes away is named the
# same!
#
# if you think about the matrix a transforming from one vector space (labeled
# y) to another (labeled x) this makes sense
a2_nd = nddata(random(10*2048),[2048,10],['y','z']).setaxis('y','#').setaxis('z','#')
a2 = a2_nd.data

# multiply two different matrices

time1 = time.time()
b = a @ a2
time2 = time.time()
b_nd = a_nd @ a2_nd
# %%
# the previous is unambiguous b/c only 'y' is shared between the two,
# but I can do the following for clarity:
# b_nd = a_nd.along('y') @ a2_nd
# 
# Note that "along" gives the dimension along which the sum is performed -- and
# so this dimension goes away upon matrix multiplication.
# If only one dimension is shared between the matrices, then we know to take
# the sum along the shared dimension.
# For example, here a2_nd transforms from a space called "z" into a space called "y",
# while a_nd transforms from "y" into "x" -- so it's obvious that a_nd @ a2_nd should
# transform from "z" into "y".
time3 = time.time()
assert b_nd.dimlabels == ['x','z'], b_nd.dimlabels
assert all(isclose(b,b_nd.data))
print("total time",(time3-time2),"time/(time for raw)",((time3-time2)/(time2-time1)))
assert ((time3-time2)/(time2-time1))<1

# %%
# calculate a projection matrix

time1 = time.time()
b = a @ a.T
time2 = time.time()
# %%
# note that here, I have to rename the column space
b_nd = a_nd.along('y',('x','x_new')) @ a_nd
time3 = time.time()
assert b_nd.dimlabels == ['x_new','x'], b_nd.dimlabels
assert all(b_nd.getaxis('x_new') == b_nd.getaxis('x'))
assert (id(b_nd.getaxis('x_new')) != id(b_nd.getaxis('x')))
assert all(isclose(b,b_nd.data))
if time2-time1>0:
    print("total time",(time3-time2),"time/(time for raw)",((time3-time2)/(time2-time1)))
    assert ((time3-time2)/(time2-time1))<1.1

# %%
# now, a standard dot product note how I don't need `along` here, since it's
# unambiguous

a_nd = nddata(random(10),[10],['myaxis']).setaxis('myaxis','#')
b_nd = nddata(random(10),[10],['myaxis']).setaxis('myaxis','#')
a = a_nd.data
b = b_nd.data
assert all(isclose(a.dot(b),(a_nd @ b_nd).data))

# %%
# Finally, let's show what happens when we multiply a matrix by itself and
# *don't* rename one of the dimensions
#
# By doing this, we indicate that we're not interested in transforming from one
# vector space to another (as a projection matrix does), but rather just have
# two sets of vectors and are interested in finding the dot products between
# the two sets
# 
# This will take the dot product of our 10 2048-long vectors, and present them
# 10-long array

a_nd = nddata(random(10*2048),[10,2048],['x','y']).setaxis('x','#').setaxis('y','#')
a = a_nd.data
b_nd = a_nd.along('y') @ a_nd
b = matmul(a_nd.data.reshape(10,1,2048),
        a_nd.data.reshape(10,2048,1)).reshape(-1)
assert all(isclose(b,b_nd.data))
assert len(b.data) == 10
