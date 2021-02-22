from pylab import *
from pyspecdata import *

# Generate fake data
t_axis = nddata(r_[0:2:2048j],'t2')
s = exp(1j*2*pi*5*t_axis - t_axis/800e-3 )
s += exp(1j*2*pi*-30*t_axis - t_axis/800e-3)
ph1 = nddata(r_[0:8]/4.,'ph1')
ph2 = nddata(r_[0,2]/4.,'ph2')
# this cannot start at 0 since we multiply s by it
repeats = nddata(r_[1:6],'repeats')
s *= repeats/repeats
s.add_noise(0.3)
s *= exp(1j*2*pi*ph1)
s *= exp(1j*2*pi*ph2)
#s['t2',0] *= 0.5
#s.ft('t2',shift=True)
#s.ft(['ph1','ph2'])
s.reorder(['repeats','t2'],first=False)
print(ndshape(s))

A = ndshape(s)
A.ndim = len(shape(A))
print(A)
print(A.ndim)

# begin new code
# determine list of divisions
# still needs to be put into algorithm for more than 3 dim
div_list = [1]
#div_list = (div_list + [2])*(A.shape[-1*A.ndim] - 2) + div_list
div_list = (div_list )*(A.shape[-1*A.ndim + 1] - 2) + div_list
div_list = (div_list + [2])*(A.shape[-1*A.ndim] - 1) + div_list

# determine num of axes objects needed
num_axes_obj = 1
for dim_idx in range(A.ndim):
    counter = -1*(A.ndim - dim_idx)
    if counter == -1:
        break
    elif counter == -2:
        break
    else:
        num_axes_obj *= A.shape[counter]
print(num_axes_obj)
axes_list = [[1.] * 4] * num_axes_obj
#axes_list = zeros(num_axes_obj,4)
print(axes_list)
print(len(axes_list))
axes_list = np.array(axes_list)
print(axes_list.shape)
reshape_tuple = list(A.shape[-1*A.ndim+x] for x in range(0,A.ndim-2))
reshape_tuple.append(-1)
reshape_tuple = tuple(reshape_tuple)
axes_list = axes_list.reshape(reshape_tuple)

top_border = 0.1
bottom_border = 0.1
left_border = 0.1
right_border = 0.1

#div_list.insert(0,0)
division_scale = 0.01
division_space = sum(div_list) * division_scale

height = (1. - (top_border+bottom_border+division_space))/float(num_axes_obj)
width = 1. - (left_border+right_border)

div_counter = 0
print(div_list)


for outer_index in range(A.shape[-1*A.ndim]):
    for inner_index in range(A.shape[-1*A.ndim + 1]):
        if div_counter == 0:
            axes_list[outer_index,inner_index] = [ left_border,
                    bottom_border,
                    width,
                    height ]
            div_counter += 1
        else:
            axes_list[outer_index,inner_index] = [ left_border,
                    bottom_border + div_counter*height + sum(div_list[:div_counter])*division_scale,
                    width,
                    height ]
            div_counter += 1
# end newer code

for outer_index in range(A.shape[-1*A.ndim]):
    for inner_index in range(A.shape[-1*A.ndim + 1]):
        temp = list(axes_list[outer_index,inner_index])
        figure(1);
        image(s['ph1',outer_index]['ph2',inner_index],ax=axes(temp))
show();quit()
