from pylab import *
from pyspecdata import *

# Generate fake data
t_axis = nddata(r_[0:2:2048j],'t2')
s = exp(1j*2*pi*5*t_axis - t_axis/800e-3 )
s += exp(1j*2*pi*-30*t_axis - t_axis/800e-3)
ph1 = nddata(r_[0:4]/4.,'ph1')
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
# may need work to become fully 'algorithmic'
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

division_scale = 0.03
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

print(axes_list[0][0])

print(A.dimlabels[-1*A.ndim])
print(A.dimlabels[-1*A.ndim + 1])

# try integrating code from alex_subplots

yMajorLocator = lambda: mticker.MaxNLocator(steps=[1,2,5,10])
majorLocator = lambda: mticker.MaxNLocator(min_n_ticks=4, steps=[1,2,5,10])
minorLocator = lambda: mticker.AutoMinorLocator(n=5)
for outer_index in range(A.shape[-1*A.ndim]):
    for inner_index in range(A.shape[-1*A.ndim + 1]):
        print(outer_index,inner_index)



    #thisaxes.set_ylabel('repeats')
    #thisaxes.yaxis.set_minor_locator(minorLocator())
    #thisaxes.yaxis.set_ticks_position('both')
    #if j == len(list_of_axes)-1:
    #    thisaxes.xaxis.set_major_locator(majorLocator())
    #    #for the minor ticks, use no labels; default NullFormatter
    #    thisaxes.xaxis.set_minor_locator(minorLocator())
    #    thisaxes.xaxis.tick_top()
    #    labels = [item.get_text() for item in thisaxes.get_xticklabels()]
    #    empty_string_labels = ['']*len(labels)
    #    thisaxes.set_xticklabels(empty_string_labels)
    #elif j == 0:
    #    thisaxes.set_xlabel('this is the x axis')
    #    thisaxes.xaxis.set_major_locator(majorLocator())
    #    thisaxes.xaxis.set_minor_locator(minorLocator())
    #else:
    #    thisaxes.set_xticks([])
for outer_index in range(A.shape[-1*A.ndim]):
    for inner_index in range(A.shape[-1*A.ndim + 1]):
        temp = list(axes_list[outer_index,inner_index])
        this_ax = axes(temp)
        figure(1);
        image(s['ph1',outer_index]['ph2',inner_index],ax=axes(temp))
        axes(temp).set_ylabel(A.dimlabels[-2])
        axes(temp).yaxis.set_minor_locator(minorLocator())
        axes(temp).yaxis.set_ticks_position('both')
        axes(temp).set_xlabel(None)
        axes(temp).xaxis.set_ticks([])
        if (outer_index == 0) and (inner_index == 0):
            axes(temp).set_xlabel(A.dimlabels[-1])
            axes(temp).xaxis.set_major_locator(majorLocator())
            axes(temp).xaxis.set_minor_locator(minorLocator())
        

        
show();quit()
