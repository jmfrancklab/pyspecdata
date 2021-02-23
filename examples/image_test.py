from pylab import *
from pyspecdata import *
import matplotlib.lines as lines

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

grid_bottom = 0.2
grid_top = 0.8
total_spacing = 0.2
a_shape = ndshape(s)
divisions = []
for j,thisdim in enumerate(a_shape.dimlabels[:-2]):
    old = [j/2.0 for j in divisions]
    divisions = (old + [1])*(a_shape[thisdim]-1)+old
    print("for",thisdim,"I get",divisions)
divisions = [j*total_spacing/sum(divisions) for j in divisions]
axes_height = (grid_top-grid_bottom-total_spacing)/prod(a_shape.shape[:-2])
axes_bottom = np.cumsum([axes_height+j for j in divisions]) # becomes ndarray
axes_bottom = r_[0,axes_bottom]
axes_bottom += grid_bottom
axes_top = grid_bottom + grid_top
fig = figure()
ax_list = []
for j,b in enumerate(axes_bottom):
    ax_list.append(axes([0.2,b,0.7,axes_height], figure=fig)) # lbwh
    ax_list[-1].get_xaxis().set_visible(False)
# to drop into ax_list, just do
# A.smoosh(a_shape.dimlabels, 'smooshed', noaxis=True)
# in ax_list[0] put A['smooshed',0], etc
idx = nddata(r_[0:prod(a_shape.shape[:-2])],[-1],['smooshed'])
idx.chunk('smooshed',a_shape.dimlabels[:-2],a_shape.shape[:-2])

show();quit()
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
left_border = 0.2
right_border = 0.2

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
for outer_index in range(A.shape[-1*A.ndim]):
    for inner_index in range(A.shape[-1*A.ndim + 1]):
        temp = list(axes_list[outer_index,inner_index])
        this_ax = axes(temp)
        fig = figure(1)
        fig;
        image(s['ph1',outer_index]['ph2',inner_index],ax=axes(temp))
        axes(temp).set_ylabel(A.dimlabels[-2])
        axes(temp).yaxis.set_minor_locator(minorLocator())
        axes(temp).yaxis.set_ticks_position('both')
        axes(temp).set_xlabel(None)
        axes(temp).xaxis.set_ticks([])

        ## Put inner phase cycle labels, according to dimension and value
        # A little tricker than the outer phase cycle labels
        x1,y1 = axes(temp).transAxes.transform(r_[0,1])
        x2,y2 = axes(temp).transAxes.transform(r_[0,0])
        x1-=40
        x2-=40
        x_text = x1-50
        x1-=49
        x2-=49
        x1,y1 = fig.transFigure.inverted().transform(r_[x1,y1])
        x_text,_ = fig.transFigure.inverted().transform(r_[x_text,0])
        x2,y2 = fig.transFigure.inverted().transform(r_[x2,y2])
        axes_obj_label = A.dimlabels[-1*A.ndim+1]+'=%d'%(inner_index)
        text(x_text, temp[1]+0.04, axes_obj_label, va='center', ha='right', rotation=90,
                transform = fig.transFigure, color='k')
        line_inner = lines.Line2D([x1,x2],[y1,y2],linewidth=3, color='k',
                transform=fig.transFigure,clip_on=False)
        fig.add_artist(line_inner)

        # Put outer phase cycle labels, according to dimension and value
        # Set these each time you reach the end of the inner loop
        if (inner_index == A.shape[-1*A.ndim+1]-1):
            x1,y1 = axes(temp).transAxes.transform(r_[0,1])
            x2,y2 = axes(temp).transAxes.transform(r_[0,0])
            x1-=40
            x2-=40
            x_text = x1-100
            x1-=90
            x2-=90
            y1-=25
            y2-=25
            x1,y1 = fig.transFigure.inverted().transform(r_[x1,y1])
            x_text,_ = fig.transFigure.inverted().transform(r_[x_text,0])
            x2,y2 = fig.transFigure.inverted().transform(r_[x2,y2])
            axes_obj_label = A.dimlabels[-1*A.ndim]+'=%d'%(outer_index)
            text(x_text, temp[1], axes_obj_label, va='center', ha='right', rotation=90,
                    transform = fig.transFigure, color='k')
            line_outer = lines.Line2D([x1,x2],[y1,y2],linewidth=3, color='k',
                    transform=fig.transFigure,clip_on=False)
            fig.add_artist(line_outer)


        # Put x-axis labels and ticks on bottom-most axes object
        if (outer_index == 0) and (inner_index == 0):
            axes(temp).set_xlabel(A.dimlabels[-1])
            axes(temp).xaxis.set_major_locator(majorLocator())
            axes(temp).xaxis.set_minor_locator(minorLocator())
        # Put x-axis ticks on top-most axes object
        if (outer_index == A.shape[-1*A.ndim]-1) and (inner_index == A.shape[-1*A.ndim+1]-1):
            axes(temp).xaxis.set_major_locator(majorLocator())
            #for the minor ticks, use no labels; default NullFormatter
            axes(temp).xaxis.set_minor_locator(minorLocator())
            axes(temp).xaxis.tick_top()
            labels = [item.get_text() for item in axes(temp).get_xticklabels()]
            empty_string_labels = ['']*len(labels)
            axes(temp).set_xticklabels(empty_string_labels)


show();quit()
