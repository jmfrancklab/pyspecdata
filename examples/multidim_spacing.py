from pylab import *
import numpy as np
from pyspecdata import *

# a_shape can be created from
# a_shape = ndshape(A)
# a_shape.pop(a_shape.dimlabels(-1))
# a_shape.pop(a_shape.dimlabels(-1))
grid_bottom = 0.2
grid_top = 0.8
total_spacing = 0.2
a_shape = ndshape([('ph1',2),('ph2',2),('ph3',4)])
divisions = []
for j,thisdim in enumerate(a_shape.dimlabels[::-1]):
    old = [j/2.0 for j in divisions]
    divisions = (old + [1])*(a_shape[thisdim]-1)+old
    print("for",thisdim,"I get",divisions)
divisions = [j*total_spacing/sum(divisions) for j in divisions]
axes_height = (grid_top-grid_bottom-total_spacing)/prod(a_shape.shape)
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
idx = nddata(r_[0:prod(a_shape.shape)],[-1],['smooshed'])
idx.chunk('smooshed',a_shape.dimlabels,a_shape.shape)
print(idx)
show()
