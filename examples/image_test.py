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
yMajorLocator = lambda: mticker.MaxNLocator(steps=[1,2,5,10])
majorLocator = lambda: mticker.MaxNLocator(min_n_ticks=4, steps=[1,2,5,10])
minorLocator = lambda: mticker.AutoMinorLocator(n=5)
for j,b in enumerate(axes_bottom):
    ax_list.append(axes([0.2,b,0.7,axes_height])) # lbwh
    if j == 0:
        #ax_list[-1].set_xlabel(a_shape.dimlabels[-1])
        ax_list[-1].xaxis.set_major_locator(majorLocator())
        ax_list[-1].xaxis.set_minor_locator(minorLocator())
    elif (j == len(axes_bottom)-1):
        ax_list[-1].xaxis.set_major_locator(majorLocator())
        ax_list[-1].set_xlabel(None)
        #for the minor ticks, use no labels; default NullFormatter
        ax_list[-1].xaxis.set_minor_locator(minorLocator())
        ax_list[-1].xaxis.tick_top()
        labels = [item.get_text() for item in ax_list[-1].get_xticklabels()]
        empty_string_labels = ['']*len(labels)
        ax_list[-1].set_xticklabels(empty_string_labels)
        ax_list[-1].set_xlabel(None)
    else:
        ax_list[-1].xaxis.set_ticks([])
        ax_list[-1].get_xaxis().set_visible(False)
        ax_list[-1].set_xlabel(None)
    ax_list[-1].set_ylabel(a_shape.dimlabels[-2])
    ax_list[-1].yaxis.set_minor_locator(minorLocator())
    ax_list[-1].yaxis.set_ticks_position('both')

A = s.smoosh(a_shape.dimlabels[:-2],'smooshed',noaxis=True)
A.reorder('smooshed',first=True)
for j in range(len(ax_list)):
    image(A['smooshed',j],ax=ax_list[j])
    if not j == 0:
        ax_list[j].set_xlabel(None)


# to drop into ax_list, just do
# A.smoosh(a_shape.dimlabels, 'smooshed', noaxis=True)
# in ax_list[0] put A['smooshed',0], etc
idx = nddata(r_[0:prod(a_shape.shape[:-2])],[-1],['smooshed'])
idx.chunk('smooshed',a_shape.dimlabels[:-2],a_shape.shape[:-2])

def draw_span(ax1, ax2, label, allow_for_text=10, allow_for_ticks=40):
    x1,y1 = ax1.transAxes.transform(r_[0,1])
    x2,y2 = ax2.transAxes.transform(r_[0,0])
    x1-=allow_for_ticks
    x_text = x1-allow_for_ticks
    x2-=allow_for_ticks
    x1,y1 = fig.transFigure.inverted().transform(r_[x1,y1])
    x_text,_ = fig.transFigure.inverted().transform(r_[x_text,0])
    x2,y2 = fig.transFigure.inverted().transform(r_[x2,y2])
    lineA = lines.Line2D([x1,x2],[y1,y2],
            linewidth=3, color='r', transform=fig.transFigure,
            clip_on=False)
    text(x_text, (y2+y1)/2, label, va='center', ha='right', rotation=90, transform=fig.transFigure, color='r')
    fig.add_artist(lineA)

for thisdim in a_shape.dimlabels[:-2]:
    # generate labels for the dimensions, outside in
    # use definition of idx in code
    for j in range(a_shape[thisdim]):
        first_axes = ax_list[idx[thisdim,j].data.ravel()[0]]
        last_axes = ax_list[idx[thisdim,j].data.ravel()[-1]]
        print(first_axes)
        print(last_axes)
        draw_span(first_axes,last_axes,"%s=%d"%(thisdim,j))
show();quit()
