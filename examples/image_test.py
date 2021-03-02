from pylab import *
from pyspecdata import *
import matplotlib.lines as lines

# Generate fake data
t_axis = nddata(r_[0:2:2048j],'t2')
s = exp(1j*2*pi*5*t_axis - t_axis/800e-3 )
s += exp(1j*2*pi*-30*t_axis - t_axis/800e-3)
ph1 = nddata(r_[0:4]/4.,'ph1')
ph2 = nddata(r_[0,2]/4.,'ph2')
ph3 = nddata(r_[0,1,2]/4.,'ph3')
# this cannot start at 0 since we multiply s by it
repeats = nddata(r_[1:6],'repeats')
s *= repeats/repeats
s.add_noise(0.3)
s *= exp(1j*2*pi*ph1)
s *= exp(1j*2*pi*ph2)
s *= exp(1j*2*pi*ph3)
#s['t2',0] *= 0.5
#s.ft('t2',shift=True)
#s.ft(['ph1','ph2'])
s.reorder(['repeats','t2'],first=False)
s.reorder('ph2',first=True)
print(ndshape(s))

fl = figlist_var()

def image_new(this_nddata,this_fig_obj):
    grid_bottom = 0.0
    bottom_pad = 0.15
    grid_bottom += bottom_pad
    grid_top = 1.0
    top_pad = 0.05
    grid_top -= top_pad
    total_spacing = 0.2
    a_shape = ndshape(s)
    num_dims = len(a_shape.dimlabels[:-2])
    divisions = []
    # should be looping in backward order from printed shape
    for j,thisdim in enumerate(a_shape.dimlabels[::-1][2:]):
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

    #labels_space = num_dims*
    LHS_pad = 0.05
    RHS_pad = 0.05
    LHS_labels = 0.08*num_dims
    width = 1.-(LHS_pad+RHS_pad+LHS_labels)

    for j,b in enumerate(axes_bottom):
        ax_list.append(axes([LHS_labels+LHS_pad,b,width,axes_height])) # lbwh
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
        #ax_list[-1].set_ylabel(a_shape.dimlabels[-2])
        #ax_list[-1].yaxis.set_minor_locator(minorLocator())
        #ax_list[-1].yaxis.set_ticks_position('both')

    if len(a_shape.dimlabels) > 3:
        A = s.smoosh(a_shape.dimlabels[:-2],'smooshed',noaxis=True)
        A.reorder('smooshed',first=True)
    else:
        A = this_nddata.C
        A.rename(a_shape.dimlabels[:-2][0],'smooshed')

    def draw_span(ax1, ax2, label, this_label_num, allow_for_text=10, allow_for_ticks=100):
        x1,y1 = ax1.transAxes.transform(r_[0,1])
        x2,y2 = ax2.transAxes.transform(r_[0,0])
        x1-=allow_for_ticks
        x_text = x1-allow_for_text
        x2-=allow_for_ticks
        # following line to create an offset for different dimension labels
        label_spacing = this_label_num*70
        x1,y1 = fig.transFigure.inverted().transform(r_[x1-label_spacing,y1])
        x_text,_ = fig.transFigure.inverted().transform(r_[x_text-label_spacing,0])
        x2,y2 = fig.transFigure.inverted().transform(r_[x2-label_spacing,y2])
        lineA = lines.Line2D([x1,x2],[y1,y2],
                linewidth=3, color='k', transform=fig.transFigure,
                clip_on=False)
        text(x_text, (y2+y1)/2, label, va='center', ha='right', rotation=90, transform=fig.transFigure, color='k')
        fig.add_artist(lineA)

    label_placed = zeros(num_dims)

    def place_labels(ax1, label, label_placed, this_label_num, check_for_label_num = True,
            allow_for_text=10, allow_for_ticks=100, y_adjustment = 55):
        if check_for_label_num:
            if not label_placed[this_label_num]:
                x1,y1 = ax1.transAxes.transform(r_[0,1])
                x1-=allow_for_ticks
                x_text = x1-allow_for_text
                label_spacing = this_label_num*65
                y1 -= y_adjustment
                #x1,y1 = fig.transFigure.inverted().transform(r_[x1-label_spacing,y1])
                x_text,y1 = fig.transFigure.inverted().transform(r_[x_text-label_spacing,y1])
                text(x_text, y1, label, va='center', ha='right', rotation=45, transform=fig.transFigure, color='k')
                label_placed[this_label_num] = 1
        else:
            x1,y1 = ax1.transAxes.transform(r_[0,2])
            x1-=allow_for_ticks
            x_text = x1-allow_for_text
            label_spacing = this_label_num*65
            y1 -= y_adjustment
            #x1,y1 = fig.transFigure.inverted().transform(r_[x1-label_spacing,y1])
            x_text,y1 = fig.transFigure.inverted().transform(r_[x_text-label_spacing,y1])
            text(x_text, y1, label, va='center', ha='right', rotation=45, transform=fig.transFigure, color='k')
            labels = [item.get_text() for item in ax1.get_xticklabels()]
            empty_string_labels = ['']*len(labels)
            ax1.set_xticklabels(empty_string_labels)
            ax1.set_xlabel(None)
            ax1.tick_params(bottom=False)
            ax1.set_ylabel(None)
            ax1.set_yticklabels(empty_string_labels)

    for j in range(len(ax_list)):
        image(A['smooshed',j],ax=ax_list[j])
        ax_list[j].set_ylabel(None)
        if not j == 0:
            ax_list[j].set_xlabel(None)

    # to drop into ax_list, just do
    # A.smoosh(a_shape.dimlabels, 'smooshed', noaxis=True)
    # in ax_list[0] put A['smooshed',0], etc
    idx = nddata(r_[0:prod(a_shape.shape[:-2])],[-1],['smooshed'])
    idx.chunk('smooshed',a_shape.dimlabels[:-2],a_shape.shape[:-2])
    remaining_dim = a_shape.dimlabels[:-2]
    depth = num_dims
    def decorate_axes(idx,remaining_dim,depth):
        thisdim=remaining_dim[0]
        print("This is remaining dim",remaining_dim)
        print("This dim is",thisdim)
        print(ndshape(idx))
        depth -= 1
        for j in range(a_shape[thisdim]):
            idx_slice = idx[thisdim,j]
            print("For",thisdim,"element",j,idx_slice.data.ravel())
            first_axes = ax_list[idx_slice.data.ravel()[0]]
            last_axes = ax_list[idx_slice.data.ravel()[-1]]
            draw_span(last_axes,first_axes,"%d"%(j),
                    this_label_num=depth)
            place_labels(ax_list[0],"%s"%(thisdim), label_placed,
                    this_label_num=depth)
            new_remaining_dim = remaining_dim[1:]
            if len(remaining_dim) > 1:
                decorate_axes(idx_slice,new_remaining_dim,depth)
    print("call recursive function")
    decorate_axes(idx,remaining_dim,depth)
    place_labels(axes([LHS_labels+LHS_pad,axes_bottom[0],width,0]),"%s"%(a_shape.dimlabels[-2]), label_placed,
            this_label_num=0, check_for_label_num = False, allow_for_text = -75, y_adjustment=55)
    return 

image_new(s,fl)
fl.show()


show();quit()
