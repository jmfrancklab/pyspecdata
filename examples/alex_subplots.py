from pyspecdata import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as lines
from matplotlib.path import Path
import matplotlib.patches as patches

bottoms = [0.13, 0.35, 0.63, 0.8]
verts = [
   (0., 0.),   # start
   (0.1, 0.),  # ctrl
   (0.1, 0.5),  # end
   (0.1, 1.0),  # ctrl
   (0.2, 1.0),  # end
   (0.25, 1.0),  # ctrl
   (0.25, 1.0),  # end
]

codes = [
    Path.MOVETO,# start
    Path.CURVE3, # ctrl
    Path.CURVE3, # end
    Path.CURVE3, # ctrl
    Path.CURVE3, # end
    Path.CURVE3, # ctrl
    Path.CURVE3, # end
]

height = 0.1
fig=plt.figure(1)
list_of_axes = []
depth_of_labels = 2 # here we have two nested dimensions
pixels_per_depth = 70 
l_pad,_ = fig.transFigure.inverted().transform(r_[pixels_per_depth*(1+depth_of_labels),0])# 30 pixels per depth
r_pad,_ = fig.transFigure.inverted().transform(r_[pixels_per_depth,0])# 30 pixels per depth
print("l_pad is",l_pad)
for thisbottom in bottoms:
    list_of_axes.append(plt.axes([l_pad,thisbottom,1-l_pad-r_pad,height]))
dimension_sizes = [2,2] # your loops should be based off of this
#}}}
#{{{adjusting tick marks and labels

yMajorLocator = lambda: mticker.MaxNLocator(steps=[1,2,5,10])
majorLocator = lambda: mticker.MaxNLocator(min_n_ticks=4, steps=[1,2,5,10])
minorLocator = lambda: mticker.AutoMinorLocator(n=5)
for j,thisaxes in enumerate(list_of_axes):
    thisaxes.set_ylabel('repeats')
    thisaxes.yaxis.set_minor_locator(minorLocator())
    thisaxes.yaxis.set_ticks_position('both')
    if j == len(list_of_axes)-1:
        thisaxes.xaxis.set_major_locator(majorLocator())
        #for the minor ticks, use no labels; default NullFormatter
        thisaxes.xaxis.set_minor_locator(minorLocator())
        thisaxes.xaxis.tick_top()
        labels = [item.get_text() for item in thisaxes.get_xticklabels()]
        empty_string_labels = ['']*len(labels)
        thisaxes.set_xticklabels(empty_string_labels)
    elif j == 0:
        thisaxes.set_xlabel('this is the x axis')
        thisaxes.xaxis.set_major_locator(majorLocator())
        thisaxes.xaxis.set_minor_locator(minorLocator())
    else:
        thisaxes.set_xticks([])
#}}}

# for the following, use the transforms

#{{{labeling ph1 for independent plots
fig.text(0.14,0.145,'ph1=0',ha='center',va='bottom',rotation='vertical')
fig.text(0.14,0.37,'ph1=1',ha='center',va='bottom',rotation='vertical')
fig.text(0.14,0.65,'ph1=0',ha='center',va='bottom',rotation='vertical')
fig.text(0.14,0.82,'ph1=1',ha='center',va='bottom',rotation='vertical')
#}}}
#{{{adding outer dimension labels
fig.text(0.07,0.25,'ph2=0',ha='center',va='bottom',rotation='vertical')
fig.text(0.07,0.725,'ph2=1',ha='center',va='bottom',rotation='vertical')
#}}}
#{{{Adding lines to separate ph1 and repeats
x1,y1 = list_of_axes[3].transAxes.transform(r_[0,1])
x15,y15 = list_of_axes[3].transAxes.transform(r_[0,0])
x2,y2 = list_of_axes[2].transAxes.transform(r_[0,1])
x25,y25 = list_of_axes[2].transAxes.transform(r_[0,0])
x3,y3 = list_of_axes[1].transAxes.transform(r_[0,1])
x35,y35 = list_of_axes[1].transAxes.transform(r_[0,0])
x4,y4 = list_of_axes[0].transAxes.transform(r_[0,1])
x45,y45 = list_of_axes[0].transAxes.transform(r_[0,0])
x1-=40
x_text = x1-10
x15-=40
x2-=40
x25-=40
x1,y1 = fig.transFigure.inverted().transform(r_[x1,y1])
x_text,_ = fig.transFigure.inverted().transform(r_[x_text,0])
x15,y15 = fig.transFigure.inverted().transform(r_[x15,y15])
x2,y2 = fig.transFigure.inverted().transform(r_[x2,y2])
x25,y25 = fig.transFigure.inverted().transform(r_[x25,y25])
x3-=40
x3_text = x3-10
x35-=40
x4-=40
x45-=40
x3,y3 = fig.transFigure.inverted().transform(r_[x3,y3])
x35,y35 = fig.transFigure.inverted().transform(r_[x35,y35])
x4_text,_ = fig.transFigure.inverted().transform(r_[x3_text,0])
x4,y4 = fig.transFigure.inverted().transform(r_[x4,y4])
x45,y45 = fig.transFigure.inverted().transform(r_[x45,y45])

lineA = lines.Line2D([x1,x15],[y1,y15],linestyle='--',
        linewidth=1.5, color='k', transform=fig.transFigure,
        clip_on=False)
lineA1 = lines.Line2D([x2,x25],[y2,y25],linestyle='--',
        linewidth=1.5,color='k',transform=fig.transFigure,
        clip_on=False)
lineC = lines.Line2D([x3,x35],[y3,y35],linestyle='--',
        linewidth=1.5,color='k',transform=fig.transFigure,
        clip_on=False)
lineC1 = lines.Line2D([x4,x45],[y4,y45],linestyle='--',
        linewidth=1.5,color='k',transform=fig.transFigure,
        clip_on=False)
#}}}
#{{{separating ph2 and ph1
x11,y11 = list_of_axes[3].transAxes.transform(r_[0,1])
x12,y12 = list_of_axes[2].transAxes.transform(r_[0,0])
x13,y13 = list_of_axes[1].transAxes.transform(r_[0,1])
x14,y14 = list_of_axes[0].transAxes.transform(r_[0,0])
x11-=65
x1_text = x11-10
x12-=65
x11,y11 = fig.transFigure.inverted().transform(r_[x11,y11])
x1_text,_ = fig.transFigure.inverted().transform(r_[x1_text,0])
x12,y12 = fig.transFigure.inverted().transform(r_[x12,y12])
x13-=65
x13_text = x13-10
x14-=65
x13,y13 = fig.transFigure.inverted().transform(r_[x13,y13])
x14_text,_ = fig.transFigure.inverted().transform(r_[x13_text,0])
x14,y14 = fig.transFigure.inverted().transform(r_[x14,y14])
lineB = lines.Line2D([x11,x12],[y11,y12], linestyle='--',linewidth=1.5, color='k', transform=fig.transFigure,
        clip_on=False)
lineD = lines.Line2D([x13,x14],[y13,y14], linestyle='--',
        linewidth=1.5, color='k',transform=fig.transFigure,
        clip_on=False)
#}}}
#{{{functions for brackets
def gen_curly(verts):
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', lw=2)
    fig.add_artist(patch)
    xs, ys = zip(*verts)
    # following shown for a guide -- don't use in final version
    fig.add_artist(lines.Line2D(xs, ys, linestyle='--', lw=2, color='black', ms=10,alpha=0.1))
def draw_bracket(x_left,y_middle,y_height,x_width=0.25):
    gen_curly([(x/0.25*x_width+x_left,y*y_height/2+y_middle) for x,y in verts])
    gen_curly([(x/0.25*x_width++x_left,-y*y_height/2+y_middle) for x,y in verts])
def add_label(x,y,thetext):
    # {{{ push out by 12 pts
    x,y = fig.transFigure.transform((x,y))
    x -= 12 
    x,y = fig.transFigure.inverted().transform((x,y))
    # }}}
    fig.text(x,y,thetext,
            va='center',ha='center',
            rotation=90,
            transform=fig.transFigure)
#}}}
#{{{drawing brackets on figure
draw_bracket(0.08,y4+0.05,0.315,x_width=0.04)
draw_bracket(0.08,y2+0.03,0.28,x_width=0.04)
draw_bracket(0.12,y45+0.05,0.12,x_width=0.02)
#}}}
show()

