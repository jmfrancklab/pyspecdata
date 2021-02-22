from pyspecdata import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as lines
#import matplotlib.patches
##import matplotlib.patches.ArrowStyle
#import matplotlib.patches.ArrowStyle.get_styles
#import matplotlib.axes.Axes.annotate
#{{{generating subplots with space between
bottoms = [0.13, 0.35, 0.63, 0.8]
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
fig.text(0.12,0.145,'ph1=0',ha='center',va='bottom',rotation='vertical')
fig.text(0.12,0.37,'ph1=1',ha='center',va='bottom',rotation='vertical')
fig.text(0.12,0.65,'ph1=0',ha='center',va='bottom',rotation='vertical')
fig.text(0.12,0.82,'ph1=1',ha='center',va='bottom',rotation='vertical')
#}}}
#{{{adding outer dimension labels
fig.text(0.05,0.26,'ph2=0',ha='center',va='bottom',rotation='vertical')
fig.text(0.05,0.75,'ph2=1',ha='center',va='bottom',rotation='vertical')
#}}}
plt.annotate(r"$\{$",fontsize=74,
        xy=(0.125,0.15),xycoords='figure fraction')
plt.annotate(r"$\{$",fontsize=74,
        xy=(0.125,0.365),xycoords='figure fraction')
plt.annotate(r"$\{$",fontsize=160,
        xy=(0.05,0.22),xycoords='figure fraction')
plt.annotate(r"$\{$",fontsize=74,
        xy=(0.125,0.64),xycoords='figure fraction')
plt.annotate(r"$\{$",fontsize=74,
        xy=(0.125,0.815),xycoords='figure fraction')
plt.annotate(r"$\{$",fontsize=150,
        xy=(0.05,0.7),xycoords='figure fraction')
show()

