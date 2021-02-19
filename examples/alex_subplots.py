from pyspecdata import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as lines
#{{{generating subplots with space between
fig=plt.figure(1)
ax1 = plt.subplot(611)
ax2 = plt.subplot(612)
ax3 = plt.subplot(614)
ax4 = plt.subplot(615)
fig.subplots_adjust(hspace=0.8)
list_of_axes = [ax1,ax2,ax3,ax4]
#}}}
#{{{adjusting tick marks and labels

yMajorLocator = lambda: mticker.MaxNLocator(steps=[1,2,5,10])
majorLocator = lambda: mticker.MaxNLocator(min_n_ticks=4, steps=[1,2,5,10])
minorLocator = lambda: mticker.AutoMinorLocator(n=5)
for j,thisaxes in enumerate(list_of_axes):
    thisaxes.set_ylabel('repeats')
    thisaxes.yaxis.set_minor_locator(minorLocator())
    thisaxes.yaxis.set_ticks_position('both')
    if j == 0:
        thisaxes.xaxis.set_major_locator(majorLocator())
        #for the minor ticks, use no labels; default NullFormatter
        thisaxes.xaxis.set_minor_locator(minorLocator())
        thisaxes.xaxis.tick_top()
        labels = [item.get_text() for item in thisaxes.get_xticklabels()]
        empty_string_labels = ['']*len(labels)
        thisaxes.set_xticklabels(empty_string_labels)
    elif j == len(list_of_axes)-1:
        thisaxes.set_xlabel('this is the x axis')
        thisaxes.xaxis.set_major_locator(majorLocator())
        thisaxes.xaxis.set_minor_locator(minorLocator())
    else:
        thisaxes.set_xticks([])
#}}}

# for the following, use the transforms

#{{{labeling ph1 for independent plots
fig.text(0.07,0.25,'ph1=0',ha='center',va='bottom',rotation='vertical')
fig.text(0.07,0.38,'ph1=1',ha='center',va='bottom',rotation='vertical')
fig.text(0.07,0.67,'ph1=0',ha='center',va='bottom',rotation='vertical')
fig.text(0.07,0.8,'ph1=1',ha='center',va='bottom',rotation='vertical')
fig.add_artist(lines.Line2D([0.08,0.08],[0.24,0.32],color='k'))
fig.add_artist(lines.Line2D([0.08,0.08],[0.37,0.45],color='k'))
fig.add_artist(lines.Line2D([0.08,0.08],[0.66,0.74],color='k'))
fig.add_artist(lines.Line2D([0.08,0.08],[0.79,0.87],color='k'))
#{{{adding outer dimension labels
fig.text(0.05,0.31,'ph2=0',ha='center',va='bottom',rotation='vertical')
fig.text(0.05,0.75,'ph2=1',ha='center',va='bottom',rotation='vertical')
fig.add_artist(lines.Line2D([0.06,0.06],[0.47,0.23],color='k'))
fig.add_artist(lines.Line2D([0.06,0.06],[0.87,0.66],color='k'))
show()

