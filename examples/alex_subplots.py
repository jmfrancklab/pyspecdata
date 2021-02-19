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
#}}}
#{{{adjusting tick marks and labels
labels = [item.get_text() for item in ax1.get_xticklabels()]
empty_string_labels = ['']*len(labels)
ax1.set_xticklabels(empty_string_labels)
ax2.set_xticks([])
ax3.set_xticks([])
defaultMajorLocator = lambda: mticker.MaxNLocator(min_n_ticks=4, steps=[1,2,5,10])
defaultMinorLocator = lambda: mticker.AutoMinorLocator(n=5)
ax4.xaxis.set_major_locator(defaultMajorLocator())
#for the minor ticks, use no labels; default NullFormatter
ax4.xaxis.set_minor_locator(defaultMinorLocator())

y1MajorLocator = lambda: mticker.MaxNLocator(min_n_ticks=4, steps=[1,2,5,10])
y1MinorLocator = lambda: mticker.AutoMinorLocator(n=5)
ax1.yaxis.set_major_locator(y1MajorLocator())
ax1.yaxis.set_minor_locator(y1MinorLocator())
ax1.yaxis.set_ticks_position('both')
y2MajorLocator = lambda: mticker.MaxNLocator(min_n_ticks=4, steps=[1,2,5,10])
y2MinorLocator = lambda: mticker.AutoMinorLocator(n=5)
ax2.yaxis.set_major_locator(y2MajorLocator())
ax2.yaxis.set_minor_locator(y2MinorLocator())
ax2.yaxis.set_ticks_position('both')
y3MajorLocator = lambda: mticker.MaxNLocator(min_n_ticks=4, steps=[1,2,5,10])
y3MinorLocator = lambda: mticker.AutoMinorLocator(n=5)
ax3.yaxis.set_major_locator(y3MajorLocator())
ax3.yaxis.set_minor_locator(y3MinorLocator())
ax3.yaxis.set_ticks_position('both')
y4MajorLocator = lambda: mticker.MaxNLocator(min_n_ticks=4, steps=[1,2,5,10])
y4MinorLocator = lambda: mticker.AutoMinorLocator(n=5)
ax4.yaxis.set_major_locator(y4MajorLocator())
ax4.yaxis.set_minor_locator(y4MinorLocator())
ax4.yaxis.set_ticks_position('both')
majorLocator = lambda: mticker.MaxNLocator(min_n_ticks=4, steps=[1,2,5,10])
minorLocator = lambda: mticker.AutoMinorLocator(n=5)
ax1.xaxis.set_major_locator(majorLocator())
#for the minor ticks, use no labels; default NullFormatter
ax1.xaxis.set_minor_locator(minorLocator())
ax1.xaxis.tick_top()


#}}}
#{{{x axis label is same for all plots
ax4.set_xlabel('this is the x axis')
#}}}
#{{{labeling repeats for independent plots
ax4.set_ylabel('repeats')
ax3.set_ylabel('repeats')
ax2.set_ylabel('repeats')
ax1.set_ylabel('repeats')

#}}}
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

