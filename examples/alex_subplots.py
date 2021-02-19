from pyspecdata import *
import matplotlib.pyplot as plt
import numpy as np
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
ax1.yaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')
ax3.yaxis.set_ticks_position('both')
ax4.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('top')
#}}}
#{{{x axis label is same for all plots
ax4.set_xlabel('this is the x axis')
#}}}
#{{{labeling inner dimensions for independent plots
ax4.set_ylabel('ph1=0\nrepeats')
ax3.set_ylabel('ph1=1\nrepeats')
ax2.set_ylabel('ph1=0\nrepeats')
ax1.set_ylabel('ph1=1\nrepeats')

#}}}
#{{{adding outer dimension labels
fig.text(0.05,0.25,'ph2=0\n------------------------',ha='center',va='bottom',rotation='vertical')
fig.text(0.05,0.65,'ph2=1\n------------------------',ha='center',va='bottom',rotation='vertical')

show()

