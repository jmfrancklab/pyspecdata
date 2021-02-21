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
fig.text(0.08,0.25,'ph1=0',ha='center',va='bottom',rotation='vertical')
fig.text(0.08,0.38,'ph1=1',ha='center',va='bottom',rotation='vertical')
fig.text(0.08,0.67,'ph1=0',ha='center',va='bottom',rotation='vertical')
fig.text(0.08,0.8,'ph1=1',ha='center',va='bottom',rotation='vertical')
#}}}
#{{{adding outer dimension labels
fig.text(0.06,0.31,'ph2=0',ha='center',va='bottom',rotation='vertical')
fig.text(0.06,0.75,'ph2=1',ha='center',va='bottom',rotation='vertical')
#}}}
#{{{Adding lines to separate ph1 and repeats
x1,y1 = ax1.transAxes.transform(r_[0,1])
x15,y15 = ax1.transAxes.transform(r_[0,0])
x2,y2 = ax2.transAxes.transform(r_[0,1])
x25,y25 = ax2.transAxes.transform(r_[0,0])
x3,y3 = ax3.transAxes.transform(r_[0,1])
x35,y35 = ax3.transAxes.transform(r_[0,0])
x4,y4 = ax4.transAxes.transform(r_[0,1])
x45,y45 = ax4.transAxes.transform(r_[0,0])
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
lineA = lines.Line2D([x1,x15],[y1,y15],
        linewidth=3, color='k', transform=fig.transFigure,
        clip_on=False)
lineA1 = lines.Line2D([x2,x25],[y2,y25],
        linewidth=3,color='k',transform=fig.transFigure,
        clip_on=False)
lineC = lines.Line2D([x3,x35],[y3,y35],
        linewidth=3,color='k',transform=fig.transFigure,
        clip_on=False)
lineC1 = lines.Line2D([x4,x45],[y4,y45],
        linewidth=3,color='k',transform=fig.transFigure,
        clip_on=False)
#}}}
#{{{separating ph2 and ph1
x11,y11 = ax1.transAxes.transform(r_[0,1])
x12,y12 = ax2.transAxes.transform(r_[0,0])
x13,y13 = ax3.transAxes.transform(r_[0,1])
x14,y14 = ax4.transAxes.transform(r_[0,0])
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
lineB = lines.Line2D([x11,x12],[y11,y12], linewidth=3, color='k', transform=fig.transFigure,
        clip_on=False)
lineD = lines.Line2D([x13,x14],[y13,y14],
        linewidth=3, color='k',transform=fig.transFigure,
        clip_on=False)
#}}}
fig.add_artist(lineA)
fig.add_artist(lineA1)
fig.add_artist(lineB)
fig.add_artist(lineC)
fig.add_artist(lineC1)
fig.add_artist(lineD)

show()

