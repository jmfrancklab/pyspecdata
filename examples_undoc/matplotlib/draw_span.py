"""
draw a span
===========
"""
from pylab import *
from matplotlib import lines
from matplotlib.transforms import IdentityTransform
fig, (ax1,ax2,ax3,ax4) = subplots(4,1, figsize=(4,4))
subplots_adjust(left=0.2)
def draw_span(ax1, ax2, label, allow_for_text=10, allow_for_ticks=40):
    x1,y1 = ax1.transAxes.transform(r_[0,1])
    x2,y2 = ax2.transAxes.transform(r_[0,0])
    x1-=allow_for_ticks
    x_text = x1-allow_for_text
    x2-=allow_for_ticks
    x1,y1 = fig.transFigure.inverted().transform(r_[x1,y1])
    x_text,_ = fig.transFigure.inverted().transform(r_[x_text,0])
    x2,y2 = fig.transFigure.inverted().transform(r_[x2,y2])
    lineA = lines.Line2D([x1,x2],[y1,y2],
            linewidth=3, color='r', transform=fig.transFigure,
            clip_on=False)
    text(x_text, (y2+y1)/2, label, va='center', ha='right', rotation=90, transform=fig.transFigure, color='r')
    fig.add_artist(lineA)
draw_span(ax1,ax3,"a label")
show()
