"""
transformation math
===================
"""
from pylab import *
from matplotlib import lines
from matplotlib.transforms import IdentityTransform
fig, (ax1,ax2) = subplots(2,1, figsize=(4,4))
subplots_adjust(left=0.2)
x1,y1 = ax1.transAxes.transform(r_[0,1])
x2,y2 = ax2.transAxes.transform(r_[0,0])
x1-=40
x_text = x1-10
x2-=40
x1,y1 = fig.transFigure.inverted().transform(r_[x1,y1])
x_text,_ = fig.transFigure.inverted().transform(r_[x_text,0])
x2,y2 = fig.transFigure.inverted().transform(r_[x2,y2])
lineA = lines.Line2D([x1,x2],[y1,y2],
        linewidth=3, color='r', transform=fig.transFigure,
        clip_on=False)
lineB = lines.Line2D([0,1],[0,1], linewidth=3, color='b', transform=ax1.transAxes,
        clip_on=False)
text(x_text, 0.5, "a label", va='center', ha='right', rotation=90, transform=fig.transFigure, color='r')
ax1.add_line(lineA)
ax1.add_line(lineB)
show()
