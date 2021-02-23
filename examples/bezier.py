import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.lines as lines

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

fig = plt.figure()
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
draw_bracket(0.2,0.5,0.6,x_width=0.1)
add_label(0.2,0.5,r'$\Delta p_1=0$')
draw_bracket(0.4,0.5+0.3,0.2,x_width=0.05)
add_label(0.4,0.5+0.3,r'$\Delta p_2=0$')
draw_bracket(0.4,0.5-0.3,0.2,x_width=0.05)
add_label(0.4,0.5-0.3,r'$\Delta p_2=1$')
plt.show()
