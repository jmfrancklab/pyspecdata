The Figure List
===============

So that we can place figures either in an organized lab
notebook,
or in a graphical display (*e.g.* Qt-based),
or in simple figures that show up at the command line,
we need a means of collating figures
(as well as intervening text)
into an organized
list -- this is the figure list.

Previously, this was simply a list of matplotlib or
mayavi figures
(typically implemented with ``fl.next('figure name')`` calls),
but that strategy proved to be insufficient for various
reasons.
Practically, this leads to complicated code when
implementing, *e.g.* twinx or twiny plots,
and it also obscures the underlying matplotlib, mayavi,
etc, code from the user.
Philosophically, as with nddata, we want to take maximal
advantage of operator overloading and other object
oriented concepts, including seamless switching between
output to, *e.g.* latex, html, bokeh, matplotlib, Qt
GUIs, *etc.*

Therefore, we build the new system around a list of
*visualization objects*.
By breaking the plotting process into several methods
which are defined in **create new document with definitions**,
we can automate the process,
while also keeping it maximally flexible.
In the simplest example,
a visualization object corresponds to *e.g.* a single
matplotlib axis;
but it can also correspond to multiple axes,
as when we want to plot multi-dimensional data, or when
we want to show the projections of 2D data.

From the end-user's perspective, plotting is achieved
by
(1) optionally modifying an existing plotting object
class to get different properties than the default:
the default classes available are
image,
plot,
contour,
and
surface
(2) creating instances of the plotting classes inside
inside the figure list,
in order of appearance
(note that sometimes, it will make sense to first
initialize blank plots, and drop data into them if the
order of plotting and data generation are different)
(3) adding nddata objects (or numpy objects, which are
converted to nddata objects) to the plotting instances
(quite literally)
(4) if interactive plotting is employed, the nddata
themselves are updated, and then the figurelist update
method is called.
