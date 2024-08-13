from ..general_functions import strm, process_kwargs, check_ascending_axis
from numpy import r_, c_, ix_, nan, pi
import numpy as np
from ..ndshape import ndshape_base as ndshape
from pylab import gca, sca, imshow, xlabel, ylabel, title, colorbar, setp
import logging


def image(A, x=[], y=[], allow_nonuniform=True, **kwargs):
    r"Please don't call image directly anymore -- use the image method of figurelist"
    x_inverted = False
    y_inverted = False
    A = A.copy()
    A.squeeze()  # drop any singleton dimensions, which cause problems
    # {{{ pull out kwargs for imagehsv
    imagehsvkwargs = {}
    for k, v in list(kwargs.items()):
        if k in ["black", "logscale", "scaling"]:
            imagehsvkwargs[k] = kwargs.pop(k)
    # }}}
    spacing, ax, x_first, origin, renumber = process_kwargs(
        [
            ("spacing", 1),
            ("ax", gca()),
            ("x_first", False),
            ("origin", "lower"),
            ("renumber", None),
        ],
        kwargs,
        pass_through=True,
    )
    sca(ax)
    if x_first:  # then the first dimension should be the column
        # dimesion (i.e. last)
        if hasattr(
            A, "dimlabels"
        ):  # if I try to use isinstance, I get a circular import
            new_dimlabels = list(A.dimlabels)
            temp = new_dimlabels.pop(0)
            A = A.reorder(new_dimlabels + [temp])
        else:
            A = A.T
    sca(ax)
    setlabels = False
    if hasattr(A, "dimlabels"):
        if renumber is not None:
            if type(renumber) is str:
                renumber = [renumber]
        for thisaxis in A.dimlabels:
            if renumber is not None and thisaxis in renumber:
                A.setaxis(thisaxis, "#")
            try:
                check_ascending_axis(
                    A.getaxis(thisaxis), allow_descending=True
                )
            except Exception:
                if allow_nonuniform:
                    logging.debug(
                        "Automatically changed to numbered axis along %s"
                        % thisaxis
                    )
                    A.setaxis(thisaxis, "#").set_units(thisaxis, "#")
                else:
                    raise ValueError(
                        "You are not allowed to use image on data that"
                        " doesn't have a uniformly spaced axis -- it is likely a"
                        " misrepresentation of the data you are looking at."
                        " For example, if you are looking at NMR data with a set of"
                        " variable delays that are unevenly spaced, relabel this axis"
                        " by index number --> .C.setaxis('%s','#').set_units('%s','scan"
                        " #').\nThen you have an accurate representation of your data"
                        % (2 * (thisaxis,))
                    )
        setlabels = True
        templabels = list(A.dimlabels)
        if A.get_prop("x_inverted"):
            x_inverted = True
        if A.get_prop("y_inverted"):
            y_inverted = True
        x_label = templabels[-1]
        if A.getaxis(x_label) is None:
            x = r_[0, ndshape(A)[x_label]]
        else:
            x = list(A.getaxis(x_label))
        x_label = A.unitify_axis(x_label)
        templabels.pop(-1)
        y_label = ""
        if len(templabels) == 1:
            y_label = templabels[0]
            try:
                y = list(A.getaxis(y_label))
            except Exception:
                y = r_[0 : A.data.shape[A.axn(y_label)]]
            y_label = A.unitify_axis(y_label)
        else:
            these_dimsizes = [str(ndshape(A)[x]) for x in templabels]

            def axis_labeler(
                x,
            ):  # below, I turn off math mode for the axis names
                if A.getaxis(x) is not None:
                    return (
                        "[$ "
                        + A.unitify_axis(x)
                        + r" $_{{{:.3g}\rightarrow{:.3g}}}]".format(
                            *(A.getaxis(x)[r_[0, -1]])
                        )
                    )
                # elif A.get_units(x) is not None:
                #    return '[$ '+A.unitify_axis(x)+ '$]'
                else:
                    return "[$ " + A.unitify_axis(x) + "$]"

            templabels = list(map(axis_labeler, templabels))
            y_label = "\\otimes".join(templabels)
            y_label = " _{(" + r"\times".join(these_dimsizes) + ")}" + y_label
            y_label = "$" + y_label + "$"  # whole expression is in math mode
        A = A.data
    if isinstance(x, list):
        x = np.array(x)
    if isinstance(y, list):
        y = np.array(y)
    if len(x) == 0:
        x = [1, A.shape[1]]
    else:
        x = x.flatten()
    if len(y) == 0:
        y = [1, A.shape[0]]
    else:
        y = y.flatten()
    dx = (x[-1] - x[0]) / len(x)
    dy = (y[-1] - y[0]) / len(y)
    if origin == "lower":
        myext = (
            x[0] - dx / 2.0,
            x[-1] + dx / 2.0,
            y[0] - dy / 2.0,
            y[-1] + dy / 2.0,
        )
    elif origin == "upper":
        myext = (
            x[0] - dx / 2.0,
            x[-1] + dx / 2.0,
            y[-1] + dy / 2.0,
            y[0] - dy / 2.0,
        )
    elif origin == "flip":
        # {{{ need to flip
        myext = (
            x[-1] + dx / 2.0,
            x[0] - dx / 2.0,
            y[-1] + dy / 2.0,
            y[0] - dy / 2.0,
        )
        # }}}
    else:
        raise ValueError(
            "I don't understand the value you've set for the origin keyword argument"
        )
    kwargs[
        "origin"
    ] = origin  # required so that imshow now displays the image correctly
    linecounter = 0
    if A.ndim > 2:
        setp(ax.get_yticklabels(), visible=False)
        ax.yaxis.set_ticks_position("none")
    while (
        A.ndim > 2
    ):  # to substitude for imagehsvm, etc., so that we just need a ersion of ft
        # order according to how it's ordered in the memory
        # the innermost two will form the image -- first add a line to the end of the images we're going to join up
        tempsize = np.array(A.shape)  # make a tuple the right shape
        if linecounter == 0 and spacing < 1.0:
            spacing = round(
                np.prod(tempsize[0:-1])
            )  # find the length of the thing not counting the columns
        tempsize[-2] = (
            2 * linecounter + spacing
        )  # all dims are the same except the image row, to which I add an increasing number of rows
        # print "iterate (A.ndim=%d) -- now linecounter is "%A.ndim,linecounter
        linecounter += tempsize[-2]  # keep track of the extra lines at the end
        A = np.concatenate(
            (A, nan * np.zeros(tempsize)), axis=(A.ndim - 2)
        )  # concatenate along the rows
        tempsize = r_[A.shape[0:-3], A.shape[-2:]]
        tempsize[-2] *= A.shape[-3]
        try:
            A = A.reshape(np.int64(tempsize))  # now join them up
        except Exception:
            raise IndexError(
                strm(
                    "problem with tempsize",
                    tempsize,
                    "of type",
                    type(tempsize),
                    "dtype",
                    tempsize.dtype,
                )
            )
    A = A[
        : A.shape[0] - linecounter, :
    ]  # really I should an extra counter besides linecounter now that I am using "spacing", but leave alone for now, to be sure I don't cut off data
    if origin == "flip":
        # {{{ if origin is "flip", we need to manually flip the data
        A = A[::-1, :]
        A = A[:, ::-1]
        kwargs["origin"] = "lower"
        # }}}
    if np.iscomplexobj(A):  # this just tests the datatype
        A = imagehsv(A, **imagehsvkwargs)
        retval = imshow(A, extent=myext, **kwargs)
    else:
        retval = imshow(A, extent=myext, **kwargs)
        colorbar()
    if setlabels:
        xlabel(x_label)
        # print y_label
        ylabel(y_label)
    if x_inverted:
        these_xlims = ax.get_xlim()
        ax.set_xlim((max(these_xlims), min(these_xlims)))
    if y_inverted:
        these_ylims = ax.get_ylim()
        ax.set_ylim((max(these_ylims), min(these_ylims)))
    return retval


def imagehsv(A, logscale=False, black=False, scaling=None):
    "This provides the HSV mapping used to plot complex number"
    # compare to http://www.rapidtables.com/convert/color/hsv-to-rgb.htm
    A = A.copy()
    n = 256
    mask = np.isnan(A)
    A[mask] = 0
    if scaling is None:
        A /= abs(A).max()
    else:
        A /= scaling
        mask |= abs(A) > 1.0 + 1e-7
        A[mask] = 0
    mask = mask.reshape(-1, 1)
    intensity = abs(A).reshape(-1, 1)
    if logscale:
        raise ValueError(
            "logscale is deprecated, use the cropped_log function instead"
        )
    # theta = (n-1.)*np.mod(np.angle(A)/pi/2.0,1)# angle in 255*cycles
    if black:
        if black is True:
            V = intensity
        else:
            V = intensity * black + (1.0 - black)
        S = 1.0  # always
    else:
        S = intensity
        V = 1.0  # always
    C = V * S
    H = (
        (np.angle(-1 * A).reshape(-1, 1) + pi) / 2.0 / pi * 6.0
    )  # divide into 60 degree chunks -- the -1 is to rotate so red is at origin
    X = C * (1 - abs(np.mod(H, 2) - 1))
    m = V - C
    colors = np.ones(list(A.shape) + [3])
    origshape = colors.shape
    colors = colors.reshape(-1, 3)
    rightarray = c_[C, X, np.zeros_like(X)]
    # http://en.wikipedia.org/wiki/HSL_and_HSV#From_HSV,
    # except that the order was messed up
    thismask = np.where(H < 1)[0]
    # C X 0
    colors[ix_(thismask, [0, 1, 2])] = rightarray[ix_(thismask, [0, 1, 2])]
    thismask = np.where(np.logical_and(H >= 1, H < 2))[0]
    # X C 0
    colors[ix_(thismask, [1, 0, 2])] = rightarray[ix_(thismask, [0, 1, 2])]
    thismask = np.where(np.logical_and(H >= 2, H < 3))[0]
    # X 0 C
    colors[ix_(thismask, [1, 2, 0])] = rightarray[ix_(thismask, [0, 1, 2])]
    thismask = np.where(np.logical_and(H >= 3, H < 4))[0]
    # 0 X C
    colors[ix_(thismask, [2, 1, 0])] = rightarray[thismask, :]
    thismask = np.where(np.logical_and(H >= 4, H < 5))[0]
    # 0 C X
    colors[ix_(thismask, [2, 0, 1])] = rightarray[thismask, :]
    thismask = np.where(H > 5)[0]
    # C 0 X
    colors[ix_(thismask, [0, 2, 1])] = rightarray[thismask, :]
    colors += m
    colors *= n - 1
    if black:
        # if the background is black, make the separators white
        # here, we have to remember that we're already scaled up to a scale of 0--255
        colors[mask * r_[True, True, True]] = 255.0
    else:
        # if the background is white, make the separators black
        colors[mask * r_[True, True, True]] = 0.0
    colors = colors.reshape(origshape)
    return np.uint8(colors.round())


def fl_image(self, A, **kwargs):
    r"""Called as `fl.image()` where `fl` is the `figlist_var`
    object

    Note that this code just wraps the figlist properties, and
    the heavy lifting is done by the `image(` function.
    Together, the effect is as follows:

    - `check_units` converts to human-readable units, and
      makes sure they match the units already used in the plot.
    - if `A` has more than two dimensions, the final dimension in
      `A.dimlabels` is used as the column dimension, and a
      direct-product of all non-column dimensions (a Kronecker
      product, such that the innermost index comes the latest in
      the list `A.dimlabels`) is used as the row dimension. A
      white/black line is drawn after the innermost index used to
      create the direct product is finished iterating.
    - If `A` consists of complex data, then an HSV plot
      (misnomer, actually an HV plot) is used:
      - convert to polar form: :math:`z=\rho \exp(i \phi)`
      - :math:`\phi` determines the color (Hue)
        - Color wheel is cyclical, like :math:`\exp(i \phi)`
        - red is taken as :math:`\phi=0`, purely real and positive
        - green-blue is :math:`pi` radians out of phase with red and
          therefore negative real
      - :math:`\rho` determines the intensity (value)
        - Depending on whether or not `black` is set (either as a
          keyword argument, or `fl.black`, the background will be
          black with high :math:`\rho` values "lit up" (intended for
          screen plotting) or the background will be white with
          the high :math:`\rho` values "colored in" (intended for
          printing)
    - If the data type (`dtype`) of the data in `A` is real
      (typically achieved by calling `abs(A)` or
      `A.runcopy(real)`), then `A` is plotted with a colormap and
      corresponding colorbar.
    - If no title has been given, it's set to the name of the
      current plot in the figurelist

    Attributes
    ----------
    A : nddata or numpy array
    x : Optional[double] or Optional[scalar]
        If `A` is a numpy array, then this gives the values along
        the x axis (columns).
        Defaults to the size of the array.
        Not used if `A` is `nddata`.
    y : Optional[double] or Optional[scalar]
        If `A` is a numpy array, then this gives the values along
        the y axis (columns).
        Defaults to the size of the array.
        Not used if `A` is `nddata`.
    x_first : boolean
        Since it's designed to represent matrices, an image plot
        by defaults is "transposed" relative to all other plots.
        If you want the first dimension on the x-axis (*e.g.*, if
        you are plotting a contour plot on top of an image), then set
        `x_first` to `True`.
    spacing : integer
        Determines the size of the white/black line drawn
        Defaults to 1
    ax : matplotlib Axes
        the Axis object where the plot should go.
    all remaning :
        are passed through to matplotlib `imshow`
    origin : {'upper', 'lower', 'flip'}
        upper and lower are passed to matplotlib.
        Flip is for 2D nmr, and flips the data manually.

    .. code-block:: python
        from pyspecdata import *
        fl = figlist_var()

        t = r_[-1:1:300j]
        x = nddata(t,[-1],['x']).labels('x',t)
        y = nddata(t,[-1],['y']).labels('y',t)

        z = x**2 + 2*y**2
        print "dimlabels of z:",z.dimlabels

        fl.next('image with contours')
        fl.image(z,x_first = True) #  x_first is needed to align
        #                             with the contour plot
        z.contour(colors = 'w',alpha = 0.75)

        fl.next('simple plot') #  just to show that x is the same
        #                         here as well
        fl.plot(z['y':(0,0.01)])

        fl.show('compare_image_contour_150911.pdf')
    """
    interpolation, ax, human_units, row_threshold = process_kwargs(
        [
            ("interpolation", None),
            ("ax", gca()),
            ("human_units", True),
            ("row_threshold", 500),
        ],
        kwargs,
        pass_through=True,
    )
    sca(ax)
    if human_units:
        firstarg = self.check_units(
            A, -1, 0
        )  # check units, and if need be convert to human units, where x is the last dimension and y is the first
    else:
        firstarg = A
    if self.black and "black" not in list(kwargs.keys()):
        kwargs.update({"black": self.black})
    retval = image(
        firstarg, **kwargs
    )  # just a placeholder for now, will later keep units + such
    if ax.get_title() is None or len(ax.get_title()) == 0:
        title(self.current)
    if interpolation is not None:
        if interpolation == "auto":
            if (
                np.prod([firstarg.shape[j] for j in firstarg.dimlabels[:-1]])
                > row_threshold
            ):
                interpolation = "bilinear"
            else:
                interpolation = "nearest"
        retval.set_interpolation(interpolation)
    return retval
