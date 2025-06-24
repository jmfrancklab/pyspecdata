import numpy as np
from numpy import r_, pi
from matplotlib import ticker as mticker
import matplotlib.pylab as plt
import matplotlib.transforms as mtransforms
import matplotlib.ticker as mticker
import itertools
#{{{ old grid and tick
def gridandtick(ax,rotation=(0,0),precision=(2,2),
        labelstring=('',''),gridcolor=r_[0,0,0],
        formatonly = False,fixed_y_locator = None,
        use_grid = True,
        spines = None,y = True):
    defaultMajorLocator = lambda: mticker.MaxNLocator(min_n_ticks=4, steps=[1,2,5,10])
    defaultMinorLocator = lambda: mticker.AutoMinorLocator(n=5)
    #{{{ taken from matplotlib examples
    def adjust_spines(ax,spines):
        xlabel = ax.get_xlabel()
        ylabel = ax.get_ylabel()
        for loc, spine in list(ax.spines.items()):
            if loc in spines:
                spine.set_position(('outward',5)) # outward by 5 points
                spine.set_smart_bounds(True)
            else:
                spine.set_color('none') # don't draw spine
        # turn off ticks where there is no spine
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            # no yaxis ticks
            ax.yaxis.set_ticks([],minor = False)
        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            # no xaxis ticks
            ax.xaxis.set_ticks([],minor = False)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    #}}}
    if spines is not None:
        adjust_spines(ax,spines = spines)
    if not formatonly:
        #{{{x ticks
        # determine the size
        width = abs(np.diff(ax.get_xlim()))
        if width==0:
            raise ValueError('x axis width is zero')
        widthexp = np.floor(np.log(width)/np.log(10.))-1
        scalefactor = 10**widthexp
        width /= scalefactor
        majorLocator = defaultMajorLocator()
        #majorFormatter = FormatStrFormatter('%0.'+'%d'%precision[0]+'f'+labelstring[0])# labelstring can be used, for instance, for pi
        #ax.xaxis.set_major_formatter(majorFormatter)
        minorLocator = defaultMinorLocator()
        ax.xaxis.set_major_locator(majorLocator)
        #for the minor ticks, use no labels; default NullFormatter
        ax.xaxis.set_minor_locator(minorLocator)
        #}}}
        if y:
            logarithmic = True if ax.get_yaxis().get_scale() == 'log' else False
            #{{{ y ticks
            width = abs(np.diff(ax.get_ylim()))
            if width==0:
                raise ValueError('y axis width is zero')
            widthexp = np.floor(np.log(width)/np.log(10.))-1
            scalefactor = 10**widthexp
            width /= scalefactor
            if fixed_y_locator is None:
                if logarithmic:
                    print("logarithmic")
                    majorLocator = mticker.LogLocator(10)
                    minorLocator = mticker.LogLocator(10,subs=r_[0:11])
                else:
                    majorLocator = defaultMajorLocator()
                    minorLocator = defaultMinorLocator()
            else:
                majorLocator = mticker.MultipleLocator(fixed_y_locator[4::5])
                minorLocator = mticker.FixedLocator(fixed_y_locator)
            #majorFormatter = FormatStrFormatter('%0.'+'%d'%precision[1]+'f'+labelstring[1])# labelstring can be used, for instance, for pi
            #ax.yaxis.set_major_formatter(majorFormatter)
            ax.yaxis.set_major_locator(majorLocator)
            #for the minor ticks, use no labels; default NullFormatter
            ax.yaxis.set_minor_locator(minorLocator)
            #}}}
    if use_grid:
        ax.yaxis.grid(use_grid,which='major',color=gridcolor,alpha=0.15,linestyle='-')
        ax.xaxis.grid(use_grid,which='major',color=gridcolor,alpha=0.15,linestyle='-')
        ax.yaxis.grid(use_grid,which='minor',color=gridcolor,alpha=0.075,linestyle='-')
        ax.xaxis.grid(use_grid,which='minor',color=gridcolor,alpha=0.075,linestyle='-')
    labels = ax.get_xticklabels()
    plt.setp(labels,rotation=rotation[0])
    if y:
        labels = ax.get_yticklabels()
        plt.setp(labels,rotation=rotation[1])
    fig = plt.gcf()
    fig.autofmt_xdate()
    return
def gridon(gridcolor=r_[0,0,0]):
    plt.grid(True,which='major',color=gridcolor,alpha=0.1,linestyle='-')
    plt.grid(True,which='minor',color=gridcolor,alpha=0.05,linestyle='-')
#}}}
#{{{ a better version?
def othergridandtick(ax,rotation=(0,0),precision=(2,2),labelstring=('',''),gridcolor=r_[0,0,0],y = True,x = True,spines = None):
    #{{{ taken from matplotlib examples
    def adjust_spines(ax,spines):
        xlabel = ax.get_xlabel()
        ylabel = ax.get_ylabel()
        for loc, spine in list(ax.spines.items()):
            if loc in spines:
                spine.set_position(('outward',5)) # outward by 5 points
                spine.set_smart_bounds(True)
            else:
                spine.set_color('none') # don't draw spine
        # turn off ticks where there is no spine
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            # no yaxis ticks
            ax.yaxis.set_ticks([],minor = False)
        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            # no xaxis ticks
            ax.xaxis.set_ticks([],minor = False)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    #}}}
    if spines is not None:
        adjust_spines(plt.gca(),spines = spines)
    if x:
        #{{{x ticks
        # determine the size
        ax.xaxis.set_major_locator(mticker.MaxNLocator(10)) # could use multiplelocator if it keeps try to do multiples of 2
        ax.xaxis.set_minor_locator(mticker.MaxNLocator(50))
        #}}}
    if y:
        #{{{ y ticks
        ax.yaxis.set_major_locator(mticker.MaxNLocator(10))
        ax.yaxis.set_minor_locator(mticker.MaxNLocator(50))
        #}}}
    plt.grid(True,which='major',color=gridcolor,alpha=0.2,linestyle='-')
    plt.grid(True,which='minor',color=gridcolor,alpha=0.1,linestyle='-')
    if x:
        labels = ax.get_xticklabels()
        plt.setp(labels,rotation=rotation[0])
    if y:
        labels = ax.get_yticklabels()
        plt.setp(labels,rotation=rotation[1])
    return
#}}}
def autolegend(*args,**kwargs):
    #lg = legend(legendstr,'best'),loc = 2, borderaxespad = 0.)
    match_colors = False
    if 'match_colors' in list(kwargs.keys()):
        match_colors = kwargs.pop('match_colors')
    alpha = 0.45
    if 'alpha' in list(kwargs.keys()):
        alpha = kwargs.pop('alpha')
    if 'ax' in list(kwargs.keys()):
        ax_list = [kwargs.pop('ax')]
    else:
        ax_list = [plt.gca()]
    if 'ax2' in list(kwargs.keys()):
        ax_list.append(kwargs.pop('ax2'))
    for ax in ax_list:
        if len(args)==0:
            lg = ax.legend(**kwargs)
        elif len(args)==1:
            lg = ax.legend(args[0],**kwargs)
        else:
            lg = ax.legend(args[0],args[1],**kwargs)
        if lg is None:
            raise ValueError("Warning! you called autolegend, but you don't seem to have anything labeled!!")
        else:
            lg.get_frame().set_alpha(alpha)
    if match_colors:
        for line, txt in zip(lg.get_lines(), lg.get_texts()): # from http://stackoverflow.com/questions/13828246/matplotlib-text-color-code-in-the-legend-instead-of-a-line 
                    txt.set_color(line.get_color())  
                    txt.set_alpha(line.get_alpha())  
    return lg
def autopad_figure(pad = 0.2,centered = False,figname = 'unknown'):
    #{{{ solve the axis issue --> this does just the left
    fig = plt.gcf()
    ax = plt.gca()
    labelsets = [] 
    #labelsets.append(('left',ax.get_yticklabels()))
    #labelsets.append(('left',ax.get_yticklines()))
    #labelsets.append(('right',ax.get_yticklines()))
    labelsets.append(('left',[plt.ylabel(ax.get_ylabel())]))
    #labelsets.append(('bottom',ax.get_xticklabels()))
    #labelsets.append(('bottom',ax.get_xticklines()))
    if len(ax.get_xlabel()) > 0:
        labelsets.append(('bottom',[plt.xlabel(ax.get_xlabel())]))
    #labelsets.append(('top',ax.get_xticklines()))
    if len(ax.get_title()) > 0:
        pass #labelsets.append(('top',[plt.title(ax.get_title())]))
    compto = {}
    def on_draw(event):
        # find the sum of the widths of all things labeled with a 'y'
        spkwargs = {}
        compto['bottom'] = fig.subplotpars.bottom
        compto['left'] = fig.subplotpars.left
        compto['right'] = fig.subplotpars.right
        compto['top'] = fig.subplotpars.top
        for axisn in ['left','bottom','top','right']:
            bboxes = []
            labellist = [x[1] for x in labelsets if x[0] is axisn]
            for labels in labellist:
                for label in labels:
                    if isinstance(label, plt.Line2D):
                        pass # just rely on the pad
                        #if np.any(map(lambda x: x == label.get_transform(),[ax.transData,ax.transAxes,fig.transFigure,None])):
                        #    print 'found it'
                        #else:
                        #    print 'didn not find it'
                        #bbox = label.get_window_extent(fig.canvas).inverse_transformed(ax.transData).inverse_transformed(fig.transFigure)
                    else:
                        try:
                            bbox = label.get_window_extent()
                        except Exception as e:
                            warnings.warn("I wasn't able to run autopad on figure"+figname+"\nGetting window extent throws error"+str(e))
                    # the figure transform goes from relative coords->pixels and we
                    # want the inverse of that
                    bboxes.append(bbox)
                # this is the bbox that bounds all the bboxes, again in relative
                # figure coords
            l = 0 
            if len(labellist):
                bbox = mtransforms.Bbox.union(bboxes)
                bboxi = bbox.transformed(fig.transFigure.inverted())
                if axisn in ['left','right']:
                    l = bboxi.width
                if axisn in ['top','bottom']:
                    l = bboxi.height
            l += pad
            if axisn in ['top','right']:
                l = 1-l
                if compto[axisn] > l:
                    spkwargs.update({axisn:l})
            else:
                if compto[axisn] < l:
                    spkwargs.update({axisn:l})
        if len(spkwargs) > 0:
            if centered and 'left' in list(spkwargs.keys()) and 'right' in list(spkwargs.keys()):
                big = max(r_[spkwargs['left'],1-spkwargs['right']])
                spkwargs.update({'left':big,'right':1-big})
            try:
                fig.subplots_adjust(**spkwargs) # pad a little
            except:
                raise RuntimeError('failed to adjust subplots spwargs = ',spkwargs)
            #print "adjusted to",spkwargs
            fig.canvas.draw()# recurse
        return False
    fig.canvas.mpl_connect('draw_event', on_draw)
    fig.subplots_adjust(left = 0, right = 1, top = 1, bottom =0)
    fig.canvas.draw()# it needs this to generate the 'renderers'
    fig.canvas.mpl_connect('draw_event', on_draw)
    fig.canvas.draw()
    return
    #}}}
def expand_x(*args):
    r'''expand the axes.  If an argument is passed, then it refers to the position relative to the current coordinates.  Values can be:
        :0: set this side of the axis to 0
        :None: leave this side of the axis alone
        :a double: rescale the distance from the center of the axis to this side by this number'''
    # this is matplotlib code to expand the x axis
    ax = plt.gca()
    xlims = np.array(ax.get_xlim())
    width = abs(np.diff(xlims))
    thismean = np.mean(xlims)
    if len(args) > 0:
        if len(args) == 1 and isinstance(args, tuple):
            args = args[0]
        for j in range(2):
            if args[j] is None:
                pass
            elif args[j] == 0:
                xlims[j] = 0
            else:
                xlims[j] = args[j]*(xlims[j]-thismean) + thismean
    else:
        xlims[0] -= width/10
        xlims[1] += width/10
    ax.set_xlim(xlims)
def expand_y(*args):
    r'''expand the axes.  If an argument is passed, then it refers to the position relative to the current coordinates.  Values can be:
        :0: set this side of the axis to 0
        :None: leave this side of the axis alone
        :a double: rescale the distance from the center of the axis to this side by this number'''
    # this is matplotlib code to expand the x axis
    ax = plt.gca()
    ylims = np.array(ax.get_ylim())
    width = abs(np.diff(ylims))
    thismean = np.mean(ylims)
    if len(args) > 0:
        if len(args) == 1 and isinstance(args, tuple):
            args = args[0]
        for j in range(2):
            if args[j] is None:
                pass
            elif args[j] == 0:
                ylims[j] = 0
            else:
                ylims[j] = args[j]*(ylims[j]-thismean) + thismean
    else:
        ylims[0] -= width/10
        ylims[1] += width/10
    ax.set_ylim(ylims)
def plot_label_points(x,y,labels,**kwargs_passed):
    kwargs = {'alpha':0.5,'color':'g','ha':'left','va':'center','rotation':0,'size':14}
    kwargs.update(kwargs_passed)
    for j in range(0,len(labels)):
        text(x[j],y[j],labels[j],**kwargs)
def addlabels(labelstring,x,y,labels):
    r'obsolete -- use plot_label_points'
    for j in range(0,len(labels)):
        text(x[j],y[j],labelstring%labels[j],alpha=0.5,color='g',ha='left',va='top',rotation=0)
def plot_color_counter(*args,**kwargs):
    """Try not to use this function any more -- the version-to-version support for capturing and setting color cycles in matplotlib is very very bad.  (And, the cycler object in newer versions of matplolib is confusing.) So, just import `cycle` from `itertools`, and use it to build a cycle that you directly call to set your properties.

    .. note::
        previous description:

        if passed an argument: make it so that the next line will have the properties given by the argument

        if not passed an argument: just return the current plot properties,so that I can cycle back to it"""
    ax = process_kwargs([('ax',plt.gca())],kwargs)
    if len(args)>0:
        if LooseVersion(matplotlib.__version__) >= LooseVersion("1.5"):
            # {{{ find the element before the one we want
            retval = args[0]
            penultimate = next(ax._get_lines.prop_cycler)
            j = next(ax._get_lines.prop_cycler)
            not_in_list_counter = 1000.
            while j != args[0]:
                penultimate = j
                j = next(ax._get_lines.prop_cycler)
                not_in_list_counter -= 1
                if not_in_list_counter == 0:
                    raise ValueError("the value isn't in the cycler!")
            # }}}
            # {{{ now, set to the element before
            not_in_list_counter = 1000.
            while j != penultimate:
                j = next(ax._get_lines.prop_cycler)
                not_in_list_counter -= 1
                if not_in_list_counter == 0:
                    raise ValueError("the value isn't in the cycler!")
            # }}}
        else:
            try:
                ax._get_lines.count = args[0] # set the value of the color counter
            except:
                ax._get_lines.color_cycle = args[0] # set the value of the color counter
            retval = args[0]
    else:
        if LooseVersion(matplotlib.__version__) >= LooseVersion("1.5"):
            # {{{ I want to return the current element of the cycle
            one_too_far = next(ax._get_lines.prop_cycler)
            j = next(ax._get_lines.prop_cycler)
            not_in_list_counter = 1000.
            while j != one_too_far:
                penultimate = j
                j = next(ax._get_lines.prop_cycler)
                not_in_list_counter -= 1
                if not_in_list_counter == 0:
                    raise ValueError("the value isn't in the cycler!")
            retval = penultimate
            # }}}
        else:
            try: # this is different depending on the version of.core
                retval = ax._get_lines.count
            except:
                retval = ax._get_lines.color_cycle
    return retval
def contour_plot(xvals,yvals,zvals,color = 'k',alpha = 1.0,npts = 300,**kwargs):
    if 'inline_spacing' in list(kwargs.keys()):
        inline_spacing = kwargs.pop('inline_spacing')
    else:
        inline_spacing = 20
    xi = np.linspace(xvals.min(),xvals.max(),npts)
    yi = np.linspace(yvals.min(),yvals.max(),npts)
    #{{{ show the diffusivity
    #plot(np.array(xvals),np.array(yvals),'k')# to show where everything is
    zi = scipy_griddata((xvals,yvals),
        zvals,
        (xi[None,:],yi[:,None]))
    zi_min = zi[np.isfinite(zi)].min()
    zi_max = zi[np.isfinite(zi)].max()
    levels = r_[zi_min:zi_max:40j]
    CS = plt.contour(xi,yi,zi,levels,colors = color,
            alpha = 0.25*alpha)
    oldspacing = levels[1]-levels[0]
    levels = r_[zi_min:zi_max:oldspacing*5]
    try:
        CS = plt.contour(xi,yi,zi,levels,colors = color,
            alpha = alpha,**kwargs)
    except Exception as e:
        raise Exception(strm("Is there something wrong with your levels?:",levels,"min z",zi_min,"max z",zi_max,explain_error(e)))
    plt.clabel(CS,inline = 1,
        #fmt = r'$k_\sigma/k_{\sigma,bulk} = %0.2f$',
        fmt = r'%0.2f',
        use_clabeltext = True,
        inline_spacing = inline_spacing,
        )
    #}}}
def plot_updown(data,axis,color1,color2,symbol = '',**kwargs):
    if symbol == '':
        symbol = 'o'
    change = r_[1,np.diff(data.getaxis(axis))]
    changemask = change > 0
    if 'force_color' in list(kwargs.keys()) and kwargs['force_color'] == True:
        if hasattr(data,'other_info'):
            if 'plot_color' in data.get_prop():
                data.other_info.pop('plot_color')
    plot(data[axis,changemask],color1+symbol,**kwargs)
    if len(kwargs) > 0 and 'label' in list(kwargs.keys()): kwargs.pop('label') # if I'm doing a legend, I want it on the first
    plot(data[axis,~changemask],color2+symbol,**kwargs)
    return
def nextfigure(figurelist,name):
    'obsolete -- now use class'
    if isinstance(figurelist,figlist_var):
        figurelist.next(name)
        return figurelist
    else:
        print('Boo! not a new style name!')
    logger.debug(strm(lsafe('DEBUG figurelist, called with',name)))
    if name in figurelist:
        fig = figure(figurelist.index(name)+1)
        logger.debug(strm(lsafen('in',figurelist,'at figure',figurelist.index(name)+1,'switched figures')))
    else:
        fig = figure(len(figurelist)+1)
        fig.add_subplot(111)
        logger.debug(strm(lsafen('added, figure',len(figurelist)+1,'because not in figurelist',figurelist)))
        figurelist.append(name)
    return figurelist
def figlistret(first_figure,figure_list,*args,**kwargs):
    if 'basename' in list(kwargs.keys()):
        basename = kwargs['basename']
    else:
        basename = thisjobname()
    if first_figure is None:
        figure_list.show(basename+'.pdf')
        return args
    else:
        args += (figure_list,)
        if len(args) == 1:
            return args[0]
        else:
            return args
def figlistini(first_figure):
    r"""processes a figure list argument:
    typically, you want to have a figure_list keyword argument for every function, which is by default set to None, then call this on the argument -- it always returns a figure list, creating a new one if required
    similarly, somewhere I have another guy that processes the output, so that if it's set to None, it will by default dump and show the figure list,
    and not return a figure list in the output"""
    if first_figure is None:
        return figlist_var() 
    else:
        return first_figure
def figlistini_old(first_figure):
    if isinstance(first_figure,figlist_var):
        return first_figure
    else:
        print("Boo, not a new style name! (initialize)")
    logger.debug(strm(lsafe('DEBUG: initialize figlist')))
    if first_figure is None:
        logger.debug(strm(lsafen('empty')))
        return []
    else:
        logger.debug(strm(lsafen(first_figure.figurelist)))
        return first_figure
def text_on_plot(x,y,thistext,coord = 'axes',**kwargs):
    ax = plt.gca()
    if coord == 'axes':
        newkwargs = {'transform':ax.transAxes,'size':'x-large',"horizontalalignment":'center'}
    elif coord == 'data':
        print("Yes, I am using data transform")
        newkwargs = {'transform':ax.transData,'size':'small',"horizontalalignment":'right'}
    color = None
    if 'match_data' in list(kwargs.keys()):
        if isinstance(kwargs['match_data'], list):
            color = kwargs['match_data'][-1].get_color() # get the color of the last line
        elif kwargs['match_data'].get_plot_color() is not None:
            color = kwargs['match_data'].get_plot_color() # don't know when this works, but apparently, it does!
        if color is not None:
            newkwargs.update({'color':color})
        else:
            raise ValueError('You passed match_data to text_on_plot, but I can\'t find a color in the object')
        kwargs.pop('match_data')
    newkwargs.update(kwargs)
    return text(x,y,thistext,**newkwargs)
#{{{subplot_dim
class subplot_dim():
    def __init__(self,firstdim,seconddim):
        self.num = r_[firstdim,seconddim,0]
    def set(self,args,x='',g=True,y='',t='',a=''):
        if isinstance(args, int):
            number = args
            ax = subplot(*tuple(self.num+r_[0,0,number]))
            xlabel(x)
            ylabel(y)
            plt.title(t)
            plt.grid(g)
        elif (isinstance(args, tuple)) and (len(args) == 3):
            # the second value passed is 
            whichsmall = args[2]
            break_into = args[1]
            number = args[0]
            mydims = self.num*r_[1,break_into,1]+r_[
                    0,0,break_into*(number-1)+whichsmall]
            try:
                ax = subplot(*tuple(mydims))
            except:
                print('failed trying subplots: ', mydims)
                raise
            xlabel(x)
            ylabel(y)
            plt.title(t)
            plt.grid(g)
        else:
            print("problem, need to pass either 1 or 3 arguments to set")
            print('type of args: ',type(args))
        return ax
#}}}
def spectrogram(waveform,f_start,f_stop,npoints_fdom=40,tdom_div=2):
    npoints_tdom = waveform.len/tdom_div # this seems to be more legible than above 
    resolution = np.diff(waveform.x[0:2])

    sigma = abs(f_start-f_stop)/np.double(npoints_fdom)
    #print "sigma = %f resolution = %f"%(sigma,resolution)
    if sigma<4*resolution:
        sigma = 4*resolution

    waveform.def_filter(sigma,npoints_tdom)# define the filter and number of points for the spectrogram windowing (define the filter such that the points are spaced sigma apart)

    # go through and apply the filter for some range of points

    f_axis = np.linspace(f_start,f_stop,npoints_fdom)

    specgram = np.zeros((npoints_fdom,npoints_tdom),dtype="complex128")

    for j in range(0,npoints_fdom):

        t_axis, specgram[j,:] = waveform.do_filter(f_axis[j])
        #plot(t_axis,abs(specgram[j,:])) # leave this in for testing what it does in the fdom
    #image(specgram,y=f_axis/1e6,x=t_axis*1e6) # now do an imagehsv (see if we can make imagerybw) plot of the resulting spectrogram
    imshow(abs(specgram),extent=(t_axis[0]*1e6,t_axis[-1]*1e6,f_axis[-1]/1e6,f_axis[0]/1e6)) # now do an imagehsv (see if we can make imagerybw) plot of the resulting spectrogram
    return plt.gca()
def colormap(points,colors,n=256):
    r = np.interp(np.linspace(0,1,n),points,colors[:,0].flatten())
    g = np.interp(np.linspace(0,1,n),points,colors[:,1].flatten())
    b = np.interp(np.linspace(0,1,n),points,colors[:,2].flatten())
    return np.reshape(r_[r,g,b],(3,n)).T

default_cycler = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
