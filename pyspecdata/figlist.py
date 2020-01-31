"""Contains the figure list class

The figure list gives us three things:

*   Automatically handle the display and scaling of nddata units.
*   Refer to plots by name, rather than number (matplotlib has a mechanism for this, which we ignore)
*   A "basename" allowing us to generate multiple sets of plots for different datasets -- *e.g.* 5 plots with 5 names plotted for 3 different datasets and labeled by 3 different basenames to give 15 plots total
*   Ability to run the same code from the command line or from within a python environment inside latex.
    *   this is achieved by choosing figlist (default gui) and figlistl (inherits
        from figlist -- renders to latex -- the :func:`figlist.show` method is
        changed)
    *   potential planned future ability to handle html
*   Ability to handle mayavi plots and matplotlib plots (switch to glumpy, etc.?)
    *   potential planned future ability to handle gnuplot

.. todo:: 

    Currently the "items" that the list tracks correspond to either plot formatting directives (see :func:`figlist.setprops`), text, or figures.

    We should scrap most elements of the current implementation of figlist and rebuild it

    *   currently the figlist is set up to use a context block.  We will not only keep this, but also make it so the individual axes.  Syntax (following a ``fl = figlist_var()`` should look like this: ``with fl['my plot name'] as p:`` and contents of the block would then be ``p.plot(...)``, *etc.*
    *   define an "organization" function of the figlist block.  This allows us
        to use standard matplotlib commands to set up and organize the axes, using
        standard matplotlib commands (twinx, subplot, etc.)
    *   figlist will still have a "next" function, but its purpose will be to simply:
        *   grab the current axis using matplotlib gca() (assuming the id of the axis isn't yet assigned to an existing figlist_axis -- see below)
        *   otherwise, if the name argument to "next" has not yet been called,
            call matplotlib's figure(), followed by subplot(111), then do the
            previous bullet point
        *   the next function is only intended to be called explicitly from within the organization function
    *   figlist will consist simply of a list of figlist_axis objects (a new object type), which have the following attributes:
        *   type -- indicating the type of object:
            *   axis (default)
            *   text (raw latex (or html))
            *   H1 (first level header -- translates to latex section)
            *   H2 (second level...)
        *   the name of the plot
        *   a matplotlib or mayavi axes object
        *   the units associated with the axes
        *   a collection.OrderedDict giving the nddata that are associated with the plot, by name.
            *   If these do not have a name, they will be automatically assigned a name.
            *   The name should be used by the new "plot" method to generate
                the "label" for the legend, and can be subsequently used to quickly
                replace data -- e.g. in a Qt application.
        *   a dictionary giving any arguments to the pyspecdata.core.plot (or countour, waterfall, etc) function
        *   the title -- by default the name of the plot -- can be a setter
        *   the result of the id(...) function, called on the axes object -->
            this can be used to determine if the axes has been used yet
        *   do not use check_units -- the plot method (or contour, waterfall,
            etc.) will only add the nddata objects to the OrderedDict, add the
            arguments to the argument dictionary, then exit
            *   In the event that more than one plot method is called, the name of the underlying nddaata should be changed
        *   a boolean legend_suppress attribute
        *   a boolean legend_internal attribute (to place the legend internally, rather than outside the axis)
        *   a show method that is called by the figlistl show method.  This
            will determine the appropriate units and use them to determine the
            units and scale of the axes, and then go through and call
            pyspecdata.core.plot on each dataset
            (in matplotlib, this should be done with a formatting statement rather than by manipulating the axes themselves)
            and finally call autolegend, unless the legend is supressed
    *   The "plottype" (currently an argument to the plot function) should be an attribute of the axis object
"""
from .plot_funcs import image as this_image

class figlist(object):
    r"""
    Attributes
    ----------
    basename : str
        A basename that can be changed to generate different sets of figures with different basenames.
        For example, this is useful if you are looping over different sets of data,
        and generating the same set of figures for each set of data (which would correspond to a basename).
    figurelist : list
        A list of the figure names
    figdict : dict
        A dictionary containing the figurelist and the figure numbers or objects that they correspond to.
        Keys of this dictionary must be elements of `figurelist`.
    propdict : dict
        Maintains various properties for each element in figurelist.
        Keys of this dictionary must be elements of `figurelist`.
    """
    def __init__(self,*arg,**kwargs):
        r"""Initialize a figure list, which can be used to generate a series of
        figures from the command line or prompt.  Then the same code (if
        `figlist_var` is used) can be included inside a ``python`` environment
        in a latex document.

        Parameters
        ----------
        black : double
            A fractional number giving how "black" "black" is. Typically 1.0 is
            actually too dark and makes things hard to see.
        mlab : object
            If you want to use mayavi, this should be the mlab (module?)
        file_name : str
            This is the argument passed to :func:`self.show`, and used to
            construct the file names.
        """
        self.verbose, self.black, self.env, self.mlab, self.file_name, self.line_spacing = process_kwargs([
            ('verbose',False),
            ('black',0.9),
            ('env',''),
            ('mlab','BLANK'),
            ('file_name','BLANK'),
            ('line_spacing','BLANK'),
            ],
                kwargs, pass_through=True)
        if len(kwargs) > 0:
            self.lplot_kwargs = kwargs
        if self.mlab == 'BLANK': del self.mlab
        if self.file_name == 'BLANK': del self.file_name
        if self.line_spacing == 'BLANK': del self.line_spacing
        if self.verbose: print(lsafe('DEBUG: initialize figlist'))
        if len(arg) == 0:
            self.figurelist = []
        else:
            self.figurelist = arg[0]
        if len(kwargs) > 0:
            self.figurelist.append(kwargs)
        self.units = {}
        self.autolegend_list = {}
        self.twinx_list = {}
        self.basename = None
        return
    def twinx(self,autopad = False,orig = False,color = None):
        #self.figurelist.insert(self.get_fig_number(self.current)-1,{'autopad':False}) #doesn't work because it changes the figure number; I can get the number with fig = gcf(); fig.number, but I can't set it; it would be best to switch to using a list that contains all the figure numbers to match all their names -- or alternatively, note that matplotlib allows you to give them names, though I don't know how that works
        if self.current in list(self.twinx_list.keys()):
            ax1,ax2 = self.twinx_list[self.current]
            if color is not None:
                if 'twinx_color' not in list(self.propdict[self.current].keys()):
                        ax2.tick_params(axis = 'y',colors = color)
                        ax2.yaxis.label.set_color(color)
                        ax2.spines['right'].set_color(color)
                        self.propdict[self.current]['twinx_color'] = color
                else:
                    if color != self.propdict[self.current]['twinx_color']:
                        raise ValueError("conflicting values for the twinx color have been given!!")
        else:
            if autopad: autopad_figure()
            ax1 = gca()
            twinx()
            ax2 = gca()
            self.twinx_list[self.current] = (ax1,ax2)
            if color is not None:
                ax2.tick_params(axis = 'y',colors = color)
                ax2.yaxis.label.set_color(color)
                ax2.spines['right'].set_color(color)
                self.propdict[self.current]['twinx_color'] = color
        if orig:
            sca(ax1)
            return ax1
        else:
            sca(ax2)
            return ax2
    def use_autolegend(self,value = None):
        'No argument sets to true if it\'s not already set'
        if value is None:
            if not self.current in list(self.autolegend_list.keys()):
                self.autolegend_list.update({self.current:True})
            else: #leave it alone
                return
        else: #passed an explicit value
            self.autolegend_list.update({self.current:value})
            return
    def push_marker(self):
        if not hasattr(self,'pushlist'):
            self.pushlist = []
        if not hasattr(self,'pushbasenamelist'):
            self.pushbasenamelist = []
        self.pushlist.append(self.current)
        self.pushbasenamelist.append(self.basename)
        return
    def pop_marker(self):
        self.basename = self.pushbasenamelist.pop()
        self.next(self.pushlist.pop())
        return
    def get_fig_number(self,name):
        cleanlist = [x for x in self.figurelist if isinstance(x, str)]
        try:
            return cleanlist.index(name)+1
        except ValueError:
            raise ValueError(strm("You are looking for",name,
                "which isn't in the list of figures",cleanlist))
    def next(self,input_name, legend=False, boundaries=None, twinx=None, **kwargs):
        r"""Switch to the figure given by input_name, which is used not only as
        a string-based name for the figure, but also as a default title and as
        a base name for resulting figure files.

        **In the future, we actually want this to track the appropriate axis object!**

        Parameters
        ----------
        legend : bool
            If this is set, a legend is created *outside* the figure.
        twinx : {0,1}
            :1: plots on an overlayed axis (the matplotlib twinx) whose y axis
                is labeled on the right when you set this for the first time, you
                can also set a `color` kwarg that controls the coloring of the
                right axis. 
            :0: used to switch back to the left (default) axis
        boundaries :
            **need to add description**
        kwargs : dict
            Any other keyword arguments are passed to the matplotlib (mayavi)
            figure() function that's used to switch (create) figures.
        """
        if not hasattr(self,'figdict'):
            self.figdict = {} # the dictionary of the various figures
        if not hasattr(self,'propdict'):
            self.propdict = {} # the properties belonging to those same figures
        logger.debug(strm("for plot",input_name,"basename is",self.basename))
        if (self.basename is not None #basename for groups of figures
                # I need to check that the basename hasn't already been added
                and not input_name.startswith(self.basename)):
            name = self.basename + ' ' + input_name
        else:
            logger.debug(strm("not using a basename",self.basename is not None))
            name = input_name
        if name.find('/') > 0:
            raise ValueError("don't include slashes in the figure name, that's just too confusing")
        logger.debug(strm('called with',name))
        if name in self.figurelist:# figure already exists
            if hasattr(self,'mlab'):
                # with this commit, I removed the kwargs and bgcolor, not sure why
                fig = self.mlab.figure(self.get_fig_number(name))
                fig.scene.render_window.aa_frames = 20
                fig.scene.anti_aliasing_frames = 20
            else:
                fig = figure(self.get_fig_number(name))
            self.current = name
            if self.verbose: print(lsafen('in',self.figurelist,'at figure',self.get_fig_number(name),'switched figures'))
            if boundaries is not None:
                if 'boundaries' not in list(self.propdict[self.current].keys()) or self.propdict[self.current]['boundaries'] != boundaries:
                    raise ValueError("You're giving conflicting values for boundaries")
            if legend:
                if 'legend' not in list(self.propdict[self.current].keys()) or self.propdict[self.current]['legend'] != legend:
                    raise ValueError("You're giving conflicting values for legend")
        else:# figure doesn't exist yet
            if hasattr(self,'current'):
                last_figure_number = self.get_fig_number(self.current)
            else:
                last_figure_number = 0
            self.current = name
            if self.current not in list(self.propdict.keys()):
                self.propdict[self.current] = {}
            if boundaries == False:
                self.propdict[self.current]['boundaries'] = False
                self.setprops(boundaries = False)
            if legend:
                self.propdict[self.current]['legend'] = True
                if 'figsize' not in list(kwargs.keys()):
                    kwargs.update({'figsize':(12,6)})
                if hasattr(self,'mlab'):
                    fig = self.mlab.figure(last_figure_number+1,bgcolor = (1,1,1),**kwargs)
                    fig.scene.render_window.aa_frames = 20
                    fig.scene.anti_aliasing_frames = 20
                else:
                    fig = figure(last_figure_number+1,**kwargs)
                fig.add_axes([0.075,0.2,0.6,0.7]) # l b w h
                self.use_autolegend('outside')
            else:
                self.propdict[self.current]['legend'] = False
                fig = figure(last_figure_number+1,**kwargs)
                if hasattr(self,'mlab'):
                    fig = self.mlab.figure(last_figure_number+1,bgcolor = (1,1,1),**kwargs)
                    fig.scene.render_window.aa_frames = 20
                    fig.scene.anti_aliasing_frames = 20
                else:
                    fig = figure(last_figure_number+1,**kwargs)
                if twinx is not None:
                    fig.add_subplot(111)
            logger.debug(strm('added, figure',len(self.figurelist)+1,'because not in figurelist',self.figurelist))
            self.figurelist.append(name)
            self.figdict.update({self.current:fig})
            if boundaries == False:
                self.setprops(boundaries = True)# set this back
        if twinx is not None:
            self.propdict[self.current]['twinx'] = True
            if twinx == 0:
                self.twinx(orig = True)
                fig = gcf()
            elif twinx == 1:
                self.twinx()
                fig = gcf()
            else:
                raise ValueError('If you pass twinx, pass 0 for the original or 1 for the right side')
            self.figdict.update({self.current:fig})
        return fig
    def plot(self,*args,**kwargs):
        r"""
        Parameters
        ----------
        linestyle: {':','--','.','etc.'}
            the style of the line
        plottype: {'semilogy','semilogx','loglog'}
            Select a logarithmic plotting style.
        nosemilog: True
            Typically, if you supply a log-spaced axis,
            a semilogx plot will be automatically selected.
            This overrides that behavior.
            Defaults to False.
        """
        if 'label' in list(kwargs.keys()):
            self.use_autolegend()
        human_units = True
        if 'human_units' in list(kwargs.keys()):
            human_units = kwargs.pop('human_units')
        if human_units:
            firstarg = self.check_units(args[0],0,1) # check units, and if need be convert to human units, where x is the first dimension and y is the last
        else:
            firstarg = args[0]
        if 'label' not in list(kwargs.keys()) and isinstance(args[0],nddata):
            thisname = args[0].name()
            if thisname is not None:
                kwargs['label'] = thisname
        retval = plot(*tuple((firstarg,)+args[1:]),**kwargs)#just a placeholder for now, will later keep units + such
        ax = gca()
        if ax.get_title() is None or len(ax.get_title()) == 0:
            try:
                title(self.current)
            except:
                title('untitled')
        return retval
    def phaseplot_finalize(self):
        ("Performs plot decorations that are typically desired for a manual phasing"
        " plot.  This assumes that the ``y``-axis is given in units of half-cycles"
        " ($\pi$ radians).")
        ax = gca()
        ylim(-1,1)
        gridandtick(ax)
        ylabel(r'$\phi / \pi$')
        # now show the pi/2 lines
        axhline(y = 0.5,color = 'r',alpha = 0.5,linewidth = 2)
        axhline(y = -0.5,color = 'r',alpha = 0.5,linewidth = 2)
        return
    def check_units(self, testdata, x_index, y_index,
            verbose=False):
        if verbose: print("-"*30)
        if verbose: print("called check_units for figure",self.current)
        if isinstance(testdata,nddata):
            if verbose: print("(check_units) it's nddata")
            testdata = testdata.copy().human_units()
            if len(testdata.dimlabels) > 1:
                if verbose: print("(check_units) more than one dimension")
                if not hasattr(self,'current'):
                    raise ValueError("give your plot a name (using .next()) first! (this is used for naming the PDF's etc)")
                if self.current in list(self.units.keys()):
                        theseunits = (testdata.get_units(testdata.dimlabels[x_index]),testdata.get_units(testdata.dimlabels[y_index]))
                        if theseunits != self.units[self.current] and theseunits[0] != self.units[self.current]:
                                raise ValueError("the units don't match (old units %s and new units %s)! Figure out a way to deal with this!"%(theseunits,self.units[self.current]))
                else:
                    if isinstance(testdata,nddata):
                        self.units[self.current] = (testdata.get_units(testdata.dimlabels[x_index]),testdata.get_units(testdata.dimlabels[y_index]))
            else:
                if verbose: print("(check_units) only one dimension")
                if not hasattr(self,'current'):
                    self.next('default')
                if self.current in list(self.units.keys()):
                    theseunits = (testdata.get_units(testdata.dimlabels[x_index]))
                    testunits = self.units[self.current]
                    if theseunits != testunits:
                        if isinstance(testunits, tuple) and testunits[1] is None:
                            pass
                        else:
                            raise ValueError("the units don't match (old units %s and new units %s)! Figure out a way to deal with this!"%(theseunits,self.units[self.current]))
                else:
                    self.units[self.current] = (testdata.get_units(testdata.dimlabels[x_index]))
        if verbose: print("-"*30)
        return testdata
    def adjust_spines(self,spines):
        ax = gca()
        #{{{ taken from matplotlib examples
        for loc, spine in list(ax.spines.items()):
            if loc in spines:
                spine.set_position(('outward',10)) # outward by 10 points
                spine.set_smart_bounds(True)
            else:
                spine.set_color('none') # don't draw spine

        # turn off ticks where there is no spine
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            # no yaxis ticks
            ax.yaxis.set_ticks([])

        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            # no xaxis ticks
            ax.xaxis.set_ticks([])
        #}}}
    def grid(self):
        ax = gca()
        if self.black:
            gridandtick(ax,gridcolor = r_[0.5,0.5,0.5])
        else:
            gridandtick(ax,gridcolor = r_[0,0,0])
        return
    image = this_image.fl_image
    def marked_text(self,marker,input_text="",sep='\n'):
        """Creates a named `marker` where we can place text.   If `marker`
        has been used, goes back and places text there."""
        if not hasattr(self,'textdict'):
            self.textdict = {}
        if marker in list(self.textdict.keys()):
            idx = self.textdict[marker]
            self.figurelist[idx]['print_string'] = (
                    self.figurelist[idx]['print_string']
                    + sep + input_text )
        else:
            self.setprops(print_string=input_text)
            idx = len(self.figurelist)-1
            self.textdict[marker] = idx
    def text(self,mytext):
        self.setprops(print_string = mytext)
    def setprops(self,**kwargs):
        self.figurelist.append(kwargs)
    def show_prep(self):
        for k,v in list(self.autolegend_list.items()):
            kwargs = {}
            if v:
                if isinstance(v, str):
                    if v[0:7] == 'colored':
                        kwargs.update(dict(match_colors = True))
                        v = v[7:]
                        if v == '':
                            v = True
                    if v == 'outside':
                        kwargs.update(dict(bbox_to_anchor=(1.05,1),loc = 2,borderaxespad=0.))
                self.next(k)
                logger.debug(strm("I am about to assign a legend for ",k,". Is it in the figurelist?:",k in self.figurelist))
                logger.debug(strm("print out the legend object:",gca().legend()))
                try:
                    autolegend(**kwargs)
                except:
                    try:
                        self.twinx(orig = True)
                    except Exception as e:
                        raise Exception(strm('error while trying to run twinx to place legend for',k,'\n\tfiglist is',self.figurelist,explain_error(e)))
                    try:
                        autolegend(**kwargs)
                    except Exception as e:
                        raise Exception(strm('error while trying to run autolegend function for',k,'\n\tfiglist is',self.figurelist,explain_error(e)))
    def show(self,*args,**kwargs):
        self.basename = None # must be turned off, so it can cycle through lists, etc, on its own
        if 'line_spacing' in list(kwargs.keys()): kwargs.pop('line_spacing')# for latex only
        if len(kwargs) > 0:
            raise ValueError("didn't understand kwargs "+repr(kwargs))
        logger.debug(strm("before show_prep, figlist is",self.figurelist))
        logger.debug(strm("before show_prep, autolegend list is",self.autolegend_list))
        self.show_prep()
        #{{{ just copy from fornnotebook to get the print string functionality
        kwargs = {}
        for figname in self.figurelist:
            logger.debug(strm("showing figure"+lsafen(figname)))
            if isinstance(figname, dict):
                kwargs.update(figname)
                if 'print_string' in kwargs:
                    print('\n\n')
                    print(kwargs.pop('print_string'))
                    print('\n\n')
        #}}}
        if len(args) == 1:
            if (args[0][:-4] == '.pdf') or (args[0][:-4] == '.png') or (args[0][:-4] == '.jpg'):
                print("you passed me a filename, but I'm just burning it")
        if hasattr(self,'mlab'):
            print("running mlab show!")
            self.mlab.show()
        else:
            #print "not running mlab show!"
            show()
    def label_point(self, data, axis, value, thislabel,
            show_point=True, xscale=1, **new_kwargs):
        """only works for 1D data: assume you've passed a single-point nddata, and label it

        xscale gives the unit scaling

        ..todo::

            Improve the unit scaling, so that this would also work.

            Allow it to include a format string that would use the value.
        Parameters
        ----------

        show_point : bool

            Defaults to `True`. Actually generate a point (circle), *vs.*
            just the label.
        """
        kwargs = {'alpha':0.5,'color':'k','ha':'left','va':'bottom','rotation':45,'size':14}
        kwargs.update(new_kwargs)
        y = double(data[axis:value].data)
        x_ind = argmin(abs(data.getaxis(axis)-value))
        x = data.getaxis(axis)[x_ind]
        text(x/xscale, y, thislabel, **kwargs)
        if show_point:
            plot(x/xscale, y, 'o', color=kwargs["color"],
                    alpha=kwargs["alpha"])
        return
    def header(self,number_above,input_string):
        header_list = ['\\section','\\subsection','\\subsubsection','\\paragraph','\\subparagraph']
        self.text(header_list[number_above+1]+'{%s}'%input_string)
        return number_above + 1
    def mesh(self,plotdata,Z_normalization = None,equal_scale = True,
            lensoffset = 1e-3,
            show_contours = False,
            grey_surf = False,
            **kwargs):
        plotdata = self.check_units(plotdata,0,1)
        if hasattr(self,'mlab'):
            fig = self.figdict[self.current]
            fig.scene.disable_render = True
            X,Y,Z,x_axis,y_axis = plotdata.matrices_3d(also1d = True)# return the axes, and also alter "plotdata" so it's downsampled
            X_normalization = X.max()
            X /= X_normalization
            if equal_scale:
                Y_normalization = X_normalization
            else:
                Y_normalization = Y.max()
            Y /= Y_normalization
            if Z_normalization is None:
                Z_normalization = Z.flatten().max()
            Z /= Z_normalization
            surf_kwargs = {}
            if grey_surf:
                surf_kwargs.update(color = (0.5,0.5,0.5))# opacity and the contour lines don't play well, otherwise I would like to make this transluscent
            self.mlab.surf(X,Y,Z,**surf_kwargs)
            if show_contours:
                contour_kwargs = {'line_width':24}
                contour_kwargs.update(opacity = 0.5)
                if not grey_surf:
                    contour_kwargs.update(color = (1,1,1))
                self.mlab.contour_surf(X,Y,Z+lensoffset,contours = r_[-1:1:10j].tolist(),**contour_kwargs)
                contour_kwargs.update(opacity = 0.1)
                self.mlab.contour_surf(X,Y,Z+lensoffset,contours = r_[-1:1:46j].tolist(),**contour_kwargs)# for some reason, 46 gives alignment (I think 9+1 and 9*5+1)
            if equal_scale:
                self.generate_ticks(plotdata,(x_axis,y_axis),X_normalization,Z_normalization,verbose = True)
            else:
                self.generate_ticks(plotdata,(x_axis,y_axis),X_normalization,Z_normalization,y_rescale = Y_normalization/X_normalization,verbose = True)
            fig.scene.disable_render = False
        else:
            # this should be upgraded, or rather moved to here
            plotdata.meshplot(alpha=1.0, cmap=cm.jet, **kwargs)
        return Z_normalization
    def generate_ticks(self,plotdata,axes,rescale,z_norm = None,y_rescale = 1,text_scale = 0.05,verbose = False,follow_surface = False,
            lensoffset = 0.5e-2,
            line_width = 1e-3,
            tube_radius = 1e-3,
            fine_grid = False,
            ):
        'generate 3d ticks and grid for mayavi'
        if follow_surface and z_norm is None:
            raise ValueError("if you choose to generate the mesh -- i.e. follow the surface -- then you need to pass the z normalization")
        x_axis,y_axis = axes
        x_dim = plotdata.dimlabels[0]
        y_dim = plotdata.dimlabels[1]
        def gen_list(thisaxis,desired_ticks = 7.):
            #{{{ out of the following list, choose the one that gives as close as possible to the desired ticks
            axis_span = thisaxis.max() - thisaxis.min()
            possible_iterators = r_[0.1,0.5,1,5,10,20,30,50,100,200,500,1000]
            iterator = possible_iterators[argmin(abs(axis_span/desired_ticks -
                possible_iterators))]
            #}}}
            if verbose: print('iterator is',iterator)
            return iterator,r_[ceil(thisaxis.min()/iterator):
                floor(thisaxis.max()/iterator)+1]*iterator
        #{{{ now, I need to get the list of multiples that falls inside the axis span
        xiterator,xlist = gen_list(x_axis)
        yiterator,ylist = gen_list(y_axis)
        if verbose: print('range of x ',x_axis.min(),x_axis.max())
        if verbose: print('xlist',xlist)
        if verbose: print(plotdata.unitify_axis(0))
        if verbose: print('range of y ',y_axis.min(),y_axis.max())
        if verbose: print('ylist',ylist)
        if verbose: print(plotdata.unitify_axis(1))
        #}}}
        if xiterator < 1:
            x_ticklabels = ['{:0.1f}'.format(j) for j in xlist]
        else:
            x_ticklabels = ['{:0.0f}'.format(j) for j in xlist]
        if yiterator < 1:
            y_ticklabels = ['{:0.1f}'.format(j) for j in ylist]
        else:
            y_ticklabels = ['{:0.0f}'.format(j) for j in ylist]
        #{{{ rescale absolutely everything
        xlist /= rescale
        ylist /= (rescale*y_rescale)
        x_axis /= rescale
        y_axis /= (rescale*y_rescale)
        #}}}
        x_range = r_[x_axis.min(),x_axis.max()]
        y_range = r_[y_axis.min(),y_axis.max()]
        extension_factor = text_scale * 3
        #{{{ y ticks
        if follow_surface:
            if fine_grid:
                dy = ylist[1]-ylist[0]
                finer_ylist = r_[ylist[0]-dy:ylist[-1]+dy:1j*((len(ylist)+2-1)*5+1)]
                finer_ylist = finer_ylist[finer_ylist>=y_axis.min()]
                finer_ylist = finer_ylist[finer_ylist<=y_axis.max()]
            else:
                finer_ylist = ylist
            for j,y in enumerate(finer_ylist):
                x_linedata = plotdata.getaxis(x_dim)/rescale
                z_linedata = plotdata[y_dim:(y*rescale)].data.flatten()/z_norm
                self.mlab.plot3d(x_linedata,y*ones_like(x_linedata),
                        z_linedata+lensoffset,
                        color = (0,0,0), line_width = line_width,
                        tube_radius = tube_radius)
        for j,y in enumerate(ylist):
            self.mlab.plot3d(x_range+extension_factor*r_[-1,1],
                    y*ones(2),zeros(2),
                    color = (0,0,0), line_width = line_width,
                    tube_radius = tube_radius)
            self.mlab.text3d(x_range[0]-2*extension_factor, y, 0,
                    y_ticklabels[j],color = (0,0,0),
                    scale = text_scale # in figure units
                    )
            self.mlab.text3d(x_range[1]+2*extension_factor, y, 0,
                    y_ticklabels[j],color = (0,0,0),
                    scale = text_scale # in figure units
                    )
        self.mlab.text3d(x_range[1] + 3 * extension_factor,y_range.mean(), 0,
                plotdata.unitify_axis(1), color = (0,0,0),
                scale = text_scale,
                orient_to_camera = False,
                orientation = (0,0,90))# the last angle appears to be rotaiton about z
        #}}}
        #{{{ x ticks
        if follow_surface:
            if fine_grid:
                dx = xlist[1]-xlist[0]
                finer_xlist = r_[xlist[0]-dx:xlist[-1]+dx:1j*((len(xlist)+2-1)*5+1)]
                finer_xlist = finer_xlist[finer_xlist>=x_axis.min()]
                finer_xlist = finer_xlist[finer_xlist<=x_axis.max()]
            else:
                finer_xlist = xlist
            for j,x in enumerate(finer_xlist):
                y_linedata = plotdata.getaxis(y_dim)/(rescale*y_rescale)
                z_linedata = plotdata[x_dim:(x*rescale)].data.flatten()/z_norm
                self.mlab.plot3d(x*ones_like(y_linedata),y_linedata,
                        z_linedata+lensoffset,
                        color = (0,0,0), line_width = line_width,
                        tube_radius = tube_radius)
        for j,x in enumerate(xlist):
            self.mlab.plot3d(x*ones(2),y_range+extension_factor*r_[-1,1],
                    zeros(2),
                    color = (0,0,0), line_width = line_width,
                    tube_radius = tube_radius)
            self.mlab.text3d(x, y_range[0]-2*extension_factor, 0,
                    x_ticklabels[j],color = (0,0,0),
                    scale = text_scale # in figure units
                    )
            self.mlab.text3d(x, y_range[1]+2*extension_factor, 0,
                    x_ticklabels[j],color = (0,0,0),
                    scale = text_scale # in figure units
                    )
        self.mlab.text3d(x_range.mean(), y_range[1] + 3 * extension_factor,
                0,
                plotdata.unitify_axis(0), color = (0,0,0),
                scale = text_scale,
                orient_to_camera = False,
                orientation = (0,0,180))# the last angle appears to be rotaiton about z
        #}}}
        return
    def __enter__(self):
        return self
    def __exit__(self, exception_type, exception_value, traceback):
        r'''show the plots, unless there are errors.

        Because this is executed before raising any errors, we want to avoid showing any plots if there are errors.
        Otherwise, it gets very confusing.
        '''
        if exception_type is None:
            if hasattr(self,'file_name'):
                if hasattr(self,'line_spacing'):
                    self.show(self.file_name,line_spacing = self.line_spacing)
                else:
                    self.show(self.file_name)
            else:
                self.show()
            return
