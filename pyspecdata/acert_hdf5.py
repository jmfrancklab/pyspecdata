from .core import *
from .load_files import *
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from .datadir import get_notebook_dir
from nmr import phaseopt
from sympy import var
import tables
import h5py
import warnings
import re
import os
import scipy.io

def gen_composite(basenames,collected_data,
    timeconst = 4e-9/1.8, 
    plot_title = ''):
    """process several ringdown traces acquired with "stepped" Rx amplitude and generate the composite trace of the ringdown with high dynamic range

    :param basenames: list of collected data, in order from highest power (shortest deadtime) to lowest (longest deadtime)
    :param collected_data: collected data in dictionary format
    :param timeconst: time constant of the receiver opening -- used for the correction of the initial time points
    :param plot_title: base name for the plot title
    :return: (composite trace, the list of amplifications for each input trace, the overall max -- in input units)"""
    startpoint = collected_data[basenames[0]].get_prop('execution_deadtime') # startpoint of the receiver opening -- used for the correction of the initial time points
    amplification_list = {}
    for this_basename in basenames:
        try:
            data = collected_data[this_basename]
        except:
            raise ValueError("I can't find key "+repr(this_basename)+" in "+repr(collected_data.keys()))
        attenuator_opened = False
        if data.name().find('open') > -1:
            attenuator_opened = True
        if attenuator_opened:
            prev_name = data.name().replace('_open','') # name of experiment before opening attenuator
            amplification = data['phcyc1',-1].runcopy(lambda x: abs(x)**2).sum('t2').data / collected_data[prev_name]['phcyc1',-1].runcopy(lambda x: abs(x)**2).sum('t2').data
            amplification = sqrt(amplification)
            for ampkey in amplification_list.keys():
                amplification_list[ampkey] *= amplification # because we assume this one is lower in power than all previous that were entered, all previous are this much higher in amplitude than we thought they were
        amplification_list[data.name()] = 1.0 # assume for now that we've hit the "most open" setting
    saturation_point = array([abs(v['phcyc1',-1]).run(max,'t2').data for v in collected_data.values()])
    saturation_point = saturation_point.max() # assume that I've run the experiments such that the highest peak out of all the experiments comes right up to the saturation point
    print 'saturation point:',saturation_point
    #{{{ generate the composite
    thisname = basenames[0]
    fl.next(plot_title + ' separate',legend = True)
    composite_ringdown = abs(collected_data[thisname]['t2':(None,125e-9),'phcyc1',-1]) * amplification_list[thisname] / saturation_point
    thisname = basenames[-1]
    temp = abs(collected_data[thisname]['t2':(None,125e-9),'phcyc1',-1]) * amplification_list[thisname] / saturation_point
    composite_ringdown['t2':(43e-9,None)] = temp['t2':(43e-9,None)]
    #{{{ adjust by the rise time
    timeslice = (startpoint+1.5e-9,startpoint + timeconst*5)
    composite_ringdown['t2':timeslice] /=\
            composite_ringdown['t2':timeslice].\
            fromaxis('t2',lambda x: 1.0 - exp(-(x-startpoint)/timeconst))
    composite_ringdown = composite_ringdown['t2':(timeslice[0],None)]
    lines = fl.plot(composite_ringdown,'k',alpha = 0.5,linewidth = 3)#,alpha = 0.75)
    #}}}
    for thisname in basenames:
        fl.next(plot_title + ' separate')
        attenuator_opened = False
        if thisname.find('open') > -1:
            attenuator_opened = True
        thisdata = collected_data[thisname].copy()
        linestyle = '-'
        if not attenuator_opened:
            linestyle = '--'
        fl.next(plot_title + ' $\Delta c = 0$',legend = True)
        fl.plot(abs(thisdata['t2':(None,125e-9),'phcyc1',0]),alpha = 0.5)
        fl.next(plot_title + ' $\Delta c = +1$',legend = True)
        fl.plot(abs(thisdata['t2':(None,125e-9),'phcyc1',1]),alpha = 0.5)
        fl.next(plot_title + ' $\Delta c = \pm 2$',legend = True)
        fl.plot(abs(thisdata['t2':(None,125e-9),'phcyc1',2]),alpha = 0.5)
        fl.next(plot_title + ' separate')
        thisdata = abs(thisdata['t2':(None,125e-9),'phcyc1',-1]) * amplification_list[thisname] / saturation_point
        if attenuator_opened:
            fl.plot(thisdata,linestyle,color = prevcolor,alpha = 0.75)
            axvline(x = thisdata.get_prop('execution_deadtime')/1e-9,color = lines[-1].get_color(),alpha = 0.3,linewidth = 3)
        else:
            lines = fl.plot(thisdata,linestyle,alpha = 0.75)
            prevcolor = lines[-1].get_color()
        fl.next(plot_title + ' semilog',legend = True)
        title('semilog')
        thisdata = thisdata['t2':(thisdata.get_prop('execution_deadtime'),125e-9)] # stop at 125, just so I'm not skewed by the baseline
        plot(thisdata.getaxis('t2')/1e-9,log10(thisdata.data),'.',label = thisname)
    return composite_ringdown,amplification_list,saturation_point
def plot_comparison(input_data,
        date, fl, phase=True,
        discard_error=True,
        scale=1.0):
    "Compare a series of cw experiments"
    # {{{ load initial info
    list_of_cw_files, mod_amp = zip(*input_data)
    single_date = True
    if date is None:
        mod_amp,date = zip(*mod_amp)
        single_date = False
    # }}}
    # {{{ compare the noise levels
    for complex_type in ['real','imag']:
        fl.next('compare noise level %s'%complex_type,
                legend=True)
    ax = gca()
    color_cycle = ax._get_lines.color_cycle
    normalization = 1e-7 # a standard number
    for j in range(0,len(list_of_cw_files)):# I end up looping back over and reloading because they might have different axes, though in an older format, I did this all in one batch
        next_color = next(color_cycle)
        short_basename = list_of_cw_files[j]
        if single_date:
            file_regexp = ('%s.*'%date + short_basename + '\\.')
        else:
            print "note single date -- searching for"
            file_regexp = ('%s.*'%date[j] + short_basename + '\\.')
        data = find_file(file_regexp, exp_type='cw', postproc=acert.postproc_blank)['repeats',0]
        for complex_type in ['real','imag']:
            fl.next('compare noise level %s'%complex_type)
            if complex_type == 'real':
                plotdata = data.runcopy(real)
            else:
                plotdata = data.runcopy(imag)
            fl.plot(plotdata/normalization - 2*j,
                    alpha=0.5, color=next_color,
                    label=short_basename)
    # }}}
    short_basename = list_of_cw_files[0]
    collected_data = []
    for j in range(0,len(list_of_cw_files)):
        short_basename = list_of_cw_files[j]
        if single_date:
            data = find_file('%s.*'%date + short_basename + '\\.', exp_type='cw', phase=phase)
        else:
            data = find_file('%s.*'%date[j] + short_basename + '\\.', exp_type='cw', phase=phase)
        if discard_error: data.set_error(None)
        data /= mod_amp[j]
        if j == 0:
            data_shape = ndshape(data)
            data_shape.add_correctly((len(list_of_cw_files),'indirect'))
        print ndshape(data)
        collected_data.append(data)
    fl.next('cw plots -- real + imag',legend = True,boundaries = False)
    normalization = max(r_[array(
            map(lambda x:
                x.runcopy(real).run(abs).data.max(),
                    collected_data)).mean(),normalization])
    normalization /= scale
    ax = gca()
    color_cycle = ax._get_lines.color_cycle
    for j in range(0,len(list_of_cw_files)):# I end up looping back over and reloading because they might have different axes, though in an older format, I did this all in one batch
        short_basename = list_of_cw_files[j]
        if single_date:
            data = find_file('%s.*'%date + short_basename + '\\.', exp_type='cw', phase=phase)
        else:
            data = find_file('%s.*'%date[j] + short_basename + '\\.', exp_type='cw', phase=phase)
        if discard_error: data.set_error(None)
        data /= mod_amp[j]
        next_color = next(color_cycle)
        fl.plot(data.runcopy(real)/normalization - 2*j,
                alpha = 0.75,color = next_color,label = list_of_cw_files[j])
        autolegend()
        fl.plot(data.runcopy(imag)/normalization - 2*j,
                alpha = 0.25,color = next_color)
def plot_coherence_diagram(ax1,list_of_DeltaC,list_of_sizes,plotidx,outof):
    #print "new coherence diagram"
    aliases_to_search = [0,-1,1,-2,2]
    evolution_width = 6
    #yscale = 1.0/(2*2.1) # everything below goes up to 2.1 -- so now it spans -0.5,0.5
    yscale = 1.0/(2*2.5) # everything below goes up to 2.1 -- so now it spans -0.5,0.5, but I want to scale down further so they are separated
    gridpos = -0.5 - plotidx # now each spans 0 to -1
    time_points = [0,1] # start with line from zero to first pulse
    for j in range(len(list_of_DeltaC)):
        time_points.append(1+(j+1)*evolution_width) # add a time point for evolution after each pulse
    time_points = array(time_points)
    coherence_levels = [0] # always start at zero
    possible_pathways = [[0]]
    #print 'list of delta c',list_of_DeltaC
    for j,thisDeltaC in enumerate(list_of_DeltaC):
        new_possible_pathways = []
        for pathway in possible_pathways:
            for thisalias in aliases_to_search:
                testDeltaC = thisDeltaC + thisalias*list_of_sizes[j]
                newcohlevel = pathway[-1]+testDeltaC
                if not abs(newcohlevel)>1:
                    new_possible_pathways.append(pathway + [newcohlevel])
        possible_pathways = new_possible_pathways
    selected_a_pathway = False
    #print 'possible pathways',possible_pathways
    for pathway in possible_pathways:
        if pathway[-1] == -1:
            show_coherence_pathway(ax1,time_points,pathway,yscale,gridpos,linewidth = 2)
            selected_a_pathway = True
    #print 'selected a pathway?',selected_a_pathway
    if not selected_a_pathway:
        # show one impossible pathway
        fallback_pathway = [0]
        for j,thisDeltaC in enumerate(list_of_DeltaC):
            #{{{ just select the one that's closest to zero, for asthetic purposes
            selectedcohlevel = 100 # a random big number
            for thisalias in aliases_to_search:
                testDeltaC = thisDeltaC + thisalias*list_of_sizes[j]
                newcohlevel = fallback_pathway[-1]+testDeltaC
                if abs(newcohlevel) < abs(selectedcohlevel):
                    selectedcohlevel = newcohlevel
            #}}}
            fallback_pathway.append(selectedcohlevel)
        show_coherence_pathway(ax1,time_points,fallback_pathway,yscale,gridpos)
    return
def show_coherence_pathway(ax1, time_points, coherence_levels,
        yscale, gridpos, linewidth = 1):
    pulse_width = 0.3
    coherence_levels = array(coherence_levels)
    coherence_levels = r_[coherence_levels,coherence_levels[-1]]
    #print time_points,coherence_levels
    #print len(time_points),len(coherence_levels)
    x = (time_points.reshape(-1,1)*r_[1,1]+r_[0,pulse_width]).flatten() # now there is a pulse width for each
    x = r_[x[1:-1]] # coherence level begins only after the pulse width, and there should be only one zero point
    y = (coherence_levels.reshape(-1,1)*r_[1,1]).flatten() # duplicate all the time points
    y = y[:-2]
    #print len(x),len(y)
    #print x,y
    #{{{ make mask for red line(invalid pathways)
    mask = abs(y)>1
    if y[-1] != -1:
        mask[-1] = True
    mask[:-1] = logical_or(mask[:-1],mask[1:]) # alwys include the preceding point
    redlines = []
    for j in range(len(mask)):
        if mask[j] and mask[j-1]:
            if j == 0:
                raise IndexError("You're trying to continue a"
                        " line when no line has been started!")
            try:
                redlines[-1] += [(x[j],y[j]*yscale+gridpos)]
            except IndexError, e:
                raise IndexError(repr(e)+"\n"+("the length of x is"
                    " {:d}, of y is {:d},"
                    " and redlines is {:s}, trying to index"
                    " {:d}").format(
                        len(x),len(y),repr(redlines),j))
        elif mask[j]:# start a new line
            redlines += [[(x[j],y[j]*yscale+gridpos)]]
    #print 'redlines:',redlines
    #print 'len(redlines):',len(redlines)
    #print 'gridpos array:',[(0,gridpos)]*len(redlines)
    #}}}
    lines = LineCollection([zip(x,           y*yscale        ),
                            zip(x[r_[0,-1]], r_[1,1]*yscale  ),     
                            zip(x[r_[0,-1]], r_[0,0]*yscale  ),
                            zip(x[r_[0,-1]], r_[-1,-1]*yscale),
                            ],
                            lw = [linewidth,0.5,0.5,0.5], alpha = 0.7, linestyles = ['-',':',':',':'],
                            color = 'k',
                            offsets = [(0,gridpos),(0,gridpos),(0,gridpos),(0,gridpos)])
    ax1.text(x[0]-0.1,0*yscale+gridpos,'0',ha = 'right',va = 'center',size = 'xx-small')
    ax1.text(x[0]-0.1,1*yscale+gridpos,'+1',ha = 'right',va = 'center',size = 'xx-small')
    ax1.text(x[0]-0.1,-1*yscale+gridpos,'-1',ha = 'right',va = 'center',size = 'xx-small')
    ax1.add_collection(lines)
    if len(redlines) > 0:
        lines = LineCollection(redlines,
                                lw = [1]*len(redlines), alpha = 1.0,
                                color = 'r')
                                #offsets = [(0,gridpos)]*len(redlines))
        ax1.add_collection(lines)
    ##{{{ going from .1 above and below out to 2.1, add the red lines
    #ax1.add_patch(Rectangle([0,gridpos - 2.1*yscale],x[-1],1*yscale,ec = 'none',fc = 'r',alpha = 0.25))#xy, width, height)
    #ax1.add_patch(Rectangle([0,gridpos + 1.1*yscale],x[-1],1*yscale,ec = 'none',fc = 'r',alpha = 0.25))#xy, width, height)
    #ax1.add_patch(Rectangle([x[-2],gridpos-0.1*yscale],x[-1]-x[-2],1.2*yscale,ec = 'none',fc = 'r',alpha = 0.25))#xy, width, height)
    ##}}}
    return
def show_pathways(ax1,plotdata,
        selected_pulses = None, # if given, this is a dictionary like {'phcyc3':(-1,2)}, where the second element of the tuple is the number of phase cycle steps selected
        verbose = False
        ):
    ("Plots the coherence pathway diagrams associated with `plotdata`."
            " It assumes that the phase-cycling dimensions of `plotdata` are"
            " named like ``phcyc1``, ``phcyc2``, *etc.*")
    if verbose: print "(show_pathways): shape of input data",ndshape(plotdata)
    pulseaxes = [j for j in plotdata.dimlabels if j[:5]=='phcyc']
    if verbose: print "(show_pathways): determined the following to be pulse axes",pulseaxes
    if len(pulseaxes) == 0 and selected_pulses is None:
        raise ValueError("I can't find any phase cycling dimensions, so I don't"
                " know what you want me to do!!")
    netshape = ones(len(plotdata.data.shape),dtype = int)
    #{{{ now, figure out the shape of the matrix where I knock out everything but the phcyc dims -- seems there should be an nddata-like way to do this
    for j in pulseaxes:
        thisindex = plotdata.dimlabels.index(j)
        netshape[thisindex] = plotdata.data.shape[thisindex]
    if verbose: print "(show_pathways): and the shapes of those axes are",netshape
    #}}}
    if selected_pulses is not None:
        pulseaxes += selected_pulses.keys()
    coherence_changes = ones((netshape.prod(),len(pulseaxes)),dtype = int)
    pulseaxes_sorted = sort(pulseaxes).tolist()
    for j,thisaxis in enumerate(pulseaxes):
        if thisaxis in plotdata.dimlabels:
            x =  plotdata.fromaxis(thisaxis,lambda x: x).data
        else:
            x = selected_pulses[thisaxis][0]
        pulseindex = pulseaxes_sorted.index(thisaxis)
        coherence_changes[:,pulseindex] = (x*ones(netshape)).flatten()
    if verbose: print "(show_pathways): coherence changes",coherence_changes
    for j in range(coherence_changes.shape[0]):
        #print "about to pass list",coherence_changes[j,:]
        plot_coherence_diagram(ax1,coherence_changes[j,:],
                [plotdata.data.shape[plotdata.axn(k)] if k in plotdata.dimlabels 
                    else selected_pulses[k][1]
                    for k in pulseaxes_sorted],
                j,coherence_changes.shape[1])
    axis('tight')
    ax1.set_ylim(-1*coherence_changes.shape[0],0)
    axis('off')
    return
def diagnostic_plot(fl,plotdata,figname = 'diagnostic',show_abs = True,selected_pulses = None,normalized = False,interpolation = 'bicubic'):
    def add_coherence_diagrams(plotdata,selected_pulses = selected_pulses):
        ax1 = axes([0.1,0.1,0.25,0.8]) # l b w h
        show_pathways(ax1,plotdata,selected_pulses = selected_pulses)
        ax2 = axes([0.35,0.1,0.6,0.8]) # l b w h
        ax2.yaxis.set_label_position('right')
        return ax1,ax2
    if show_abs:
        fl.next(figname,figsize = (8,20))
        ax1,ax2 = add_coherence_diagrams(plotdata,selected_pulses = selected_pulses)
        fl.image(abs(plotdata).cropped_log(),interpolation = interpolation)
    fl.next(figname+', with phase info',figsize = (8,20))
    ax1,ax2 = add_coherence_diagrams(plotdata,selected_pulses = selected_pulses)
    #imgplot = fl.image(plotdata.cropped_log(),interpolation = 'bicubic',ax = ax2)subtracted
    if normalized:
        newdata = plotdata.copy().cropped_log(subplot_axes = [j for j in plotdata.dimlabels if j[:5]!='phcyc'])
    else:
        newdata = plotdata.copy().cropped_log()
    imgplot = fl.image(newdata,interpolation = interpolation,ax = ax2)
    return
def open_cw_file(filename,fl = None,use_sweep = False,**kwargs):
    raise RuntimeError("open_cw_file is deprecated -- use find_file with no post-processing instead")
def load_nutation_curve(main,background = None,fl = None,max_freq = 30e6,deadtime = 0e-6,intwidth = 2e6,**kwargs):
    if fl is None:
        fl = figlist_var()
    data,fl = find_file(main,postproc = postproc_b1_fid_file,**kwargs)
    if background is not None:
        data_bg,fl = find_file(background, postproc = postproc_b1_fid_file, fl = [fl,'background'],**kwargs) 
        fl.next('subtracted')
        difference = data-data_bg
        fl.image(difference)
        print 'units are',fl.units[fl.current]
        if fl.units[fl.current][0] == 'ns':
            if fl.black:
                color = '#ffffff'
            else:
                color = 'k'
            axvline(x = deadtime/1e-9,alpha = 0.9,linewidth = 3,color = color)
        else:
            raise ValueError('units of plot were not ns')
    else:
        difference = data
    fl.next('ft')
    difference = difference['t2',lambda x: x>deadtime]
    difference.ft('t2')
    difference = difference['t2',lambda x: abs(x)<max_freq]
    difference.set_units('t2','Hz')
    difference.set_units('plen','s')
    double_ft = difference.copy()
    difference.human_units()
    fl.image(difference)
    if fl.units[fl.current][0] == 'MHz':
        if fl.black:
            color = '#ffffff'
        else:
            color = 'k'
        axvline(x = -intwidth/1e6,alpha = 0.5,linewidth = 3,color = color)
        axvline(x = +intwidth/1e6,alpha = 0.5,linewidth = 3,color = color)
    fl.next('ft')
    fl.image(difference)
    double_ft = double_ft['plen':(4e-9,)] # because the other stuff is junk
    double_ft.ft('plen',shift = True,pad = 512)
    if double_ft.get_units('plen') == 'Hz':
        x = double_ft.getaxis('plen')
        x[:] /= gammabar_e/1e4 #2.807 MHz/G
        double_ft.set_units('plen','G')
    else:
        raise ValueError('At this point, units should be in Hz, but they are in %s!!'%double_ft.get_units('plen'))
    double_ft = double_ft['plen',lambda x: abs(x)<50.]
    double_ft.human_units()
    fl.next('abs double ft (%s)'%main)
    fl.image(abs(double_ft))
    gridandtick(gca())
    fl.next('FIDs')
    difference.reorder('t2','plen')
    fl.plot(abs(difference))
    fl.next('nutation curve:\nfor signal $\pm %0.1f$ MHz about resonance'%(intwidth/1e6),figsize = (7,4))
    forplot = difference['t2',lambda x: abs(x) < intwidth].copy()
    forplot = forplot['plen':(5.,)]
    forplot_abssum = abs(forplot)
    forplot.run(sum,'t2')
    forplot_abssum.run(sum,'t2')
    forplot.data *= phaseopt(forplot['plen':(5.,)].data)
    forplot.data *= sign(forplot['plen':5.].data)
    fl.plot(forplot.runcopy(imag),'go-',alpha = 0.25,label = 'imag')
    fl.plot(forplot.runcopy(real),'bo-',alpha = 0.5,label = 'real')
    forplot.run(abs)
    fl.plot(forplot,'ko-',alpha = 0.25,label = 'abs')
    forplot_abssum /= forplot_abssum.data.max()
    forplot_abssum *= forplot.data.max()
    fl.plot(forplot_abssum,'ko-',alpha = 0.15,label = 'int(abs)')
    gridandtick(gca())
    expand_y()
    #xlims = gca().get_xlim()
    #xlims = list(xlims)
    #xlims[0] = 5. #set to 5 ns, but don't delete the data, so I can pan over
    #xlim(xlims)
    return difference,fl

def oscilloscope_data(*args):
    r"""read in a set of data that was stored on the oscilloscope using record_scope_data.pyw
    :param arguments: one or two file names.  if two are given, one should have `_bg` in the filename, indicating that it is a background scan (taken with the Tx attenuator turned all the way up)
    :type arguments: strings
    :return: a tuple of data:
        - copolarized reflection
        - modulator power
        - tuning signal (complex for I/Q)
    :rtype: tuple of 3 nddata objects
    """
    if len(args) == 2:
        file1 = getDATADIR('95GHz','agilent_scope')+args[0]
        file2 = getDATADIR('95GHz','agilent_scope')+args[1]
        a = scipy.io.loadmat(file1) # because 0.17 crashes in anaconda, needed to downgrade anaconda scipy with conda install scipy=0.16
        b = scipy.io.loadmat(file2)
        if '_bg' in args[0]:
            mdata_bg,mdata = a,b
        elif '_bg' in args[1]:
            mdata_bg,mdata = b,a
        else:
            raise RuntimeError("can't figure out which is the background")
        for variable in ['ch%d'%(j+1) for j in range(4)]:
            mdata[variable] -= mdata_bg[variable]
    elif len(args) == 1:
        mdata = scipy.io.loadmat(args[0])
    else:
        raise RuntimeError("don't know what to do with that!")
    copol_data = nddata(mdata['ch1'],[-1],['t']).labels('t',mdata['time'].flatten()).set_units('t','s')
    mod_data = nddata(mdata['ch3'],[-1],['t']).labels('t',mdata['time'].flatten()).set_units('t','s')
    tune_data = mdata['ch2'].flatten()+1j*mdata['ch4'].flatten()
    tune_data = nddata(mdata['ch2']+1j*mdata['ch4'],[-1],['t']).labels('t',mdata['time'].flatten()).set_units('t','s')
    tune_data.name('Voltage')
    if real(tune_data.data).sum() < 0:
        tune_data *= -1
    return copol_data,mod_data,tune_data
def find_attenuation(basename,
        attenuated_filename,
        unattenuated_filename,
        background_filename,
        fl):
    r'''This takes two identical oscilloscope readings that are different only by the attenuation, 
    and it finds the attenuation
    :returns: the attenuation (>1, on a linear scale)
    :param fl: it adds a plot to figure list fl that shows the
        attenuated and non-attenuated data, and allows you to check
        that the modulator and copolarized reflection are the same and
        that the attenuation correctly relates the detector channels
        for both measurements'''
    copol,mod,tune1 = oscilloscope_data(attenuated_filename,background_filename)

    fl.next('determine attenuation')
    fl.plot(copol*6,alpha = 0.5,label = 'copol. refl.',color = '#888800',linewidth = 2)
    fl.plot(mod,alpha = 0.5,label = 'modulator',color = '#8800ff',linewidth = 2)
    fl.plot(abs(tune1),'k',label='detector -- abs',alpha = 0.5,linewidth = 2)
    #fl.plot(tune1.runcopy(real),'b',label='detector -- real',alpha = 0.5,linewidth = 1)
    gridandtick(gca())

    copol,mod,tune2 = oscilloscope_data(unattenuated_filename,background_filename)

    fl.plot(copol*6,alpha = 0.5,label = 'copol. refl.',color = '#888800',linewidth = 2)
    fl.plot(mod,alpha = 0.5,label = 'modulator',color = '#8800ff',linewidth = 2)
    fl.plot(abs(tune2),'k',label='detector -- abs',alpha = 0.5,linewidth = 2)
    #fl.plot(tune2.runcopy(real),'b',label='detector -- real',alpha = 0.5,linewidth = 1)
    gridandtick(gca())

    ratio = (abs(tune2)**2).mean('t') / (abs(tune1)**2).mean('t')
    ratio = sqrt(ratio.data)

    print "scaleup for this attenuation is",ratio
    fp = tables.openFile(get_notebook_dir()+'reflection_tests.h5',
            mode = 'a',
            title = 'reflection tests')
    d = {'exp':basename,
                'ratio':ratio}
    try:
        tablenode = h5table(fp.root,'attenuation_calibration',None)
        new_table = False
    except:
        new_table = True
    if new_table:
        print "trying to add row"
        h5addrow(fp.root,
                'attenuation_calibration',
                d)
        print "created a new table for",d
    else:
        t = tablenode.read()
        if basename in t['exp']:
            mask = t['exp'] == basename 
            if t[mask] == ratio:
                print "value already entered"
            else: # remove the existing and add this one
                h5remrows(fp.root,
                        'attenuation_calibration',
                        {'index':t[mask]['index'][0]})
                h5addrow(fp.root,
                        'attenuation_calibration',
                        d,
                        verbose = True)
        else:
            h5addrow(fp.root,
                    'attenuation_calibration',
                    d,
                    verbose = True)
    fp.close()
    fl.plot(abs(tune1)*ratio,'r',label='test ratio -- rescaled',alpha = 0.3,linewidth = 2)
    return ratio
def echo_display(data,
        fl = None,
        decimation = 4,
        bw = 200e6,
        plot_title = '',
        truncate_t2 = None,
        ):
    r'''Prepare frequency domain data and display the
    time-domain decay of the echo signal as a function of :math:`t_1`.

    It a assumes that the data is given as Fourier conjugate domains of
    :math:`t_1` and :math:`t_2`, *i.e.*, ``t1`` and ``t2``.
    Likely, this means that it was created with :func:`load_formatted_data <pyspecdata.acert_hdf5.load_formatted_data>`

    It follows this procedure:

        * ift the data
        * use the last quarter of the data in order to determine the extra peak that occurs at zero frequency along the direct dimension, and subtract it out
        * reorder the dimensions to put the direct dimension first
        * filter to select only :math:`\pm \text{bw}` along :math:`t_2`
        * reduce the number of points along `t1` while improving their SNR:
            * convolve by 1ns (should later just be the sampling time) mult. by `decimation`
            * sample one every `decimation` points
        * manually rename axes as ..math::`f (direct)` and ..math::`t_{echo}`, with units of MHz and ns (respectively)
        * generate a waterfall plot

    Parameters
    ----------
    data : nddata
        input data
    bw : float, optional
        bandwidth over which to plot ..math::`t_2`
    fl : figlist_var
        the figure to which it plots
    decimation : int, optional
        the decimation (downsampling) that is used -- note that we throw
        away points without throwing away SNR
    plot_title : str, optional
        title of the plot -- usually the experiment name
    truncate_t2 : double, optional
        truncate the t2 dimension to this value (typically, the value after
        which only noise remains)

    Returns
    -------
    nddata
        the data shown in the waterfall plot
    '''
    if fl is None:
        raise ValueError("You need to pass a figure list")
    data = data.copy()
    data.ift('t1')
    #{{{ subtract the axial peak
    value = data['t2':0]
    value = value['t1':(data.getaxis('t1').max()*0.75,inf)]
    print 'mean value:',value.mean('t1').set_error(None)
    data['t2':0] -= value
    #}}}
    if truncate_t2 is not None:
        data.ift('t2')
        data = data['t2',lambda x: x<truncate_t2]
        data.ft('t2')
    data.reorder(['t2','t1'])# put t2 first
    data = data['t2',lambda x: abs(x)<bw]
    #{{{ here, I use the approximation that a gaussian is "about 0" at 6 sigma
    data.convolve('t1',1e-9/decimation)
    data = data['t1',decimation::decimation]
    #}}}
    #fl.image(data)
    x = data.getaxis('t2')
    x /= 1e6
    data.rename('t2','f (direct) / MHz')
    x = data.getaxis('t1')
    x /= 1e-9
    data.rename('t1','$t_{echo}$ / ns')
    data.name('signal (arb.u)')
    fl.next('3D view of abs data\n'+plot_title)
    abs(data).waterfall(alpha = 0.8,color = 'w',rotation = (60,20))
    plot_title_str = plot_title
    if fl.basename is not None:
        plot_title_str += ' ('+fl.basename+')'
    title(plot_title_str)
    return data
def load_and_format(list_of_exp,
        zerofill = False,
        time_shifts = (0,0),
        field_dep_shift = None,
        deadtime = 0e-9,
        SW = (None,None),
        t2_limit = 100e-9,
        verbose = True,
        show_origin = True,
        transform = 'secsy',
        phaseplot_thresholds = (0,0),#phaseplot_thresholds = (0.075,0.1),
        fl = None):
    r""" Generate formatted signal that's manually phased through selection of the appropriate :math:`t_1` and :math:`t_2` time-shifts.

    Loads the file(s) given by `list_of_exp`, applying the phase corrections indicated by `time_shifts`.
    If there is a bin dimension, it chunks it to separate out any field dimension, and then averages along any remaining bin dimension.
    Finally, formats the data (*e.g.* inhomogeneity format or SECSY format) by applying the appropriate transform.

    .. versionchanged:: beta

        Previously, this selected out a particular field and/or :math:`T` value, but it doesn't do that anymore.

    To phase a spectrum from scratch:

    1. Start by calling with ``time_shifts=(0.0,0.0)``
    2. Look at the second plot, called "check timing correction".  Adjust the
        time shifts so that the beginning of the signal is labeled as
        approximately :math:`t_1=0` :math:`t_2=0`.

        * A negative shift corresponds to moving the signal left (up).
        * Note that the :math:`x`-axis is :math:`t_2-t_1`, so it's easiest to adjust the :math:`t_1` time-shift first, before messing with :math:`t_2`.

    3. Look at the plot that shows the phased data with :math:`t_1` in the time domain, 
        and if there are :math:`>1` cycles in the phase (color)
        along :math:`\mathcal{F}[t_2]`.

        * Add a positive time-shift to correct for a
            left --> right roll of magenta --> red --> yellow (since this is
            positive increase).
        * The time-shift can be calculated as (phase shift [cyc])/(range of x-axis [cyc/s])
    4. Look at the plot that gives
        the phase *vs.* :math:`\mathcal{F}[t_2]`,
        and apply any additional corrections:

        * Positive phase corrects positive slope -- the 2 here adjusts for the y units (pi rad) vs. cycles
        * The time-shift can be calculated as (0.5 [cyc/:math:`\pi` rad])\*(phase shift [:math:`\pi` rad])/(range of x-axis [cyc/s])
        * For causal signal, the red lines show the phase where the signal should start (top/left) and end (bottom/right).  However (**why?**), the phase may turn around and come back to zero shortly after it reaches the furthest point from zero.
    5. In order to apply the same correction along :math:`t_1`, it's typically necessary to set a reasonable ``SW`` for :math:`t_1` and to set ``zero_fill`` to ``True``.

    Parameters
    ----------

    list_of_exp : list of 3-tuples
        ``[(e1,d1,b1),(e2,d2,b2),...]`` where

        :`e1`: experiment type (*eg.* ``echo_t2``)
        :`d1`: date (regexp, no leading/trailing ``.*``)
        :`b1`: file base name (regexp, gives the pattern up to the end of the file)
    time_shifts : tuple pair
        a :math:`(t_1,t_2)` pair of numbers that are used to determine
        the linear phase correction.  These are *in addition* to a
        :math:`\frac{\pi}{2}t_{pulse}` correction that is applied to correct
        for evolution during the pulse.
    field_dep_shift : double
        an additional, field-dependent phase shift
        (though it's not entirely clear to me why this should be allowed)
    SW : tuple pair of tuple pairs
        ``((f1_start,f1_stop),(f2_start,f2_stop))`` is a tuple pair of tuple pairs, where 
        ``(f1_start,f1_stop)`` gives the limits on the spectral window for
        :math:`\mathcal{F}[t_1]`
        and ``(f2_start,f2_stop)`` gives those for
        :math:`\mathcal{F}[t_2]`
    t2_limit : double
        Defaults to 100 ns.  Most spectra will die off relatively quickly, so
        only plot t2 up to a certain value.
    show_origin : bool, optional
        Show white cross-hairs at the origin.
        Defaults to true.
    deadtime : double
        Throw out the signal before `deadtime`.
    transform : str (case insensitive)
        Determines the final format of the time-domain data.  Note that `SECSY`
        will only pull one coherence pathway :math:`S_{c-}`, while
        `inhomogeneity` will pull the full hypercomplex set (:math:`S_{c+}`
        stacked on top of :math:`S_{c-}`)

        :SECSY: Apply the SECSY transform, which selects the :math:`S_{c-}` and throws out the first half of the echo.  Performs skew with :func:`skew <pyspecdata.skew>` to make sure there is no aliasing in either the current or Fourier conjugate dimension.
        :manual SECSY: Like SECSY, but rather than calling skew manually applies a frequency-dependent phase shift to accomplish the skew, potentially leading to aliasing. 
        :inh: Applies the inhomogeneity transform (a :math:`45^\circ` rotation, followed by mirroring the data that falls to the left of the axis).
    fl : figlist_var
        Figure list used to generate plots that aid in manual phasing.

        .. note::
            The following parameters are generated on the fly and used to determine these plots.  The corresponding variables store the *indices* (rather than axis labels) corresponding to the values described below:

            :fields_toplot: select 5 evenly spaced fields that span the range of the "fields" axis labels, and then select the central 3 out of those.
            :field_at_max: take the abs, mean both :math:`t_1` and :math:`t_2` and find the field with max signal.
            :echos_toplot: take the abs, mean both :math:`t_1` and :math:`fields` and find the time with max signal, then take five evenly spaced points from there to the end, and drop the last one.
            :F2_toplot: Take the mean of the absolute value along :math:`\mathcal{F}[t_1]`, then use :func:`contiguous <pyspecdata.nddata.contiguous>` to select the largest peak along :math:`\mathcal{F}[t_2]`, then then select 5 evenly spaced frequencies along that peak.  Determine this only for the first field indexed by `fields_toplot`.
            
        The plots are:

        :signal slice before timing correction: A cropped log plot of the selected coherence pathway
        :(pre-transform) check timing correction: (0,0) is aliased to center and echo is sheared.  By default (if ``show_origin==True``), shows cross-hairs that should align with the apparent time origin of the signal.

            * show signal only for `field_at_max`

        :after transform: Shows the time-domain signal after the selected transform has been applied.
        :phased data: there are several plots with this label. The first three are used to check the phase roll along :math:`t_2`:.  These are all plotted without interpolation.

            * a plot that is in the frequency domain.
                * White lines are used to mark the positions given by `F2_toplot`.
                * Selects a subset of data given by `fields_toplot`.
            * a plot that is in the time domain along :math:`t_1`
                * Selects a subset of data given by `fields_toplot`.
            * a frequency domain cropped-log plot
                * Selects a subset of data given by `fields_toplot`.

            Then there is *one* of the following, which is used to show what the final output looks like:

            * select real, for pure absorption
            * double real ft

        :plot the phase: Finally, for detailed phasing are the following two plots. These phase plots use :func:`contiguous <pyspecdata.nddata.contiguous>` to selectively plot the phase over the frequencies that are a quarter of the maximum signal:

            * plot the phase at different :math:`t_1`.
                * Selects a subset of the data as indicated by `fields_toplot` and `echos_toplot`
            * plot the phase at different :math:`\mathcal{F}[t_2]`.
                * Selects a subset of the data as indicated by `fields_toplot` and `F2_toplot`.

    """
    plt_chk_t1ph = 'phased data (lines check phasing along $t_1$)'
    plt_chk_t2ph = 'phased data -- time domain $t_1$\nto check phasing along $t_2$'
    # {{{ load the parameters
    if fl is None:
        fl = figlist_var()
    retval_list = []
    t1_timeshift, t2_timeshift = time_shifts
    transform = transform.lower()
    # }}}
    for j,info in enumerate(list_of_exp):
        # {{{ load the data and make sure that it's set up for shifted ft
        exp_type,date,short_basename = info
        fl.basename = short_basename+'\n' #almost forgot about this -- only need to do this once, and it adds to all the other names!
        if SW[1] is None:
            tempkwargs = {}
        else:
            tempkwargs = {'prefilter':SW[1]} # prefilter leaves it shifted
        data = find_file(date+'.*[\-_]%s\.'%short_basename,exp_type = exp_type,
                **tempkwargs)
        if SW[1] is None:
            data.ft('t2',shift = True)# this seems redundant with the next line, but needed to ensure that subsequent shifts give symmetric signal
        data.ift('t2')
        # }}}
        #{{{ select the appropriate coherence pathways (depending on the experiment) to generate a uniformly formatted t1 x t2 2D experiments
        if 'bin' in data.dimlabels:
            if verbose: print "taking the mean along multiple bins"
            data.chunk_auto('bin','fields')
        dropped_labels = data.squeeze()
        if 'te' in dropped_labels or 't1' in dropped_labels:
            has_indirect = False
        else:
            has_indirect = True
        if 'bin' in data.dimlabels:
            if verbose: print "taking the mean along multiple bins"
            data.run(mean,'bin')
        if 't1' in data.dimlabels:
            assert data.get_units('t1') == 's'
        if 'T' in data.dimlabels:
            assert data.get_units('T') == 's'
        if data.get_prop('description_class') in ['ELDOR','ELDOR_3D']:
            signal_slice = data['phcyc1',1,'phcyc2',-1,'phcyc3',-1]
        elif data.get_prop('description_class') == 'echo_T2':
            signal_slice = data['phcyc1',1,'phcyc2',-2]
            if has_indirect:
                signal_slice.rename('te','t1')
                signal_slice.set_units('t1','s')
        else:
            raise ValueError("I can't deal with this type of experiment!")
        #}}}
        # {{{ throw out the signal before `deadtime`
        if deadtime > 0:
            signal_slice = signal_slice['t2':(deadtime,inf)]
        # }}}
        fl.next('signal slice\nbefore timing correction (cropped log)')
        fl.image(signal_slice.copy().cropped_log(), interpolation = 'bicubic')
        ax = gca(); ax.title.set_fontsize('small')
        #{{{ set up the values of the indirect dimensions that we iterate over in the phase plots
        if has_indirect:
            time_w_max_signal = abs(signal_slice).set_error(None).mean('fields').mean('t2').argmax('t1').data
            echos_toplot = r_[time_w_max_signal:
                    signal_slice.getaxis('t1')[-1]:
                    5j][:-1]
            # need to convert to indices only after we zero fill
        else:
            echos_toplot = [None]
        if 'fields' in signal_slice.dimlabels:
            fields_toplot = signal_slice.getaxis('fields')
            if len(fields_toplot) > 3:
                fields_toplot = r_[fields_toplot[0]:fields_toplot[-1]:5j][1:-1]
            fields_toplot = signal_slice.indices('fields',fields_toplot)
            field_at_max = abs(signal_slice).mean('t1').mean('t2').argmax('fields', raw_index=True).data
        else:
            fields_toplot = [None]
            field_at_max = None
        #}}}
        #{{{ shift the time axis to account for evolution during the pulse
        if signal_slice.get_prop('t90') is not None:
            t2_timeshift -= signal_slice.get_prop('t90')*2/pi
            if verbose: print "pulse length timeshift is",signal_slice.get_prop('t90')*2/pi
        elif signal_slice.get_prop('pulse_length') is not None:
            t2_timeshift -= signal_slice.get_prop('pulse_length')*2/pi
            if verbose: print "pulse length timeshift is",signal_slice.get_prop('pulse_length')*2/pi
        else:
            raise ValueError("I can't find the length of the pulse (to apply the appropriate phase correction")
        #}}}
        # {{{ apply the t2 timeshift
        signal_slice.setaxis('t2',lambda x: x+t2_timeshift)
        # }}}
        # {{{ apply the field-dependent shift
        signal_slice.ft('t2')
        if field_dep_shift is not None:
            signal_slice *= signal_slice.fromaxis('fields',
                    lambda f: exp(-1j*2*pi*field_dep_shift*f)) # negative shift corrects negative slope
        signal_slice.ift('t2')
        # }}}
        #{{{ apply filtering and timeshift along t1
        if has_indirect:
            if verbose: print "the starting value of t1 is",signal_slice.getaxis('t1')[0]
            original_t1_max = signal_slice.getaxis('t1')[-1]
            print "DEBUG: dt for t1",diff(signal_slice.getaxis('t1')[r_[0,1]])
            if zerofill:
                signal_slice.ft('t1',shift = True,pad = 512)
                print "DEBUG 1: range of t1",signal_slice.getaxis('t1')[r_[0,-1]]
            else:
                signal_slice.ft('t1',shift = True)
                print "DEBUG 1: range of t1",signal_slice.getaxis('t1')[r_[0,-1]]
            if SW[0] is not None:
                signal_slice = signal_slice['t1':SW[0]]
            signal_slice.ift('t1')
            print "DEBUG again: dt for t1",diff(signal_slice.getaxis('t1')[r_[0,1]])
            signal_slice.setaxis('t1',lambda x: x+t1_timeshift)
            echos_toplot = signal_slice.indices('t1',echos_toplot)
        #}}}
        t_start = signal_slice.getaxis('t1')[0]
        if verbose: print 'applied a t1 timeshift: t_start=',t_start/1e-9,'manual=',t1_timeshift/1e-9
        del t_start
        #{{{ check that the timing correction looks OK
        fl.next('(pre-transform) check timing correction\n(0,0) is aliased to center and echo  is sheared')
        forplot = signal_slice.copy()
        print "t1 ft prop",forplot.get_ft_prop('t1')
        print "t2 ft prop",forplot.get_ft_prop('t2')
        forplot.shear('t2','t1',-1.0,method='linear')
        #forplot.ft('t2')
        #forplot.ft_clear_startpoints('t2',f='current')
        #forplot.ift('t2',shift = True) # get a centered echo here
        #if has_indirect:
        #    forplot.ft('t1') # ft along t1 so I can clear the startpoint, below
        #    forplot.ft_clear_startpoints('t1',f='current')
        #    forplot.ift('t1',shift = True) # generate a centered echo here
        sheared_t2_name = r'$t_2-t_1$' 
        forplot.rename('t2',sheared_t2_name)
        if field_at_max is not None:
            forplot = forplot['fields',field_at_max]
        fl.image(forplot, interpolation = 'nearest')
        ax = gca(); ax.title.set_fontsize('small')
        if show_origin:
            ax.axvline(x = 0,color = 'w',alpha = 0.25,linewidth = 3)
            ax.axhline(y = 0,color = 'w',alpha = 0.25,linewidth = 3)
        fl.next('check the sheared FT')
        forplot.ft(['t1',sheared_t2_name])
        fl.image(forplot)
        del forplot
        #}}}
        # at this point, signal_slice is in the time domain 
        # {{{ do whatever transformation we have selected
        fl.next('after transform')
        print "DEBUG: dt for t1 immediately before transform",diff(signal_slice.getaxis('t1')[r_[0,1]])
        if transform == 'manual secsy':
            signal_slice.secsy_transform_manual('t2','t1',has_indirect = has_indirect)
        elif transform == 'secsy':
            signal_slice.secsy_transform('t2','t1',has_indirect = has_indirect)
            print "DEBUG again 2: dt for t1",diff(signal_slice.getaxis('t1')[r_[0,1]])
            print "DEBUG FT startpoint",signal_slice.get_ft_prop('t1',['start','freq'])
        if has_indirect:
            fl.image(signal_slice['t1':(None,original_t1_max)]['t2':(None,t2_limit)])
        else:
            fl.image(signal_slice['t2':(None,t2_limit)])
        # }}}
        # {{{ before continuing, re-apply the filters, if needed
        # (because I fill when performing a skew, this can be necessary)
        signal_slice.ft(['t1','t2'])
        print "DEBUG 2: range of t1",signal_slice.getaxis('t1')[r_[0,-1]]
        for j,val in enumerate(SW):
            if val is not None:
                signal_slice = signal_slice['t{:d}'.format(j+1):SW[j]]
        # }}}
        #{{{ correct the zeroth order
        phase = signal_slice.copy().mean_all_but([None]).data
        assert isscalar(phase),' '.join(
                map(repr,["Error, phase is",phase])) # to make
                #    sure I'm not phasing anything differently
        phase /= abs(phase)
        signal_slice /= phase
        #}}}
        #{{{ various plots to check the phasing
        fl.next(plt_chk_t1ph)
        fl.image(signal_slice['fields',fields_toplot])
        fl.next(plt_chk_t2ph)
        if has_indirect:
            fl.image(signal_slice.copy().ift('t1')['fields',fields_toplot,'t1':(None,original_t1_max)
                ].reorder('t1',first = False))
        else:
            fl.image(signal_slice.copy().ift('t1')['fields',fields_toplot])
        fl.next('phased data -- time domain $t_2$\ncropped log -- to check phasing along $t_1$')
        fl.image(signal_slice.copy().ift('t2')['fields',fields_toplot,'t2':(None,t2_limit)
            ].reorder('t1',first = False))
        fl.grid()
        if transform == 'inh':
            fl.next('phased data\nselect real, for pure absorption')
            fl.image(signal_slice.runcopy(real),interpolation = 'bicubic')
        elif transform == 'secsy':
            fl.next('phased data\ndouble real ft')
            forplot = signal_slice.copy().ift(['t1','t2']).ft('t1').run(lambda x: complex128(real(x))).ft('t2').run(real)
            fl.image(forplot,interpolation = 'bicubic')
            del forplot
        #}}}
        # {{{ generate the linear plots that I use to adjust the timing corrections
        def quarter_of_max(arg,axis):# returns a mask for half of maximum
            return arg > 0.25*arg.runcopy(max,axis)
        #{{{ check phasing along t2: phase plot from 3553 -- modified to show the magnitude reference
        forplot_ini = signal_slice.copy()
        forplot_ini.ift('t1')
        fl.next('plot the phase at different $t_1$',legend = True)
        forplot_ini.reorder('t2') # make t2 x
        max_abs = abs(forplot_ini.data).flatten().max()
        for this_field_idx in fields_toplot:
            forplot = forplot_ini.copy()
            # {{{ select the slice out of forplot_ini that I want
            field_label = ''
            if this_field_idx is not None:
                field_label += '%.2f T'%(forplot.getaxis('fields')[this_field_idx])
                forplot = forplot['fields',this_field_idx]
            # }}}
            for echo_idx in echos_toplot:
                this_label = []
                if len(field_label) > 0: this_label.append(field_label)
                # {{{ show the magnitude reference
                magnitude_reference = abs(forplot).mean('t1').set_error(None)
                fl.plot(magnitude_reference/max(magnitude_reference.data)-0.5,
                        'k', linewidth=3,
                        alpha=0.1)
                # }}}
                if echo_idx is not None:
                    this_echo_time = forplot.getaxis('t1')[echo_idx]/1e-9
                    this_label.append('%d ns'%(this_echo_time))
                    # {{{ show where I pull the slices
                    fl.push_marker()
                    fl.next(plt_chk_t2ph)
                    gca().axvline(x = this_echo_time,color = 'w',alpha = 0.25)
                    fl.pop_marker()
                    # }}}
                fl.plot(forplot['t1',echo_idx][lambda x:
                    abs(x) > phaseplot_thresholds[1]*max_abs
                    ].runcopy(angle)/pi,'.',alpha = 0.5,
                    markersize = 3,label = ', '.join(this_label))
        fl.phaseplot_finalize()
        del forplot,forplot_ini
        #}}}
        #{{{ check phasing along t1 as a function of t2 frequency offset
        if has_indirect:# there is a t1 dimension
            #{{{ phase plot from 3553 -- modified to add amplitude
            fl.next('plot the phase at different $F_2$',legend = True)
            for j,this_field_idx in enumerate(fields_toplot):
                forplot = signal_slice
                this_label = []
                if this_field_idx is not None:
                    this_label.append('%.2f T'%(forplot.getaxis('fields')[this_field_idx]))
                    forplot = forplot['fields',this_field_idx]
                # {{{ show the magnitude reference
                magnitude_reference = abs(forplot).mean('t2').set_error(None)
                fl.plot(magnitude_reference/max(magnitude_reference.data)-0.5,
                        'k', linewidth=2,
                        alpha=0.25)
                # }}}
                sum_for_contiguous = abs(forplot).mean('t1').set_error(None)
                f_start,f_stop = sum_for_contiguous.contiguous(
                        quarter_of_max,'t2')[0,:]#  pull the first row, which is
                #                                   the largest block
                if j == 0:
                    F2_toplot = forplot.indices('t2',r_[f_start:f_stop:5j])
                forplot.reorder('t1') # make t1 x
                max_abs = abs(forplot.data).flatten().max()
                this_label.append('')
                for this_freq_idx in F2_toplot:
                    this_freq = forplot.getaxis('t2')[this_freq_idx]
                    this_label[-1] = '%d MHz'%(this_freq/1e6)
                    # {{{ show where I pull the slices
                    fl.push_marker()
                    fl.next(plt_chk_t1ph)
                    gca().axvline(x = this_freq/1e6,color = 'w',alpha = 0.25)
                    fl.pop_marker()
                    # }}}
                    fl.plot(forplot['t2',this_freq_idx][lambda x:
                        abs(x) > phaseplot_thresholds[0]*max_abs].runcopy(angle)/pi, '.',
                        alpha = 0.5, markersize = 3,
                        label = ', '.join(this_label))
            fl.phaseplot_finalize()
            #}}}
        #}}}
        # }}}
        retval_list.append(signal_slice)
    return retval_list
def plot_oned_v_field(thisdata,
        oned_plot_name = '1D spectrum from SECSY',
        field = None,
        fl = None
        ):
    r"""take a dataset and plot it as a function of field
    if no field is specified, assume this is not a field-swept experiment, and pull the field from the dataset's properties
    to avoid issues, this does create a copy of the data
    (most of this is copied from 3946)

    `thisdata` can be a single slice with dimension `t2`, or it can be:
    - a list of nddata objects
        (here, each one is passed recursively to `plot_oned_v_field`)
    - an nddata with dimensions `t2`x`fields`
        (here, each field is passed and the with the `field` kwarg set appropriately)
    """
    if fl is None:
        fl = figlist_var()
    if type(thisdata) is list:
        for j in thisdata:
            plot_oned_v_field(j,
                    oned_plot_name = oned_plot_name,
                    field = field,
                    fl = fl)
        return
    elif 'fields' in thisdata.dimlabels:
        x = thisdata.getaxis('fields')
        for j in range(ndshape(thisdata)['fields']):
            plot_oned_v_field(thisdata['fields',j],
                    oned_plot_name = oned_plot_name,
                    field = x[j],
                    fl = fl)
        return
    else:
        thisdata = thisdata.copy() # because I screw with the axis
        #{{{ axis conversion
        if field is None:
            field = thisdata.get_prop('field')
        x = thisdata.getaxis('t2')
        x[:] = field + x / gammabar_e
        #}}}
        fl.next(oned_plot_name)
        lines = fl.plot(thisdata.runcopy(abs),'-',alpha = 0.5,label = 'abs')
        try:
            color = lines[-1].get_color()
        except:
            print "lines is",lines
            print "has attributes",dir(lines[-1])
            print ndshape(thisdata),thisdata.get_prop('FT')
            raise RuntimeError("problem with getting the color -- possibly, you  are passing data with too many dimensions")
        realy = thisdata.runcopy(real)
        fl.plot(realy,'--',color = color,alpha = 0.5,label = 'Re')
        fl.plot(thisdata.runcopy(imag),':',color = color,alpha = 0.5,label = 'Im')
        #{{{ just show the field this spectrum was acquired at
        y_at_field = realy['t2':field].data
        fl.plot(field,y_at_field,'o',color = color,alpha = 0.5)
        text(field,y_at_field,'%0.4f T'%field,alpha=0.5,color = color,size = 'xx-small',ha='left',va='bottom',rotation=45)
        #}}}
        xlabel(r'($B_0$ / $T$) + $\Delta f$ / ($%0.3f\times 10^{10}$ $\frac{Hz}{T}$)'%(gammabar_e/1e10))
        return
def plot_oned_v_offset(thisdata,
        oned_plot_name = '1D offset spectrum from SECSY',
        fl = None
        ):
    r"""this is similar to plot_oned_v_field, but it plots *vs.* offset, rather than field

    like for plot_oned_v_field, `thisdata` can be a single slice with dimension `t2`, or it can be:
    - a list of nddata objects
        (here, each one is passed recursively to `plot_oned_v_field`)
    - an nddata with dimensions `t2`x`fields`
        (here, each field is passed and the with the `field` kwarg set appropriately)
    """
    if fl is None:
        fl = figlist_var()
    if type(thisdata) is list:
        for j in thisdata:
            plot_oned_v_offset(j,
                    oned_plot_name = oned_plot_name,
                    fl = fl)
        return
    elif 'fields' in thisdata.dimlabels:
        x = thisdata.getaxis('fields')
        for j in range(ndshape(thisdata)['fields']):
            plot_oned_v_offset(thisdata['fields',j],
                    oned_plot_name = oned_plot_name,
                    fl = fl)
        return
    else:
        fl.next(oned_plot_name)
        lines = fl.plot(thisdata.runcopy(abs),'-',alpha = 0.5,label = 'abs')
        color = lines[-1].get_color()
        realy = thisdata.runcopy(real)
        fl.plot(realy,'--',color = color,alpha = 0.5,label = 'Re')
        fl.plot(thisdata.runcopy(imag),':',color = color,alpha = 0.5,label = 'Im')
        axis('tight')
        expand_y()
        return
