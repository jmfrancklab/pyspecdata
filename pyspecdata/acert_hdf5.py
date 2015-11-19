from .core import *
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
def plot_comparison(input_data,date,fl,phase = True):
    "Compare a series of cw experiments"
    list_of_cw_files, mod_amp = zip(*input_data)
    single_date = True
    if date is None:
        mod_amp,date = zip(*mod_amp)
        single_date = False
    short_basename = list_of_cw_files[0]
    for j in range(0,len(list_of_cw_files)):
        short_basename = list_of_cw_files[j]
        if single_date:
            data = cw('%s.*'%date + short_basename + '\\.', phase = phase)
        else:
            data = cw('%s.*'%date[j] + short_basename + '\\.', phase = phase)
        data /= mod_amp[j]
        if j == 0:
            data_shape = ndshape(data)
            data_shape.add_correctly((len(list_of_cw_files),'indirect'))
            collected_data = data_shape.alloc()
            collected_data.set_units('field',data.get_units('field'))
            collected_data.labels(['indirect','field',],[array(list_of_cw_files),data.getaxis('field')])
        collected_data['indirect',j] = data
    fl.next('cw plots -- real + imag',legend = True,boundaries = False)
    normalization = abs(collected_data.runcopy(real)).run(max,'field').run(max,'indirect').data
    normalization = max(r_[abs(collected_data.runcopy(imag)).run(max,'field').run(max,'indirect').data,normalization])
    ax = gca()
    color_cycle = ax._get_lines.color_cycle
    for j in range(0,len(list_of_cw_files)):# I end up looping back over and reloading because they might have different axes, though in an older format, I did this all in one batch
        short_basename = list_of_cw_files[j]
        if single_date:
            data = cw('%s.*'%date + short_basename + '\\.', phase = phase)
        else:
            data = cw('%s.*'%date[j] + short_basename + '\\.', phase = phase)
        data /= mod_amp[j]
        next_color = next(color_cycle)
        fl.plot(data.runcopy(real)/normalization + 2*j,
                alpha = 0.75,color = next_color,label = collected_data.getaxis('indirect')[j])
        autolegend()
        fl.plot(data.runcopy(imag)/normalization + 2*j,
                alpha = 0.5,color = next_color)
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
def automagical_phasecycle(data,verbose = False):
    "Use the phase cycle list to determine the phase cycles, and then ift them to return coherence skips"
    if verbose: print "shape of data",ndshape(data)
    if verbose: print data.other_info['phasecycle']
    phasecyc_origindeces = r_[0:data.other_info['phasecycle'].shape[0]]
    phase_cycles = data.other_info['phasecycle']
    this_dtype = phase_cycles.dtype.descr * phase_cycles.shape[1]
    #{{{ construct the new phase cycle set
    new_fields = ['phcyc%d'%(j+1) for j in range(phase_cycles.shape[1])]
    phase_cycles =  phase_cycles.view(this_dtype).squeeze()
    phase_cycles.dtype.names = new_fields
    #}}}
    prod_of_nuniqueinfield = array([len(unique(phase_cycles[j])) for j in phase_cycles.dtype.names]).prod()
    if verbose: print 'product of unique elements in each field',prod_of_nuniqueinfield
    if verbose: print 'actual unique tuples',len(unique(phase_cycles))
    redundancy = prod_of_nuniqueinfield/len(unique(phase_cycles))
    if verbose: print 'redundancy',redundancy
    final_field = phase_cycles.dtype.names[-1]
    if verbose: print "before rotation",data['indirect',0]
    if redundancy == 2:
        mask = phase_cycles[final_field] == 90
        mask |= phase_cycles[final_field] == 270
        if verbose: print "these match",phase_cycles[mask]
        phase_cycles[final_field][mask] -= 90 
        data['phcyc',mask] *= exp(-1j*pi/2) # pulse was rotated 90 forward, which makes signal 90 backwards
        if verbose: print "after rotation",data['indirect',0]
    if verbose: print "argsort:",argsort(phase_cycles,order = phase_cycles.dtype.names)
    if verbose: print "argsorted:"
    for j in phase_cycles[argsort(phase_cycles)]:
        if verbose: print j
    if verbose: print ndshape(data)
    #{{{ actually assign and chunk the phase cycle dimension
    data.setaxis('phcyc',phase_cycles)
    sorted_indeces = argsort(phase_cycles)
    if verbose: print 'sorted indeces are',sorted_indeces
    data = data['phcyc',sorted_indeces]
    if verbose: print data.getaxis('phcyc')
    #for phcyc_i in [phase_cycles.dtype.names[0]]: # for debugging, just do the first
    for phcyc_i in phase_cycles.dtype.names[:-1]:
        if verbose: print 'chunking out',phcyc_i
        if verbose: print 'with axis label',data.getaxis('phcyc')
        data.chunk_auto('phcyc',phcyc_i,dimname = 'phcyc')
        if verbose: print 'just did',phcyc_i,'and got',ndshape(data)
    data.rename('phcyc',phase_cycles.dtype.names[-1])
    if verbose: print 'finalized and got',ndshape(data)
    #}}}
    for j in phase_cycles.dtype.names:
        x = data.getaxis(j)
        data.setaxis(j,r_[0:1:1j*len(x)])
        data.ift(j)
    return data
def cw(file_regexp,phase = True,use_sweep = False):
    r'this opens the cw data, using search_freed_file and open_cw_file, and then autophases it'
    if use_sweep:
        otherdim = 'current'
    else:
        otherdim = 'field'
    print "looking for "+lsafen(search_freed_file(file_regexp,'cw'))
    data = open_cw_file(search_freed_file(file_regexp,'cw'),use_sweep = use_sweep)
    if 'repeats' in data.dimlabels:
        print "found",ndshape(data)['repeats'],"averages"
        data.mean('repeats')
    #{{{ phase the spectrum
    baseline = (data[otherdim,0:2].mean().data + data[otherdim,-3:-1].mean().data)/2.0
    data -= baseline
    if phase:
        testphases = r_[0:2*pi:100j]
        result = zeros_like(testphases)
        for j,ph in enumerate(testphases):
            test = data.copy()
            test *= exp(-1j*ph)
            test_denom = test.runcopy(imag)
            test_denom.mean(otherdim)
            test.run(real)
            test.mean(otherdim)
            result[j] = abs(test.data)**2/abs(test_denom.data)**2
        data *= exp(-1j*testphases[argmin(result)])
        index = argmax(abs(imag(data.data)))
        data *= -1*sign(imag(data.data[index]))
    if use_sweep:
        data.set_units('current','A')
    else:
        data.set_units('field','T')
    #}}}
    return data
def open_cw_file(filename,fl = None,use_sweep = False,**kwargs):
    if fl is None:
        fl = figlist_var()
    if type(fl) is list and len(fl) == 2:
        figname_append = fl[1]
        fl = fl[0]
        if figname_append[-1] != '_':
            figname_append += '_'
    else:
        figname_append = ''
    h5 = h5py.File(filename,'r')
    #{{{ set up the complex number the hard way, for good form
    data = empty_like(h5['experiment']['data.r'],dtype = complex128)
    data_introspect = data.view([('r',double),('i',double)])
    for i_or_r in ['i','r']:
        data_introspect[i_or_r] = h5['experiment']['data.'+i_or_r]
    #}}}
    if len(data.shape) == 1:
        if use_sweep:
            data = nddata(data,data.size,['current']).labels('current',array(h5['experiment']['sweep_currents']))
        else:
            data = nddata(data,data.size,['field']).labels('field',array(h5['experiment']['fields']))
    elif len(data.shape) == 2:
        if use_sweep:
            data = nddata(data,data.shape,['repeats','current']).labels('current',array(h5['experiment']['sweep_currents']))
        else:
            data = nddata(data,data.shape,['repeats','field']).labels('field',array(h5['experiment']['fields']))
    h5.close()
    return data
def postproc_blank(data):
    return data
def postproc_generic(data):
    data = automagical_phasecycle(data)
    return data
def postproc_eldor_3d(data):
    if data.get_prop('phasecycle') is None:
        warnings.warn("Warning!  There is no phase cycling information -- you should really fix this")
        data.chunk('phcyc',['phcyc1','phcyc2','phcyc3'])
        #data.rename('indirect','t_e')
        data.ift('phcyc3')
        data.ift('phcyc2')
        data.ift('phcyc1')
        data.rename('indirect','t1')
    else:
        data = automagical_phasecycle(data)
    return data
def postproc_echo_T2(data):
    if data.get_prop('phasecycle') is None:
        warnings.warn("Warning!  There is no phase cycling information -- you should really fix this")
        data.chunk('phcyc',['phcyc2','phcyc1']) # the second pulse is on the outside, and the first pulse is on the inside
        data.setaxis('phcyc1',r_[0:1:4j])
        data.setaxis('phcyc2',r_[0:1:4j])
        data.ift('phcyc1') # remember from 2889 that I have to ift
        data.ift('phcyc2')
    else:
        data = automagical_phasecycle(data)
    return data
def postproc_eldor_old(data,**kwargs):
    data = automagical_phasecycle(data)
    data.rename('indirect','t1')
    t1_min = data.get_prop('t1_min')
    t1_inc = data.get_prop('t1_inc')
    t1_steps = data.get_prop('t1_steps')
    data.labels({'t1':r_[0:t1_steps]*t1_inc+t1_min})
    data.set_units('t1','s')
    return data
def postproc_b1_fid_file(data,fl = None,**kwargs):
    if fl is None:
        fl = figlist_var()
    if type(fl) is list and len(fl) == 2:
        figname_append = fl[1]
        fl = fl[0]
        if figname_append[-1] != '_':
            figname_append += '_'
    else:
        figname_append = ''
    data.rename('indirect','plen')
    if 'L_min' in kwargs.keys() and kwargs['L_min'] is not None:
        L_min = kwargs['L_min']
    else:
        L_min = data.get_prop('L_min')
    in_field_sweep = data.get_prop('in_field_sweep_current') 
    obsn("In field sweep current is %0.6f A"%in_field_sweep)
    if 'L_inc' in kwargs.keys() and kwargs['L_inc'] is not None:
        L_inc = kwargs['L_inc']
    else:
        L_inc = data.get_prop('L_inc')
    if 'L_steps' in kwargs.keys() and kwargs['L_steps'] is not None:
        L_steps = kwargs['L_steps']
    else:
        L_steps = data.get_prop('L_steps')
    data.labels({'plen':r_[0:L_steps]*L_inc+L_min,'t2':r_[0:ndshape(data)['t2']]*1e-9,'phcyc':r_[0:4]})
    fl.next(figname_append+'raw_data')
    fl.image(data)
    fl.next(figname_append+'shifted and phase cycled')
    data.ft('t2',shift = True)
    data *= data.fromaxis(['t2','plen'],lambda x,y: exp(2j*pi*x*(y*0.5)))
    data.ift('t2')
    data.ft('phcyc')
    fl.image(data)
    return data['phcyc',1],fl
def postproc_eldor_file(data,fl = None):
    steps = ndshape(data)['indirect']
    if steps != data.get_prop('t1_steps'):
        raise ValueError('I am confused -- the size of your indirect dimension doesn\'t match the t1_steps parameter')
    x = data.get_prop('t1_min')+data.get_prop('t1_inc')*r_[0:steps]
    data.rename('indirect','t_1')
    data.labels('t_1',x)
    data.set_units('t_1','s')
    return data
def search_freed_file(searchstring,exp_type,
        directory = dirformat(getDATADIR()+'95GHz'),
        print_result = False,
        verbose = False):
    ("this searches in the 95GHz subdirectory of the data directory in order to"
    " find the appropriate file\n"
    "Returns\n"
    "-------\n"
    " str\n"
    "     The name of the file matching the search string\n")
    directory = dirformat(directory + exp_type)
    if os.path.isdir(directory):
        files = re.findall('.*' + searchstring + '.*','\n'.join(os.listdir(directory)))
    else:
        raise RuntimeError("I can't find the directory:\n%s\nin order to get a file that matches:\n%s"%(directory,searchstring))
    if len(files) == 0:
        raise ValueError("I can't find a file matching the regular expression {:s} in {:s}".format(searchstring,directory))
    else:
        if len(files) > 1:
            warnings.warn('found multiple files:\n'+repr(files)+'\nand opening last')
        elif print_result and verbose:
            obsn("found only one file, and loading it:"+repr(files))
    return directory + files[-1]
def find_file(searchstring,
            exp_type = 'b1_fid',
            postproc = None,
            print_result = True,
            verbose = False,
            prefilter = None,
            indirect_dimlabels = None,# in case "dimlabels" is not set properly, I can manually pass the value of "indirect_dimlabels"
            **kwargs):
    r'''find the file (assumed to be an h5 format file) given by the regular
    expression "searchstring" inside the directory "directory", and postprocess
    with the function "postproc," which is fed the nddata data and the source
    text as arguments, as well as the keyword arguments "kwargs".
    If "postproc" is not set, it's chosen based on the value of
    (h5 root).experiment.description['class']'''
    #{{{ actually find one file and load it into the h5 object
    if 'subdirectory' in kwargs.keys():
        raise ValueError("The `subdirectory` keyword argument is not longer valid -- use `exp_type` instead!")
    filename = search_freed_file(searchstring,exp_type,print_result = print_result,verbose = verbose)
    h5 = h5py.File(filename,'r')
    #}}}
    #{{{ set up the complex number the hard way, for good form
    data = empty_like(h5['experiment']['data.r'],dtype = complex128)
    data_introspect = data.view([('r',double),('i',double)])
    for i_or_r in ['i','r']:
        data_introspect[i_or_r] = h5['experiment']['data.'+i_or_r]
    #}}}
    #{{{ start with the dimensions used in the HDF5 file -- "indirect" is a placeholder (see below)
    if len(data.shape) == 4:
        dimlabels = ['bin','phcyc','indirect','t2']
    elif len(data.shape) == 3:
        dimlabels = ['phcyc','indirect','t2']
    elif len(data.shape) == 2:
        dimlabels = ['phcyc','t2']
    else:
        raise ValueError("I don't know how to interpret data of shape %d"%len(data.shape))
    data = nddata(data,list(data.shape),dimlabels)
    #}}}
    #{{{ put all further information into the nddata in a way such that it can be retrieved with "get_prop"
    data.set_prop(dict([(k,v) for k,v in h5['experiment'].attrs.iteritems() if k[0] != '_' and k[0] != 'dimlabels']))
    data.set_prop(dict([('execution_'+k,v) for k,v in h5['experiment']['execution'].attrs.iteritems() if k[0] != '_']))
    data.set_prop(dict([('description_'+k,v) for k,v in h5['experiment']['description'].attrs.iteritems() if k[0] != '_']))
    #}}}
    #{{{ finish the t2 axis
    t2_steps = data.get_prop('t2_steps')
    t2_inc = data.get_prop('t2_inc')
    if t2_steps is not None:
        if ndshape(data)['t2'] != t2_steps:
            raise ValueError("There is a problem, because the dimensions of your data don't match the value given by t2_steps")
        data.labels('t2',r_[0:t2_steps]*t2_inc)
    else:
        if verbose: print "warning, I couldn't find the t2 steps parameter"
        data.labels('t2',r_[0:ndshape(data)['t2']]*1e-9)
    data.set_units('t2','s')
    if prefilter is not None:
        data.ft('t2',shift = True)
        data = data['t2':prefilter]
    #}}}
    #{{{ now, pull the dimlabels that we want
    if 'dimlabels' in h5['experiment'].attrs.keys():
        dimlabels = h5['experiment'].attrs['dimlabels'].tolist()
    else:
        print "You didn't set dimlabels in your pulse program -- the data is ambiguous and I'm taking my best guess!!!"
    if indirect_dimlabels is not None:# in case it's not stored correctly in the file and I need to hard-code it
        idx = dimlabels.index('indirect')
        dimlabels = dimlabels[:idx] + indirect_dimlabels + dimlabels[idx+1:]
    #}}}
    #{{{ assign all the dimensions
    indirect_names_in_h5 = [k for k,v in h5['experiment'].iteritems() if k not in ['data.i','data.r','bin_switch_times','fields','sweep_currents'] and type(v) == h5py.Dataset]
    if verbose: print "dimlabels are",dimlabels
    expected_dimlabels = set(dimlabels) - {'phcyc','t2','bin'} # we expect to find these axes stored in the HDF5
    if verbose: print "expected dimlabels",expected_dimlabels
    common_dimensions = expected_dimlabels & set(indirect_names_in_h5)
    #{{{ deal with the case where dimlabels and the stored axes don't match up correctly 
    if verbose: print "common dimensions",common_dimensions
    unlabeled_indirect = False
    if len(indirect_names_in_h5) == 0:
        warnings.warn('Warning!! The indirect dimension is not labeled!\nDimensions of the data:%s\nfilename: %s'%(repr(ndshape(data)),filename))
        unlabeled_indirect = True
    if len(common_dimensions) != len(indirect_names_in_h5):
        h5_but_not_dimlabels = common_dimensions ^ set(indirect_names_in_h5)
        dimlabels_but_not_h5 = common_dimensions ^ expected_dimlabels
        if len(h5_but_not_dimlabels) == 1 and len(dimlabels_but_not_h5) == 1:
            real_name = h5_but_not_dimlabels.pop()
            dimlabels_name = dimlabels_but_not_h5.pop()
            warnings.warn("Warning! The dimlabels set in the HDF5 file have a dimension called %s and the arrays given in the HDF5 file have one called %s -- I'm assuming they're the same and using the name from the HDF5 file."%(dimlabels_name,real_name))
            common_dimensions |= {real_name} # so the axis gets labeled below
            dimlabels[dimlabels.index(dimlabels_name)] = real_name
        elif dimlabels_but_not_h5 == {'indirect'}:# I need to split up the indirect dimension
            raise ValueError("You want me to split the indirect dimension, but I don't know how!:\n\tDimensions in the file: %s\n\tDimensions listed in dimlabels %s\n\tDimensions of the data: %s"%(repr(indirect_names_in_h5),repr(dimlabels),repr(ndshape(data))))
        else:
            raise ValueError("This file has dimensions that I can't just automatically sort!:\n\tDimensions in the file: %s\n\tDimensions listed in dimlabels %s\n\tDimensions of the data: %s\n\tDimensions in hdf5 but not dimlabels: %s\n\tDimensions in dimlabels but not HDF5: %s"%(repr(indirect_names_in_h5),repr(dimlabels),repr(ndshape(data)),repr(h5_but_not_dimlabels),repr(dimlabels_but_not_h5)))
    #}}}
    if verbose: print "now dimlabels are",dimlabels,"and indirect_names_in_h5 is",indirect_names_in_h5
    #{{{ chunk the "indirect" dimension up appropriately, and assign the dimensions
    if not unlabeled_indirect:
        if 'indirect' in data.dimlabels:# note that if there was only a single indirect dimension, it's already been used up
            chunk_dict = dict([(j,len(h5['experiment'][j])) for j in dimlabels if j not in ['bin','phcyc','t2']])
            if verbose: print ndshape(data)
            if verbose: print chunk_dict
            data.chunk('indirect',chunk_dict)
        for this_axis in common_dimensions:
            axis_data = array(h5['experiment'][this_axis])
            temp = ndshape(data)[this_axis]
            if len(axis_data) > temp:
                warnings.warn("Warning: the length of the axis '"+str(this_axis)+"' seems to exceed the length of the data!")
                axis_data = axis_data[0:temp]
            data.labels(this_axis,axis_data)
    #}}}
    #}}}
    #{{{ label the bin dimension, and chunk it up appropriately
    if 'bin' in data.dimlabels:
        #{{{ assign the labels for the bin dimension as a structured array
        forstruct = []
        forstruct_names = []
        for thisaxis in [x for x in 'bin_switch_times','fields','sweep_currents' if x in h5['experiment'].keys() and type(h5['experiment'][x]) is h5py.Dataset]:
            forstruct.append((h5['experiment'][thisaxis])[0])
            forstruct_names.append(thisaxis)
        print "structure is",forstruct,"names",forstruct_names,"zeros like",h5['experiment']['bin_switch_times'].shape
        x = make_rec(forstruct,forstruct_names,zeros_like = h5['experiment']['bin_switch_times'].shape)
        warnings.warn("detected a 'bin' dimension, and associating it with dtype "+repr(x.dtype))
        for thisaxis in forstruct_names:
            x[thisaxis] = h5['experiment'][thisaxis]
        data.labels('bin',x)
        #}}}
        if 'fields' in x.dtype.names:
            fields,indeces = unique(x['fields'],return_index = True)
            if verbose: print "number of fields",len(fields),"number of unique fields",len(x['fields']),'unique fields are',x['fields']
            if len(fields) != len(x['fields']):
                spacing = diff(indeces)
                if all(spacing == spacing[0]):
                    data.chunk('bin',['field','bin'],[spacing[0],len(x)/spacing[0]])
                else:
                    raise ValueError('The fields have some weird/inconsistent spacing')
        if data.getaxis('bin').dtype.names == ('bin_switch_times',):
            data.setaxis('bin',data.getaxis('bin')['bin_switch_times'])
            data.set_units('bin','s')
    #}}}
    #{{{ it always makes sense to put the phase cycle on the outside
    oldorder = list(data.dimlabels)
    data.reorder(['phcyc'])
    #}}}
    h5.close()
    temp = ndshape(data)
    if 'indirect' in temp.dimlabels and temp['indirect'] == 1:
        data = data['indirect',0]
    if postproc is None:
        if data.get_prop('description_class') == 'ELDOR':
            data = postproc_eldor_old(data)
        elif data.get_prop('description_class') == 'ELDOR_3D':
            data = postproc_eldor_3d(data)
        elif data.get_prop('description_class') == 'FID':
            data = postproc_generic(data)
        elif data.get_prop('description_class') == 'echo_T2':
            data = postproc_echo_T2(data)
        else:
            raise ValueError('postprocessing not defined for file with description-->class'+repr(data.get_prop('description_class')))
        return data
    else:
        return postproc(data,**kwargs)
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
    fl.next('abs ft')
    fl.image(abs(difference))
    double_ft = double_ft['plen':(4e-9,)] # because the other stuff is junk
    double_ft.ft('plen',shift = True,pad = 512)
    if double_ft.get_units('plen') == 'Hz':
        x = double_ft.getaxis('plen')
        x[:] /= 2.807e6 #2.807 MHz/G
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
        a = scipy.io.loadmat(getDATADIR('agilent_scope',args[0]))
        b = scipy.io.loadmat(getDATADIR('agilent_scope',args[1]))
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
    fp = tables.openFile(get_notebook_dir('reflection_tests.h5'),
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
        ):
    r'''Prepare frequency domain data for time-domain display:

        * ift the data
        * use the last quarter of the data in order to determine the extra peak that occurs at zero frequency along the direct dimension, and subtract it out
        * reorder the dimensions to put the direct dimension first
        * filter to select only ..math::`\pm` `bw`
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
            the decimation (downsampling) that is used -- note that we throw away points without throwing away SNR
        plot_title : str, optional
            title of the plot -- usually the experiment name

        Returns
        -------
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
    fl.next('3D view of abs data')
    abs(data).waterfall(alpha = 0.8,color = 'w',rotation = (60,20))
    title(plot_title)
    return data
def secsy_format(list_of_exp,
        zerofill = False,
        time_shifts = (0,0),
        field_dep_shift = None,
        deadtime = 29e-9,
        SW = (None,None),
        verbose = True,
        phaseplot_thresholds = (0,0),#phaseplot_thresholds = (0.075,0.1),
        fl = None):
    ("""Generate the "SECSY-format" signal that's manually phased -- """
     "`time_shifts` is a (t1,t2) pair of numbers\n"
     """`field_dep_shift` is an additional, field-dependent phase shift
     though it's not entirely clear to me why this should be allowed\n"""
     "`SW`=(`f1_limits`,`f2_limits`) is a tuple pair of tuple pairs and gives the limits on "
     "the spectral window for $\mathcal{F}[t_1]$ and $\mathcal{F}[t_2]$, respectively\n"
     "If there is a bin dimension, it chunks it to separate out any field "
     "dimension, and then averages along any remaining bin dimension\n\n"
     "loads the file(s) given by `list_of_exp`, which consists of a list of experiment type, date (regexp), basename (regexp) tuples, and where basename is assumed to give the pattern up to the end of the filename\n"
     "(upgrade note): previously, this selected out a particular field and/or T "
     "value, but it doesn't do that anymore")
    if fl is None:
        fl = figlist_var()
    retval_list = []
    t1_timeshift, t2_timeshift = time_shifts
    for j,info in enumerate(list_of_exp):
        exp_type,date,short_basename = info
        fl.basename = short_basename+'\n' #almost forgot about this -- only need to do this once, and it adds to all the other names!
        if SW[1] is None:
            tempkwargs = {}
        else:
            tempkwargs = {'prefilter':SW[1]} # prefilter leaves it shifted
        data = find_file(date+'.*[\-_]%s\.'%short_basename,exp_type = exp_type,
                **tempkwargs)
        if SW[1] is None:
            data.ft('t2',shift = True)
        data.ift('t2')
        #{{{ take various experiments, and convert them uniformly into t1 x t2 2D experiments
        if 'bin' in data.dimlabels:
            if verbose: print "taking the mean along multiple bins"
            data.chunk_auto('bin','fields')
        dropped_labels = data.squeeze()
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
            if 'te' not in dropped_labels:
                signal_slice.rename('te','t1')
        else:
            raise ValueError("I can't deal with this type of experiment!")
        #}}}
        if deadtime > 0:
            signal_slice = signal_slice['t2':(deadtime,inf)]
        fl.next('signal slice before timing correction (cropped log)')
        fl.image(abs(signal_slice).cropped_log(), interpolation = 'bicubic')
        #{{{ set up the values of the indirect dimensions that we iterate over in the phase plots
        if 'te' not in dropped_labels:
            echos_toplot = r_[signal_slice.getaxis('t1')[0]:
                    signal_slice.getaxis('t1')[-1]:
                    5j][1:-1]
        else:
            echos_toplot = [None]
        if 'fields' in signal_slice.dimlabels:
            fields_toplot = signal_slice.getaxis('fields')
        else:
            fields_toplot = [None]
        #}}}
        #{{{ phase correct and ft
        #{{{ shift back by the pulse length
        if signal_slice.get_prop('t90') is not None:
            t2_timeshift -= signal_slice.get_prop('t90')*2/pi
            if verbose: print "pulse length timeshift is",signal_slice.get_prop('t90')*2/pi
        elif signal_slice.get_prop('pulse_length') is not None:
            t2_timeshift -= signal_slice.get_prop('pulse_length')*2/pi
            if verbose: print "pulse length timeshift is",signal_slice.get_prop('pulse_length')*2/pi
        else:
            raise ValueError("I can't find the length of the pulse (to apply the appropriate phase correction")
        #}}}
        signal_slice.setaxis('t2',lambda x: x+t2_timeshift)
        signal_slice.ft('t2')
        #{{{ shift back by the echo time
        if 'te' in dropped_labels:
            echo_time = signal_slice.get_prop('te')
            signal_slice *= signal_slice.fromaxis('t2',lambda t2: exp(1j*2*pi*echo_time*t2)) # positive time shift corrects positive slope
        else:
            signal_slice *= signal_slice.fromaxis(['t1','t2'],lambda t1,t2: exp(1j*2*pi*t1*t2))
            t_start = signal_slice.getaxis('t1')[0]
        #}}}
        if field_dep_shift is not None:
            signal_slice *= signal_slice.fromaxis('fields',
                    lambda f: exp(-1j*2*pi*field_dep_shift*f)) # negative shift corrects negative slope
        if not 'te' in dropped_labels:
            #{{{ apply filtering and timeshift along t1
            if verbose: print "the starting value of t1 is",signal_slice.getaxis('t1')[0]
            if zerofill:
                signal_slice.ft('t1',shift = True,pad = 512)
            else:
                signal_slice.ft('t1',shift = True)
            if SW[0] is not None:
                signal_slice = signal_slice['t1':SW[0]]
            if verbose: print 'applied a t1 timeshift: t_start=',t_start/1e-9,'manual=',t1_timeshift/1e-9
            signal_slice.ift('t1')
            signal_slice.setaxis('t1',lambda x: x+t1_timeshift)
            #}}}
        #{{{ check that the timing correction looks OK
        fl.next('check timing correction')
        forplot = signal_slice.copy()
        if not 'te' in dropped_labels:
            forplot.ft('t1') # ft along t1 so I can clear the startpoint, below
        forplot.ft_clear_startpoints('t2')
        if not 'te' in dropped_labels:
            forplot.ft_clear_startpoints('t1')
            forplot.ift('t1') # start at zero here
        forplot.ift('t2',shift = True) # get a centered echo here
        fl.image(forplot, interpolation = 'bicubic')
        #ax = gca()
        #ax.axvline(x = 0,color = 'w',alpha = 0.25,linewidth = 1)
        #ax.axhline(y = 0,color = 'w',alpha = 0.25,linewidth = 1)
        #}}}
        #{{{ correct the zeroth order
        phase = signal_slice.copy().mean_all_but([None]).data
        assert isscalar(phase),' '.join(
                map(repr,["Error, phase is",phase])) # to make
                #    sure I'm not phasing anything differently
        phase /= abs(phase)
        signal_slice /= phase
        #}}}
        #}}}
        #{{{ various plots to check the phasing
        fl.next('phased data')
        fl.image(signal_slice, interpolation = 'bicubic')
        fl.next('cropped log -- to check phasing along $t_2$')
        fl.image(signal_slice.copy().cropped_log(), interpolation = 'bicubic')
        fl.grid()
        fl.next('and select real')
        fl.image(signal_slice.runcopy(real),interpolation = 'bicubic')
        #}}}
        #{{{ phase plot from 3553
        fl.next('SECSY:\nplot the phase at different $t_1$',legend = True)
        signal_slice.reorder('t2') # make t2 x
        max_abs = abs(signal_slice.data).flatten().max()
        for this_field in fields_toplot:
            for this_echotime in echos_toplot:
                forplot = signal_slice
                this_label = []
                if this_field is not None:
                    forplot = forplot['fields':(this_field)]
                    this_label.append('%.2f T'%(this_field))
                if this_echotime is not None:
                    forplot = forplot['t1':(this_echotime)]
                    this_label.append('%d ns'%(this_echotime/1e-9))
                fl.plot(forplot[lambda x:
                    abs(x) > phaseplot_thresholds[1]*max_abs
                    ].runcopy(angle)/pi,'.',alpha = 0.5,
                    markersize = 3,label = ', '.join(this_label),
                    human_units = False) # set human_units to false because
                #                          the units don't match
        fl.phaseplot_finalize()
        #}}}
        if not 'te' in dropped_labels:# there is a t1 dimension
            #{{{ check phasing along the indirect dimension
            signal_slice.ft('t1')
            fl.next('secsy mode')
            fl.image(signal_slice, interpolation = 'bicubic')
            #{{{ phase plot from 3553
            fl.next('secsy:\nplot the phase at different $F_2$',legend = True)
            def quarter_of_max(arg,axis):# returns a mask for half of maximum
                return arg > 0.25*arg.runcopy(max,axis)
            for this_field in fields_toplot:
                forplot = signal_slice
                this_label = []
                if this_field is not None:
                    forplot = forplot['fields':(this_field)]
                    this_label.append('%.2f T'%(this_field))
                sum_for_contiguous = abs(forplot).mean('t1').set_error(None)
                f_start,f_stop = sum_for_contiguous.contiguous(
                        quarter_of_max,'t2')[0,:]#  pull the first row, which is
                #                                   the largest block
                forplot.reorder('t1') # make t1 x
                max_abs = abs(forplot.data).flatten().max()
                this_label.append('')
                for this_freq in r_[f_start:f_stop:5j]:
                    this_label[-1] = '%d MHz'%(this_freq/1e6)
                    fl.plot(forplot['t2':(this_freq)][lambda x:
                        abs(x) > phaseplot_thresholds[0]*max_abs].runcopy(angle)/pi, '.',
                        alpha = 0.5, markersize = 3,
                        label = ', '.join(this_label))
                    #fl.plot(forplot['t2':(this_freq)].runcopy(angle)/pi,'.',alpha = 0.5,markersize = 3,label = '%d MHz'%(this_freq/1e6))
            fl.phaseplot_finalize()
            #}}}
            #}}}
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
        x[:] = field + x / 2.807e10
        #}}}
        fl.next(oned_plot_name, legend = True)
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
        xlabel(r'($B_0$ / $T$) + $\Delta f$ / ($2.807\times 10^{10}$ $\frac{Hz}{T}$)')
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
