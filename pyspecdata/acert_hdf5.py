from .matlablike import *
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from .nmr import phaseopt
from sympy import var
import tables
import re
import os

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
def show_coherence_pathway(ax1,time_points,coherence_levels,yscale,gridpos,linewidth = 1):
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
            redlines[-1] += [(x[j],y[j]*yscale+gridpos)]
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
        selected_pulses = None # if given, this is a dictionary like {'phcyc3':(-1,2)}, where the second element of the tuple is the number of phase cycle steps selected
        ):
    pulseaxes = [j for j in plotdata.dimlabels if j[:5]=='phcyc']
    netshape = ones(len(plotdata.data.shape),dtype = int)
    #for thisaxis in pulseaxes:
    #    x = plotdata.getaxis(thisaxis)
    #    axis_size = plotdata.shape[plotdata.axn(thisaxis)]
    #    #x[x>axis_size/2] = x[x>2]-4
    #{{{ now, figure out the shape of the matrix where I knock out everything but the phcyc dims -- seems there should be an nddata-like way to do this
    for j in pulseaxes:
        thisindex = plotdata.dimlabels.index(j)
        netshape[thisindex] = plotdata.data.shape[thisindex]
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
    #print "coherence changes:",coherence_changes
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
    h5 = tables.openFile(filename)
    #{{{ set up the complex number the hard way, for good form
    data = empty_like(h5.root.experiment._v_children['data.r'],dtype = complex128)
    data_introspect = data.view([('r',double),('i',double)])
    for i_or_r in ['i','r']:
        data_introspect[i_or_r] = h5.root.experiment._v_children['data.'+i_or_r]
    #}}}
    if len(data.shape) == 1:
        if use_sweep:
            data = nddata(data,data.size,['current']).labels('current',array(h5.root.experiment._v_children['sweep_currents']))
        else:
            data = nddata(data,data.size,['field']).labels('field',array(h5.root.experiment._v_children['fields']))
    elif len(data.shape) == 2:
        if use_sweep:
            data = nddata(data,data.shape,['repeats','current']).labels('current',array(h5.root.experiment._v_children['sweep_currents']))
        else:
            data = nddata(data,data.shape,['repeats','field']).labels('field',array(h5.root.experiment._v_children['fields']))
    h5.close()
    return data
def postproc_blank(data):
    return data
def postproc_generic(data):
    data = automagical_phasecycle(data)
    return data
def postproc_eldor_3d(data):
    if data.get_prop('phasecycle') is None:
        print "Warning!  There is no phase cycling information -- you should really fix this"
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
        print "Warning!  There is no phase cycling information -- you should really fix this"
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
    data.ift('t2',shift = True)
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
def search_freed_file(searchstring,exptype,directory = getDATADIR()):
    directory = dirformat(directory + exptype)
    files = re.findall('.*' + searchstring + '.*','\n'.join(os.listdir(directory)))
    if len(files) == 0:
        raise ValueError("I can't find a file matching the regular expression"+searchstring)
    else:
        if len(files) > 1:
            print 'found:',files,'and opening last\n\n'
    return directory + files[-1]
def find_file(searchstring,
            base_directory = getDATADIR(),
            subdirectory =  'b1_fid/',
            postproc = None,
            print_result = True,
            verbose = False,
            indirect_dimlabels = None,# in case "dimlabels" is not set properly, I can manually pass the value of "indirect_dimlabels"
            **kwargs):
    r'''find the file (assumed to be an h5 format file) given by the regular
    expression "searchstring" inside the directory "directory", and postprocess
    with the function "postproc," which is fed the nddata data and the source
    text as arguments, as well as the keyword arguments "kwargs".
    If "postproc" is not set, it's chosen based on the value of
    (h5 root).experiment.description['class']'''
    #{{{ actually find one file and load it into the h5 object
    directory = base_directory + subdirectory
    filename_re = re.compile('.*' + searchstring + '.*')
    if os.path.isdir(directory):
        files = filename_re.findall('\n'.join(os.listdir(directory)))
    else:
        raise RuntimeError("I can't find the directory:\n%s\nin order to get a file that matches:\n%s"%(directory,searchstring))
    if len(files) == 0:
        raise ValueError("I can't find a file matching the regular expression "+searchstring+" in "+directory)
    if verbose:
        if len(files) > 1:
            print 'found:',files,'and opening last'
        elif print_result:
            obsn("found only one file, and loading it:"+repr(files))
    filename = directory + files[-1]
    h5 = tables.openFile(filename)
    #}}}
    #{{{ set up the complex number the hard way, for good form
    data = empty_like(h5.root.experiment._v_children['data.r'],dtype = complex128)
    data_introspect = data.view([('r',double),('i',double)])
    for i_or_r in ['i','r']:
        data_introspect[i_or_r] = h5.root.experiment._v_children['data.'+i_or_r]
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
    #{{{ now, pull the dimlabels that we want
    if 'dimlabels' in h5.root.experiment._v_attrs.__dict__.keys():
        dimlabels = h5.root.experiment._v_attrs.__dict__['dimlabels'].tolist()
    else:
        print "You didn't set dimlabels in your pulse program -- the data is ambiguous and I'm taking my best guess!!!"
    if indirect_dimlabels is not None:# in case it's not stored correctly in the file and I need to hard-code it
        idx = dimlabels.index('indirect')
        dimlabels = dimlabels[:idx] + indirect_dimlabels + dimlabels[idx+1:]
    #}}}
    #{{{ put all further information into the nddata in a way such that it can be retrieved with "get_prop"
    data.set_prop(dict([(k,v) for k,v in h5.root.experiment._v_attrs.__dict__.iteritems() if k[0] != '_' and k[0] != 'dimlabels']))
    data.set_prop(dict([('execution_'+k,v) for k,v in h5.root.experiment.execution._v_attrs.__dict__.iteritems() if k[0] != '_']))
    data.set_prop(dict([('description_'+k,v) for k,v in h5.root.experiment.description._v_attrs.__dict__.iteritems() if k[0] != '_']))
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
    #}}}
    #{{{ assign all the dimensions
    indirect_names_in_h5 = [k for k,v in h5.root.experiment._v_children.iteritems() if k not in ['data.i','data.r','bin_switch_times','fields','sweep_currents'] and isinstance(v,tables.array.Array)]
    if verbose: print "dimlabels are",dimlabels
    expected_dimlabels = set(dimlabels) - {'phcyc','t2','bin'} # we expect to find these axes stored in the HDF5
    if verbose: print "expected dimlabels",expected_dimlabels
    common_dimensions = expected_dimlabels & set(indirect_names_in_h5)
    #{{{ deal with the case where dimlabels and the stored axes don't match up correctly 
    if verbose: print "common dimensions",common_dimensions
    unlabeled_indirect = False
    if len(indirect_names_in_h5) == 0:
        print lsafen('Warning!! The indirect dimension is not labeled!\nDimensions of the data:%s\nfilename: %s'%(repr(ndshape(data)),filename))
        unlabeled_indirect = True
    if len(common_dimensions) != len(indirect_names_in_h5):
        h5_but_not_dimlabels = common_dimensions ^ set(indirect_names_in_h5)
        dimlabels_but_not_h5 = common_dimensions ^ expected_dimlabels
        if len(h5_but_not_dimlabels) == 1 and len(dimlabels_but_not_h5) == 1:
            real_name = h5_but_not_dimlabels.pop()
            dimlabels_name = dimlabels_but_not_h5.pop()
            print lsafen("Warning! The dimlabels set in the HDF5 file have a dimension called %s and the arrays given in the HDF5 file have one called %s -- I'm assuming they're the same and using the name from the HDF5 file."%(dimlabels_name,real_name))
            common_dimensions |= {real_name} # so the axis gets labeled below
            dimlabels[dimlabels.index(dimlabels_name)] = real_name
        elif dimlabels_but_not_h5 == {'indirect'}:# I need to split up the indirect dimension
            raise ValueError("You want me to split the indirect dimension, but I don't know how!:\n\tDimensions in the file: %s\n\tDimensions listed in dimlabels %s\n\tDimensions of the data: %s"%(repr(indirect_names_in_h5),repr(dimlabels),repr(ndshape(data))))
        else:
            raise ValueError("This file has dimensions that I can't just automatically sort!:\n\tDimensions in the file: %s\n\tDimensions listed in dimlabels %s\n\tDimensions of the data: %s\n\tDimensions in hdf5 but not dimlabels: %s\n\tDimensions in dimlabels but not HDF5: %s"%(repr(indirect_names_in_h5),repr(dimlabels),repr(ndshape(data)),repr(h5_but_not_dimlabels),repr(dimlabels_but_not_h5)))
    #}}}
    if verbose: print "now dimlabels are",dimlabels,"and indirect_names_in_h5 is",indirect_names_in_h5
    if not unlabeled_indirect:
        if 'indirect' in data.dimlabels:# note that if there was only a single indirect dimension, it's already been used up
            #{{{ chunk the "indirect" dimension up appropriately, and assign the dimensions
            chunk_dict = dict([(j,len(h5.root.experiment._v_children[j])) for j in dimlabels if j not in ['bin','phcyc','t2']])
            if verbose: print ndshape(data)
            if verbose: print chunk_dict
            data.chunk('indirect',chunk_dict)
        for this_axis in common_dimensions:
            axis_data = array(h5.root.experiment._v_children[this_axis])
            data.labels(this_axis,axis_data)
        #}}}
    #}}}
    #{{{ label the bin dimension, and chunk it up appropriately
    if 'bin' in data.dimlabels:
        #{{{ assign the labels for the bin dimension as a structured array
        forstruct = []
        forstruct_names = []
        for thisaxis in [x for x in 'bin_switch_times','fields','sweep_currents' if x in h5.root.experiment._v_children.keys()]:
            forstruct.append((h5.root.experiment._v_children['bin_switch_times'])[0])
            forstruct_names.append(thisaxis)
        x = make_rec(forstruct,forstruct_names,zeros_like = h5.root.experiment._v_children['bin_switch_times'].shape)
        print "detected a 'bin' dimension, and associating it with dtype",x.dtype
        for thisaxis in forstruct_names:
            x[thisaxis] = array(h5.root.experiment._v_children[thisaxis])
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
    difference.ft('t2',shift = True)
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
