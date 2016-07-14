# just make this a library of all NMR reading software
from .core import *
from .nmrfit import *
from .general_functions import process_kwargs
import re
import string
import os
import fornotebook
from scipy.io import loadmat
from .load_files import *

def OUTPUT_notebook():
    return True
#{{{ general, non file-format specific functions
def dbm_to_power(dbm,cavity_setup = 'newcnsi'):
    if cavity_setup == 'cnsi':
        attenuation = 30.0 # starting 6/20/11 this is actually what it uses, because I removed the extra -10, etc, in h5nmr.py
    if cavity_setup == 'newcnsi':
        attenuation = 35.16 # plus or minus a full dB --> this was measured
    elif cavity_setup == 'te102':
        attenuation = 40.0
    elif cavity_setup == 'dielectric':
        attenuation = 20.0
    return 1e-3*10.0**((dbm+attenuation)/10.0) #20 db for each atten
def power_to_dbm(power):
    return 10.0*log(power/1e-3)/log(10.0)-40 #20 db for each atten
#{{{ auto_steps
def show_powerlog(filename,threshold = -35,extra_t1_problem = False,figure_list = None,**kwargs):
    #{{{ generate the powers for the T1 series
    figure_list = fornotebook.figlistini(figure_list)
    #figure_list.text('\n\nCheck the $T_1$ powers:\n\n')
    figure_list.next('powerlog_raw')
    figure_list.next('checkpowerlog')
    figure_list.next('powerlog')
    t1_dbm,figure_list.figurelist = auto_steps(filename,
        threshold = threshold,
        first_figure = figure_list.figurelist,# note that auto_steps is old-style, so I pull the list out of the figlist class
        **kwargs)
    t1mask = ones(len(t1_dbm),dtype='bool')
    if extra_t1_problem == True:
        t1mask[-1] = 0 # this is the line you sometimes want to uncomment
    figure_list = check_autosteps(threshold,t1_dbm,figure_list = figure_list,mask = t1mask)
    figure_list.text('\n\nresulting array:'+lsafen(t1_dbm))
    #}}}
    return figure_list,t1_dbm,t1mask
def check_autosteps(threshold,values,figure_list = None,mask = None):
    if figure_list == None:
        figure_list = figlistl()
    figure_list.next('checkpowerlog')
    x = r_[r_[0:len(values)],r_[0:len(values)]+0.99]
    y = r_[values,values]
    s = argsort(x)
    plot(x[s],y[s],'b',linewidth=1)
    if mask is not None:
        if sum(mask) > 0:
            x = r_[r_[0:len(values)][~mask],r_[0:len(values)][~mask]+0.99,r_[0:len(values)][~mask]+0.991]
            y = r_[values[~mask],values[~mask],NaN*ones(shape(values[~mask]))]
            s = argsort(x)
            plot(x[s],y[s],'r',linewidth=2)
    expand_y()
    expand_x()
    ax = gca()
    ylims = list(ax.get_ylim())
    ylims[0] = threshold
    ax.set_ylim(ylims)
    gridandtick(ax)
    ylims = list(ax.get_ylim())
    ylims[0] = threshold
    ax.set_ylim(ylims)
    title('Interpretation:\n(red values are not used)')
    xlabel('experiment number')
    ylabel('power (dBm)')
    return figure_list
def auto_steps(filename,threshold = -35, upper_threshold = 30.0, t_minlength = 0.5*60,minstdev = 0.1,showplots = True, showdebug = False,t_start=0,t_stop=60*1000,tolerance = 2,t_maxlen = inf,return_lastspike = False,first_figure = None):
    r'Plot the raw power output in figure 1, and chop into different powers in plot 2'
    figurelist = figlistini_old(first_figure)
    v = loadmat(filename)
    p_ini = v['powerlist']
    t_ini = v['timelist']
    #{{{ plot the raw power log, then just pull out the times we're interested in
    if showplots:
        figurelist = nextfigure(figurelist,'powerlog_raw')
        plot(t_ini/60,p_ini)
        ax = gca()
        ylim = list(ax.get_ylim())
        ylim[0] = -35
        ax.set_ylim(ylim)
        gridandtick(ax)
        title(filename)
    mask = t_ini < t_stop
    minsteps = sum(t_ini<t_minlength)
    maxsteps = sum(t_ini<t_maxlen)
    mask = logical_and(mask,p_ini < upper_threshold)
    mask = logical_and(mask,t_ini > t_start)
    mask = logical_and(mask,isfinite(p_ini))
    mask = logical_and(mask,~(p_ini==-40.0))
    p_ini = p_ini[mask]
    t_ini = t_ini[mask]
    #}}}
    #a = whereblocks(logical_or(p_ini>threshold,~(isfinite(p_ini))))
    a = whereblocks(logical_or(p_ini>threshold,~(isfinite(p_ini))))
    a = a[argmax(map(len,a))]
    t = t_ini[a]
    p = p_ini[a]
    flattened = NaN*zeros(len(p)+minsteps)
    flattenedstd = NaN*zeros(len(p)+minsteps)
    track_curmean = NaN*zeros(len(p)+minsteps)
    track_threshold = NaN*zeros(len(p)+minsteps)
    maxp = max(p)
    minp = min(p)
    spikes = minp*ones(len(p)+minsteps)
    plotdev = zeros(len(p)+minsteps)
    plotstd = zeros(len(p)+minsteps)
    powerlist = []
    stdlist = []
    nextpos = 0
    while nextpos < len(t)-2:
        nextpos += 1 # just to skip a couple points so I don't grab the rise 
        figurelist = nextfigure(figurelist,'powerlog')
        # grab the stdev and average for the minimum number of steps
        blockstart = nextpos
        subset= p[nextpos:minsteps+nextpos]
        curmean = mean(subset[isfinite(subset)])
        stdev = std(subset[isfinite(subset)])
        if stdev < minstdev:
            stdev = minstdev
        # now that we have a decent stdev, just test to see that every step is less than the stdev
        track_curmean[nextpos] = curmean
        nextpos += minsteps
        #{{{ iterate over blocks
        while (nextpos < len(t)-1) and (abs(p[nextpos]-curmean)<tolerance*stdev or ~isfinite(p[nextpos])):
            subset= p[blockstart:nextpos]
            curmean = mean(subset[isfinite(subset)])
            stdev = std(subset[isfinite(subset)])
            if stdev < minstdev:
                stdev = minstdev
            if showdebug:
                plotdev[nextpos] = p[nextpos]-curmean
                plotstd[nextpos] = tolerance*stdev 
            track_curmean[nextpos] = curmean
            track_threshold[nextpos] = tolerance*stdev
            nextpos += 1
            ##{{{ if we get up to the maximum length, I need to break twice the maximum length into two jumps
            if t[nextpos]-t[blockstart] > t_maxlen:
                if t_minlength > t_maxlen:
                    raise CustomError('for auto_steps, minlength can\'t be greater than maxlen')
                #print 'DEBUG: maxlen triggered, nextpos=',nextpos
                biggestjump = blockstart+minsteps+1+argmax(abs(diff(t[blockstart+minsteps:blockstart+2*maxsteps]))) # find the biggest jump over the interval of two lengths
                if biggestjump-blockstart < 2*minsteps: # if I can't make two steps before the biggest jump
                    nextpos = blockstart + minsteps + 1 + argmax(abs(diff(t[blockstart+minsteps:blockstart+maxsteps])))
                    break
                if t[biggestjump]-t[blockstart] > t_maxlen: # otherwise, the biggest jump happens after another jump
                    #print 'DEBUG: greater than'
                    try:
                        nextpos = blockstart + minsteps + 1 + argmax(abs(diff(t[blockstart+minsteps:biggestjump-minsteps])))
                    except:
                        #figlisterr(figurelist,basename = pdfstring)
                        raise CustomError("I don't have room to make two minimum length steps of ",minsteps,"between the start of the block at ",blockstart," and the biggest jump at",biggestjump)
                    break
                else: # and sometimes the biggest jump happens before another jump, but longer than twice the minlen
                    #print 'DEBUG: less than'
                    nextpos = biggestjump
                    break
            ##}}}
        #}}}
        #{{{ need to recalculate the mean
        subset= p[blockstart:nextpos]
        curmean = mean(subset[isfinite(subset)])
        #}}}
        try:
            track_curmean[nextpos] = curmean
        except:
            raise CustomError('track\_curmean is only',len(track_curmean),'and t is only',len(t),'but nextpos is',nextpos)
        #uptohere = flattened[0:nextpos]
        #uptohere[uptohere==0] = cursum/curnum
        #flattened[nextpos-1] = cursum/curnum
        subset = flattened[0:nextpos] # subset of all the flattened stuff up to here
        subset[isnan(subset)] = curmean
        spikes[nextpos] = maxp # show a spike to see clearly where the data is divided
        try:
            lastspike = t[nextpos-minsteps]-t[-1]
        except:
            raise CustomError('len(spikes)=',len(spikes),'len(t)+minsteps=',len(t)+minsteps)
        subset = flattenedstd[0:nextpos] # subset of all the flattened stuff up to here
        subset[isnan(subset)] = stdev
        powerlist += [curmean]
        stdlist += [stdev]
        #flattened[nextpos-1] = curmean
    if showplots:
        plot(t/60,spikes[0:len(t)],'r')
        plot(t_ini/60,p_ini.flatten(),'x-b')
        plot(t/60,flattened[0:len(t)]+2*flattenedstd[0:len(t)],'g',alpha=0.5)
        plot(t/60,track_curmean[0:len(t)],'y',alpha=0.5)
        plot(t/60,track_curmean[0:len(t)]-track_threshold[0:len(t)],'k',alpha=0.5)
        plot(t/60,track_curmean[0:len(t)]+track_threshold[0:len(t)],'k',alpha=0.5)
        plot(t/60,flattened[0:len(t)]-2*flattenedstd[0:len(t)],'g',alpha=0.5)
        title(filename)
    if showdebug:
        if showplots:
            twinx()
        plot(t/60,plotdev,'k')
        plot(t/60,plotstd,'r')
        title('Power meter log')
    gridandtick(gca())
    retval = [array(powerlist)]
    if return_lastspike == True:
        retval += [lastspike]
    if first_figure is not None:
        retval += [figurelist]
    if len(retval) > 1:
        return tuple(retval)
    else:
        return retval[0]
#}}}
#{{{ error plot
def error_plot(*arg):
    width=6
    dpi=200
    fname = 'error_plot.png'
    fname = 'auto_figures/'+fname
    if grid:
        gridandtick(gca())
    savefig(fname,dpi=dpi)
    if figure:
        print r"""
        \begin{figure}[h]
        \end{figure}
        """
    print r'\includegraphics[width=%0.2fin]{%s}'%(width,fname)
    clf()
    raise CustomError(*arg)
#}}}
#}}}
#{{{ wrappers/generic functions to load acq and data files
b0 = r'$B_0$'
def show_acqu(vars):
    print '\\begin{verbatim}',vars.__repr__().replace(',','\n'),'\\end{verbatim}\n\n'
#{{{ load an nddata structure for a 2d set -- give the data needed to load
def bruker_det_rg(a):
    '''determine the actual voltage correction from the value of rg for a bruker NMR file'''
    return a
#}}}
#}}}
#{{{ lower level functions
#{{{ routines specific to Bruker
#{{{ Load T1 axis
def bruker_load_t1_axis(files):
    wait_time = []
    if type(files) is str:
        files = [files]
    #print 'DEBUG: trying to load files: ',files
    for thisfile in files:
        thisfile = dirformat(thisfile)
        thisfiletype = det_type(thisfile)
        if thisfiletype[0] == 'prospa':
            print 'need to copy over code from prospa'
        if thisfiletype[0] == 'bruker':
            wait_time += [bruker.load_vdlist(thisfile)]
        else:
            print 'couldn\'t determine thisfile type'
    return array(wait_time).flatten()
#}}}
#{{{ calculate decimation correction
def bruker_det_phcorr(v):
    if v['DIGMOD']==1:
        gdparray=array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[179,201,533,709,1097,1449,2225,2929,4481,5889,8993,11809,18017,23649,36065,47329,72161,94689,144353,189409,288737],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[184,219,384,602,852,1668,2292,3368,4616,6768,9264,13568,18560,27392,36992,55040,73856,110336,147584,220928,295040]])
        decimarray=array([2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024]) # the -1 is because this is an index, and copied from matlab code!!!
        dspfvs = v['DSPFVS']
        decim = v['DECIM']
        try:
            retval = gdparray[dspfvs,where(decimarray==decim)[0]]/2/decim
        except:
            print r'\begin{tiny}'
            print r'\begin{verbatim}'
            print CustomError('Problem returning',dspfvs,where(decimarray==decim)[0],'from gdparray of size',shape(gdparray),'because decimarray is of size',shape(decimarray))
            print r'\end{verbatim}'
            print r'\end{tiny}'
            retval = 0
        return retval
    else:
        return array([0])
#}}}
#{{{ load acqu
def bruker_match_line(line,number_re,string_re,array_re):
    m = number_re.match(line)
    if m:
        retval = (0,m.groups()[0],double(m.groups()[1]))
    else:
        m = string_re.match(line)
        if m:
            retstring = m.groups()[1]
            if retstring[-1]=='>':
                retstring = retstring[:-1]
            retval = (1,m.groups()[0],0,retstring)
        else:
            m = array_re.match(line)
            if m:
                retval = (2,m.groups()[0],(double(m.groups()[1]),double(m.groups()[2])),m.groups()[3])
            else:
                retval = (3,line)
    return retval
#{{{ bruker_load_acqu
def load_title(file):
    if det_type(file)[0] == 'bruker':
        return bruker.load_title(file)
    else:
        return ''
#}}}
#}}}
#}}}
#}}}
#{{{ higher level functions
# standard_epr moved to h5nmr
#{{{ optimize the zero-order phase the right way
def phaseopt(curve):
    curve = curve.copy()
    #{{{ find bestindex once
    phases = linspace(-pi/2,pi/2,100).reshape(1,-1) # this should work w/out altering the sign
    rotated_data = (curve.reshape(-1,1))*exp(-1j*phases)
    success = (real(rotated_data)**2).sum(axis=0)/((imag(rotated_data)**2).sum(axis=0)) #optimize the signal power to noise power
    bestindex = where(success==max(success))[0][0]
    #}}}
    #{{{ find the bestindex based on that 
    if (bestindex>phases.shape[1]-2):
        bestindex = phases.shape[1]-2
    phases = linspace(
            phases[0,bestindex-1],
            phases[0,bestindex+1],
            100).reshape(1,-1)
    rotated_data = (curve.reshape(-1,1))*exp(-1j*phases)
    success = (real(rotated_data)**2).sum(axis=0)/((imag(rotated_data)**2).sum(axis=0)) #optimize the signal power to noise power
    bestindex = where(success==max(success))[0][0]
    #}}}
    return exp(-1j*phases[0,bestindex])
#}}}
#{{{ fitting functions
#{{{ filter fitting 
def exp_fit(x,y):
    '''fit an exponential function, disregarding the first point'''
    testpoint = int(len(x)/10)
    startpoint = y[1:6].mean() #because of the filter, sometimes, the very first point can start at 0
    endpoint = y[-testpoint-6:-testpoint].mean()
    initial_slope = (y[testpoint]-startpoint)/(x[testpoint]-x[3])
    # initial slope = -startpoint/t1
    p_ini = [startpoint,-startpoint/initial_slope,endpoint]
    p_out,success = leastsq(exp_errfunc, p_ini,
            args = (x[1:],y[1:]),maxfev=3000,ftol = 1e-4, xtol = 1e-6)
    if (success<0) or (success>4):
        p_out = r_[1.,len(x)*4.,0]
        clf()
        plot(x,exp_fitfunc(p_ini,x))
        plot(x,exp_fitfunc(p_out,x))
        plot(x,y)
        legend(['initial guess','after fit','data'])
        error_plot(success,'There was an error attempting to fit')
    else:
        return p_out
def exp_fitfunc(p,x):
    return p[0]*exp(-x/p[1])+p[2]
def exp_errfunc(p,x,y):
    fit = exp_fitfunc(p,x)
    return fit-y
#}}}
#}}}
#{{{ routines for processing emax and t1 type curves
#{{{ rg_check --> check the receiver gain
def rg_check(file,expno,number_of_samples = 75,dynamic_range = 19,first_figure = None,show_complex = True):
    r'''show the fid, plotted vs. the max possible value to check the max of the receiver gain in figure 1
    in figure 2, plot in the complex plain to check digitization limit'''
    fl = figlistini_old(first_figure)
    # the following was load_emax, and I have not yet checked it
    data = load_file(file,expno,dimname = 'power',printinfo = False) # load the data
    listoffiles = format_listofexps([file,expno])
    params = load_acqu(listoffiles[0]) # load the parameters
    maxfloat = 3.4028234e38
    maxint = 2147483648
    max_for_this_dynamic_range= 2**(dynamic_range-1) # the -1 is the bit used for the sign
    data.reorder(['t2','power'])
    spec_type = det_type(listoffiles[0])[0]
    if spec_type == 'bruker':
        data.data *= params['RG'] # scale back up by RG, so we get raw numbers
    data.data /= max_for_this_dynamic_range # scale down by the max int number, so we test digitization
    nextfigure(fl,'rg1')
    if normalize == True:
        data['t2',:] /= data['t1',1:2]
    plot(data['t2',0:number_of_samples],'k')
    plot(data['t2',0:number_of_samples]*(-1j),'b',alpha=0.3)
    title('receiver gain upper value check')
    axis('tight')
    if show_complex:
        nextfigure(fl,'rg2')
        OLDplot(real(data['t2',0:number_of_samples].data),imag(data['t2',0:number_of_samples].data),'.',markersize=0.6)
        #polar(angle(data['t2',0:number_of_samples].data),abs(data['t2',0:number_of_samples].data),'.',markersize=0.1)
        title('receiver gain minimum value check')
        axis('tight')
    if first_figure == None:
        return
    else:
        return fl
#}}}
#{{{ integrate --> new integration function
def integrate(file,expno,
        integration_width=1e3,
        dimname = 'power',
        intpoints=None,
        show_image=True,
        filter_edge = 10,
        center_peak=True,
        use_baseline=True,
        plot_check_baseline=False,
        filter_direct = False,
        return_noise=True,
        show_integral = False,
        indiv_phase = False,
        scale_plot = False,
        peak_within = 1e3,
        bandpass = None,
        abs_image = False,
        max_drift = 1e3,
        first_figure = None,
        pdfstring = '',
        phnum = [],
        phchannel = [],
        offset_corr = 0):
    r'''new integration function, which replaces integrate_emax, and is used to integrate data, as for Emax and T1 curves'''
    #print lsafen("DEBUG: yes, integrate was called")
    figurelist = figlistini_old(first_figure)
    if type(plot_check_baseline) is bool:
        if plot_check_baseline:
            plot_check_baseline = 0 
        else:
            plot_check_baseline = -1
    phcycdims = ['phcyc%d'%j for j in range(1,len(phnum)+1)]
    data = load_file(file,expno,dimname = dimname, add_sizes = phnum,add_dims = phcycdims) # load the data
    #{{{ offset correction
    if type(offset_corr) is list:
        offset_corr = array(offset_corr)
    if type(offset_corr) is ndarray:
        data.data -= data['t2',offset_corr].copy().mean('t2').mean(dimname).data
    elif offset_corr > 0: # number of points to use for digitizer offset (zero glitch) correction
        offset_verbose = False
        if offset_verbose:
            nextfigure(figurelist,'beforecorr')
            image(abs(data))
            colorbar()
            print 'DEBUG: attempting offset correction by',offset_corr
        if array(offset_corr).dtype is dtype('float64'):
            data.data -= data['t2',lambda x: x > offset_corr].copy().mean('t2').mean(dimname).data
            if offset_verbose: print 'which is a double'
        else:
            data.data -= data['t2',-offset_corr:].copy().mean('t2').mean(dimname).data
        if offset_verbose:
            print '\n\n'
            nextfigure(figurelist,'aftercorr')
            image(abs(data))
            colorbar()
    #}}}
    # see_fid obsolete by rg_check
    # also remove all debug statements
    #print 'DEBUG: before phcyc, figlist is',lsafen(figurelist)
    data,figurelist = phcyc(data,names = phcycdims,selections = phchannel, show_plot = ['t2',(lambda x:abs(x)<peak_within)],first_figure = figurelist,pdfstring = pdfstring,bandpass = bandpass) # ft along t2, applying phase cycle where necessary
    data_shape = ndshape(data) # this is used to shape the output
    #{{{ abs w/ max SNR, so we can pick the peak
    data_abs = data.copy()
    data_abs['t2',(lambda x: abs(x)<peak_within)].data *= 0
    #{{{ apply the matched filter to maximize our SNR while picking the peak
    data_abs.ift('t2',shiftornot=True) # ft along t2
    filter = matched_filter(data_abs,'t2',decay_rate = 3)
    data_abs *= filter
    data_abs.ft('t2',shiftornot=True) # ft along t2
    #}}}
    data_abs = abs(data_abs)
    #}}}
    #{{{ generate topavg --> the index at the top of the average
    data_mean = data_abs.copy()
    data_mean.mean(dimname) # note that we are taking the mean over the abs here, which would not be great for the noise, but the center should still be in the right place
    data_mean.run(argmax,'t2') # put the index of the top peak there
    topavg = int32(data_mean.data)
    #}}}
    #{{{ generate center --> an array with the value at the center of each scan
    data_center = data_abs.copy() # since data_abs is currently not used, but I want to use it to do matched filtered integration, really need to make a separate variable here
    f = data_center.getaxis('t2')
    data_center['t2',abs(f-f[topavg])>max_drift] = 0# we need to keep the indeces in place, but don't want to pick anything too far out of the way
    test_drift_limit = False
    if test_drift_limit:
        plot(data_center.reorder(['t2',dimname]))
        return
    data_center_sum = data_center.copy()
    data_center_sum.sum_nopop('t2')
    data_center.data[:] /= data_center_sum.data # now the it sums to one so we can get a weighted average over the indeces
    f_ind = r_[0.0:double(size(f))].reshape(data_center.getaxisshape('t2')) # make a list of the indeces along the frequency axis
    #data_center.data[:] *= f_ind # multiply by that list
    #data_center.sum('t2') # sum, so that we would return the mean of the list if the spectrum were flat
    data_center.run(argmax,'t2')
    #center = int32(round(data_center.data))
    center = int32(array(map(round,data_center.data)))
    #}}}
    #{{{ if integration points unspec, pull out the integration width
    if intpoints==None:
        df = data.getaxis('t2').copy()
        df = df[1]-df[0]
        intpoints = floor(integration_width/(df))
    #}}}
    data_shape['t2'] = intpoints*2+1
    newdata = []
    newnoise = []
    center[center<intpoints] = intpoints # prevent a bug where the integration range exceeds the spectrum
    # plotting here removed, since this is done by phcyc
    #{{{baseline correction
    if use_baseline:
        #data.data += 0.1 #for debug baseline
        for j in range(0,ndshape(data)[dimname]):
            if plot_check_baseline == j: # if I've passed this index, show this baseline
                baseline_data = baseline_spectrum(data[dimname,j].copy(),center[j],intpoints,showplots = True) # call baseline_spectrum with center and intpoints, which should already be defined
                error_plot('Wanted to check baseline on scan ',j)
            else:
                baseline_data = baseline_spectrum(data[dimname,j].copy(),center[j],intpoints) # call baseline_spectrum with center and intpoints, which should already be defined
                data[dimname,j].data[:] -= baseline_data.data.flatten()
                if any(isnan(data[dimname,j].data[:])):
                    print 'isnan!!'
                if any(isinf(data[dimname,j].data[:])):
                    print 'isinf!!'
    #}}}
    #{{{ actually pull the integral points and the standard deviation of the noise
    #print 'DEBUG: center points are',center
    plot_noise = [] # array where we can put the noise for plotting
    other_info_retval = data.other_info # hack to compensate for the fact that other info isn't carried through the slice
    for j in range(0,data_shape[dimname]):
        newdata += [data[dimname,j,'t2',center[j]-intpoints:center[j]+intpoints+1]]
        if return_noise:
            #newnoise += [data[dimname,j,'t2',10:10+intpoints]]
            #{{{ grab intpoints to the left of the spectrum 
            temp = center[j]+r_[0:intpoints]-2*intpoints
            temp = temp[temp>0]
            #}}}
            #{{{ grab intpoints to the right of the spectrum 
            list_of_noise_indeces = center[j]+1+r_[0:intpoints]+intpoints
            list_of_noise_indeces = list_of_noise_indeces[list_of_noise_indeces<ndshape(data)['t2']]
            #}}}
            list_of_noise_indeces = int32(r_[temp,list_of_noise_indeces])
            #{{{ pull the noise data and calculate the standard deviation of the noise
            #print 'DEBUG: shape of data',ndshape(data),j,list_of_noise_indeces
            temp = data[dimname,j,'t2',list_of_noise_indeces]
            if show_image:
                plot_noise += [temp.copy()]
            temp.data = abs(temp.data-temp.data.mean()) # need to explicitly do the abs, since the data is complex
            temp.data **= 2
            temp.mean('t2') # I DO NOT take the sqrt here, because it's taken at the very end
            #}}}
            newnoise += [temp.copy()] # find the standard deviation of the noise which we have pulled --> it should be independent of the number of points that we're using
    newdata = concat(newdata,dimname)
    newdata.other_info = other_info_retval # hack to compensate for the fact that other info isn't carried through the slice
    if return_noise:
        newnoise = concat(newnoise,dimname)
        if show_image:
            plot_noise = concat(plot_noise,dimname)
    #}}}
    if show_image:
        plot_newdata = newdata.copy() # make a backup for plotting
    newdata.sum('t2') # integrate --> note that I converted this to a sum!
    #{{{ autophase
    if not indiv_phase:
        phaseoptval = phaseopt(newdata.data)
        newdata.data *= phaseoptval
        if show_image:
            plot_newdata.data *= phaseoptval
            if return_noise:
                plot_noise.data *= phaseoptval
        # do NOT rotate the noise data to be returned, since it's real!
    else:
        for j in range(0,len(newdata.data)):
            phcorr =  newdata[dimname,j]
            phcorr /= abs(phcorr)
            try:
                if show_image:
                    plot_newdata[dimname,j] *= phcorr
                newdata[dimname,j] *= phcorr
                if return_noise:
                    #newnoise[dimname,j] *= phcorr # again, don't rotate the real noise
                    if show_image:
                        plot_noise[dimname,j] *= phcorr
            except:
                print 'shape of newdatacopy',ndshape(newdatacopy)
                print 'shape of newdata',ndshape(newdata)
    #}}}
    #{{{ show what we're integrating
    if show_image:
        figurelist = nextfigure(figurelist,'intpeaks' + pdfstring)
        #print "DEBUG intpeaks figurelist =",figurelist,"gcf = ",gcf().number
        plot_newdata.reorder(['t2',dimname])
        def maybescale(x):
            if scale_plot:
                return x/newdata
            else:
                return x
        if scale_plot:
            plot_newdata /= newdata
            plot_noise /= newdata
        plot(plot_newdata,alpha=0.5)
        title('Peaks, zoomed in to integration region')
        if return_noise:
            plot_noise.reorder(['t2',dimname])
            plot_color_counter(0)
            plot(plot_noise['t2',0:intpoints],'-',alpha=0.1)
            plot_color_counter(0)
            plot(plot_noise['t2',intpoints:],'-',alpha=0.1)
        if show_integral:
            #{{{this does work to plot the integral
            plot_newdata.integrate_cumulative('t2') # this is apparently a function to do integral with all the correct bells and whistles
            #gridandtick(gca())
            ax = gca()
            myxlim = copy(ax.get_xlim())
            twinx()
            plot(plot_newdata,'k',alpha=0.1)
            gca().set_xlim(myxlim)
            #}}}
    #}}}
    #print lsafen("DEBUG: ready to return from integrate")
    #{{{ return integral and the noise of the 
    if return_noise:
        number_of_integral_points = 2.*intpoints+1. # the integral is a sum over this many points
        newdata.set_error(sqrt(newnoise.data.flatten()*number_of_integral_points))
    if first_figure == None:
        return newdata
    else:
        return newdata,figurelist
    #}}}
#}}}
#}}}
#{{{ load the data from a emax series based on input array
def print_info(filename,also = {}):
    'print the info for a file: to add other parameters, call with the "also" dictionary, with the description as a key, and the variable name or variable name, index number pair as the key'
    filetype,twod = det_type(filename)
    if filetype == 'bruker':
        v = bruker.load_acqu(dirformat(filename))
        f = open(dirformat(filename)+'pulseprogram','r')
        ppginfo = f.read()
        f.close()
        for m in re.finditer(r'\b([pd])([0-9]+)\b',ppginfo):
            also.update({m.group():[m.group(1).upper(),int(m.group(2))]})
        if OUTPUT_notebook():
            print r'\fn{%s}: sfo1:%0.5f aq:%0.3f swh:%0.3f ns:%d ds: %d rg:%0.1f d1:%0.1f p1:%0.2f pl1:%0.1f'%(filename,v['SFO1'],v['TD']/v['SW_h']/2.0,v['SW_h'],v['NS'],v['DS'],v['RG'],v['D'][1],v['P'][1],v['PL'][1])
            if len(also) > 0:
                for k,val in also.iteritems():
                    if type(val) is list or type(val) is tuple:
                        try:
                            print k,':',lsafe(v[val[0]][val[1]])
                        except:
                            print "(Can't find",k,val,map(type,val),"!)"
                    else:
                        print k,':',lsafe(v[val])
            data = fornotebook.save_data()
            pptvalue = v['SFO1']/data['current_frequency']
            if abs(pptvalue-data['current_ppt'])>1e-4:
                print ('WARNING! ppt value is off!! (',pptvalue,' ppt)')
            else:
                print '%0.4f $ppt$'%pptvalue
            print '\n\n'
#}}}
#{{{ deal with manual phase cycling
def phcyc(data,names=[],selections=[],remove_zeroglitch=None,show_plot = False,first_figure = None,pdfstring = '',bandpass = None):
    figurelist = figlistini_old(first_figure)
    data.ft('t2',shift=True)
    if (bandpass is not None):
        data = data['t2',lambda x: abs(x)<bandpass]
    if len(names)>0:
        data.ft(names)
        if remove_zeroglitch == None:
            remove_zeroglitch = True
    if remove_zeroglitch:
        index = argmin(abs(data.getaxis('t2')-0)) # need to incorporate this functionality into the abstracted function indexer
        indexlist = ['t2',index]
        for j,name in enumerate(names):
            indexlist += [name,0]
        data[tuple(indexlist)] = 0
    if show_plot:
        nextfigure(figurelist,'phcycchan' + pdfstring)
        allindexes = list(data.dimlabels)
        for name in names:
            if name in allindexes:
                allindexes.remove(name)
                allindexes = [name]+allindexes
        allindexes.remove('t2')
        data.reorder(allindexes+['t2'])
        if show_plot == True: # if passed something other than just "true", use that to subscript data
            image(data)
        else:
            image(data[tuple(show_plot)])
        titlestr = 'Raw data'
        if len(names)>0:
            titlestr += ' by phase cycle channel'
        title(titlestr)
    other_info_retval = data.other_info
    for j,name in enumerate(names):
        data = data[name,selections[j]]
    data.other_info = other_info_retval # this and the line that saves other_info_retval above are a hack because nddata doesn't do this correctly
    if first_figure == None:
        return data
    else:
        return data,figurelist
#}}}
#{{{ process_t1
def process_t1(file,expno,usebaseline = None,showimage = None,plotcheckbaseline = None,saturation = False,first_figure = None,pdfstring = '',t1_offset_corr = None,verbose = False,showlinear = False,**kwargs):
    #{{{ hack it, since it only actually takes a single file 
    file = format_listofexps([file,expno])
    if len(file) > 1: raise CustomError('I don\'t think this can handle more than one file at a time')
    expno = []
    #}}}
    figurelist = figlistini_old(first_figure)
    #{{{ legacy kwargs
    if showimage != None:
        show_image = showimage
    if usebaseline != None:
        use_baseline = usebaseline
    #}}}
    if type(file) is str:
        file = [file]
    titlestr = load_title(file[0])
    wait_time = load_t1_axis(file[0])
    if t1_offset_corr is not None:
        kwargs.update({'offset_corr':t1_offset_corr})
    integral,figurelist = integrate(file,expno,first_figure = figurelist,pdfstring = pdfstring,**kwargs)
    t1name = r'$t_1$'
    integral.rename('power',t1name)
    integral = t1curve(integral,fit_axis = t1name) # make this into an integral class, which fits along the dimension t1
    if ndshape(integral)[t1name] < len(wait_time):
        print '\n\nNote: ',t1name,'axis shorter than list of delays'
        wait_time = wait_time[0:ndshape(integral)[t1name]]
    integral.labels([t1name],[wait_time]) # before, I had to sort them manually, but now, I don't
    if verbose: print 'DEBUG wait times:',integral.getaxis(t1name)
    integral.sort(t1name)
    #{{{ finally, show the fit  
    figurelist = nextfigure(figurelist,'t1'+pdfstring)
    taxis = wait_time
    integral.data *= phaseopt(integral.data)
    plot(integral.runcopy(real),'ko')
    plot(integral.runcopy(imag),'yo')
    integral.makereal() # otherwise, it won't fit
    integral.fit()
    if verbose: print 'DEBUG: after fit, fit coeff is',integral.fit_coeff
    plot(integral.eval(300)) # evaluate the fit function on the axis taxis
    #{{{ for now, do not plot the modified versions
    #plot(taxis,t1_fitfunc(r_[p[0:2],p[2]*1.2],taxis),'y')
    #plot(taxis,t1_fitfunc(r_[p[0:2],p[2]*0.8],taxis),'y')
    #}}}
    ax = gca()
    text(0.5,0.75,integral.latex(),transform = ax.transAxes,size = 'x-large', horizontalalignment = 'center',color = 'r')
    title(titlestr)
    #}}}
    if showlinear:
        #{{{ and the straight line plot
        figurelist = nextfigure(figurelist,'t1straight'+pdfstring)
        #print 'linear data:',integral.linear().data
        plot(integral.linear(),'o')
        plot(integral.linear(taxis))
        #print ndshape(integral.linear())
        #}}}
    if first_figure == None:
        return integral # there is never a return_fit, since the fit is stored in the class itsself
    else:
        return integral,figurelist # there is never a return_fit, since the fit is stored in the class itsself
#}}}
#{{{ process cpmg 
def process_cpmg(file,dimname=''):
    r'this just takes an FT along T2 and pulls the signal at the frequency where it\'s maxed'
    data = load_file(file,dimname=dimname)
    data.ft('t2')
    findmax = abs(data)
    findmax.mean_all_but(['t2'])
    findmax = findmax.run(argmax,'t2').data
    data = data['t2',findmax]
    data.mean_all_but(['echo','t1'])
    data.data *= phaseopt(data.data) # I added this in, not sure why it was gone!
    return data
#}}}
#{{{regularization
def regularize1d(b,t,tau,alpha):
    # for Lx=b
        if size(b) != size(t):
            print "ERROR, size of b doesn't match size of t"
        tau = tau.reshape(1,-1)
        t = t.reshape(-1,1)
        L = exp(-t/tau)
        U,s,V = svd(L,full_matrices=0)
        rms = zeros(len(alpha),dtype='double')
        coeff = zeros((size(tau),size(alpha)),dtype='double')
        fit = zeros((size(t),size(alpha)),dtype='double')
        for j in range(0,len(alpha)):
            S = diag(s / (s**2 + alpha[j]**2))
            x = dot(
                    dot(
                        conj(transpose(V)),
                        dot(S,conj(transpose(U))))
                    ,b)# was b
            fit[:,j] = dot(L,x)
            try:
                coeff[:,j] = x.flatten()
            except:
                print 'shape(coeff)',shape(coeff),'shape(x)',shape(x)
                print 'first',shape(coeff[:,j]),'second',shape(x.reshape(-1,1))
                raise
            rms[j] = linalg.norm(fit[:,j]-b)
        return (coeff,fit,rms)
#}}}
#{{{ matched filter
def matched_filter(data,along_dim,decay_rate = 1,return_fit=False):
    r'''take ift'd data, and apply the matched filter to it'''
    #{{{ actually find the filter   
    data_abs = abs(data)
    timeaxis = data_abs.getaxis('t2')
    labels = list(data.dimlabels)
    labels.pop(labels.index(along_dim))
    for thisdim in labels:
        data_abs.mean(thisdim)
    p = exp_fit(timeaxis,data_abs.data)
    #}}}
    #{{{ actually apply the filter
    filter = ndshape(data)
    for thisdim in labels:
        filter.pop(thisdim)
    filter = filter.alloc()
    if (not return_fit): # don't mess with it if we want to check the fit
        p[1] /= decay_rate
    filter.data = exp_fitfunc(p,data.getaxis(along_dim).copy())
    #print 'for matched filter, the x axis is ',data.getaxis(along_dim).copy()
    if not return_fit:
        filter.data[:] -= filter.data.flatten()[-1] # drop so that end is at zero (since we have a noise baseline)
        filter.data[:] /= filter.data.flatten()[0] # normalize
    filter.labels(['t2'],[data_abs.getaxis('t2')])
    return filter
    #}}}
#}}}
def __baseline_gen_L(data):
    x = data.getaxis('t2').copy()
    x_norm = max(abs(x))
    x /= x_norm # normalize, otherwise we get ridiculously large higher order terms
    L = array([ones(shape(x)),x,x**2/2,x**3/6,x**4/24,x**5/120]).T
    #L = array([ones(shape(x)),x,x**2/2]).T
    return x,x_norm,L
def baseline_spectrum_peakpick(data,showplots=False,threshold=10,check_filter=False,set_errorlevel=False):
    #print 'diagnose: start baseline'
    data = data.copy()
    prelim_offset = mean(r_[data.data[0],data.data[1],data.data[-2],data.data[-1]])
    data.data[:] -= prelim_offset
    data.ift('t2',shiftornot=True)
    if check_filter:
        clf()
        #plot(abs(data))
        plot(matched_filter(data,'t2',return_fit=True))
        legend(['fits to']) #legend(['data','fits to'])
        twinx()
        plot(matched_filter(data,'t2'))
        legend(['filter with'])
        return
    #{{{ make the abs of a broadened spectrum   
    data_widefilt = data.copy() * matched_filter(data,'t2',decay_rate=10)
    data_widefilt.ft('t2',shiftornot=True)
    #data_widefilt.data *= phaseopt(data_widefilt.data)
    #data_widefilt.data = real(data_widefilt.data)
    #data_widefilt.data /= sign(data_widefilt.data[argmax(abs(data_widefilt.data))])
    data_widefilt.data -= (data_widefilt.data[0]+data_widefilt.data[-1])/2
    data_widefilt = abs(data_widefilt)
    mask = (data_widefilt.data<(data_widefilt.data.max()/threshold)) # generate the mask according to the threshold
    #}}}
    #{{{ ft the data
    data.ft('t2',shiftornot=True)
    #}}}
    if sum(mask)==0:
        erroronnopeak = False
        if erroronnopeak:
            legendstring = []
            plot(abs(data))
            legendstring += ['mf data']
            plot(data_widefilt)
            legendstring += ['wide filter']
            legend(legendstring)
            error_plot("Error -- not able to find any non-baseline data")
        else:
            print "Warning, fit entire spectrum to baseline.\n\n"
        mask = (data_widefilt.data<(data_widefilt.data.max()/1.5))
        mask[:] = True
    data_baseline = data.copy()
    data_baseline.data = data_baseline.data[mask]
    data_baseline.axis_coords[data_baseline.dimlabels.index('t2')] = data_baseline.getaxis('t2')[mask]
    legendstring = []
    x,x_norm,L = __baseline_gen_L(data_baseline) # return a normalized x axis, the value used to normalize it, and the array of normalized polynomials
    try:
        fit_coeff = dot(pinv(L,rcond=1e-5),data_baseline.data) # L * fit_coeff = data
    except:
        raise CustomError(maprep('problem inverting:',shape(L),shape(data_baseline.data)))
    #print 'diagnose: inverted'
    if set_errorlevel:
        if any(abs(dot(L,fit_coeff))>set_errorlevel):
            showplots = True
    if showplots:
        #plot(abs(data))
        #legendstring += ['data']
        plot(data_widefilt)
        legendstring += ['wide filter']
        plot(abs(data_baseline))
        legendstring += ['baseline portion']
        show_L = False
        if show_L:
            plot(x*x_norm,L)
            legendstring += map(
                    (lambda x:'L'+str(x)),
                    range(1,1+L.shape[1])
                    )
        plot(x*x_norm,abs(dot(L,fit_coeff)))
        legendstring += ['fit to baseline']
    x,x_norm,L = __baseline_gen_L(data)
    if showplots:
        plot(x*x_norm,dot(L,fit_coeff))
        legendstring += ['entire fit']
        legend(legendstring,'best')
    baseline_data = nddata(prelim_offset+dot(L,fit_coeff),[size(x),1],['t2','power'])
    baseline_data.labels(['t2'],[x])
    #print 'diagnose: shape of baseline ',ndshape(baseline_data)
    return baseline_data
def baseline_spectrum(data,center,points,showplots=False):
    data = data.copy()
    #{{{ here, I should define a mask in the same way I do in integrate, just inverted, so that I DON'T use the spectrum
    mask = bool8(zeros(shape(data.data)))
    mask[center-points:center+points+1] = True # this is an exact copy of the indeces used to create newdata in integrate
    mask = ~mask
    #}}}
    #{{{ splice out data_baseline --> the data which makes up the baseline
    data_baseline = data.copy()
    data_baseline.data = data_baseline.data[mask]
    data_baseline.axis_coords[data_baseline.dimlabels.index('t2')] = data_baseline.getaxis('t2')[mask]
    #}}}
    #{{{ perform the leastsquare fit
    x,x_norm,L = __baseline_gen_L(data_baseline) # return a normalized x axis, the value used to normalize it, and the array of normalized polynomials
    try:
        fit_coeff_Re = dot(pinv(L,rcond=1e-5),real(data_baseline.data)) # L * fit_coeff_Re = real(data)
    except:
        raise CustomError(maprep('problem inverting:',shape(L),shape(data_baseline.data)))
    try:
        fit_coeff_Im = dot(pinv(L,rcond=1e-5),imag(data_baseline.data)) # L * fit_coeff_Im = imag(data)
    except:
        raise CustomError(maprep('problem inverting:',shape(L),shape(data_baseline.data)))
    fit_coeff = fit_coeff_Re + 1j * fit_coeff_Im
    #}}}
    # here, i deleted set_errorlevel, since I can't remember what it does, so it must not be important
    #{{{ generate the matrices that span the full dataset and show the plots with all the info
    if showplots:
        clf() # this is in case I'm not running in the notebook, and want a decent error plot
        legendstring = []
        plot(abs(data),alpha=0.5)
        legendstring += ['data']
        plot(abs(data_baseline))
        legendstring += ['baseline portion']
        show_L = False
        if show_L:
            plot(x*x_norm,L)
            legendstring += map(
                    (lambda x:'L'+str(x)),
                    range(1,1+L.shape[1])
                    )
        plot(x * x_norm,abs(dot(L,fit_coeff)))
        legendstring += ['fit to baseline']
    x,x_norm,L = __baseline_gen_L(data)
    if showplots:
        plot(x*x_norm,dot(L,fit_coeff))
        legendstring += ['entire fit']
        data_forplot = data.copy()
    #}}}
    #{{{ generate and return the baseline curve
    baseline_data = nddata(dot(L,fit_coeff),[size(x),1],['t2','power'])
    baseline_data.labels(['t2'],[x])
    #}}}
    #{{{ show what the resulting data should look like
    if showplots:
        data.data[:] -= baseline_data.data.flatten()
        plot(abs(data))
        legendstring += ['baseline corrected data']
        legend(legendstring,'best')
    #}}}
    return baseline_data

#{{{ plot_noise
def plot_noise(path,j,calibration,mask_start,mask_stop,rgmin=0,k_B = None,smoothing = False, both = False, T = 293.0,plottype = 'semilogy',retplot = False):
    '''plot noise scan as resistance'''
    filename = r'%s%d'%(path,j)
    try:
        data = load_file(filename,calibration=calibration)
    except:
        raise CustomError('error loading file'+filename)
    k_B = 1.3806504e-23
    data.ft('t2',shift = True)
    newt2 = r'F2 / $Hz$'
    data.rename('t2',newt2)
    v = bruker.load_acqu(r'%s%d/'%(path,j))
    dw = 1/v['SW_h']
    dwov = dw/v['DECIM']
    rg = v['RG']
    de = v['DE']
    aq = v['TD']*dw
    if rg>rgmin:
        plotdata = abs(data)
        plotdata.data **= 2
        johnson_factor = 4.0*k_B*T
        plotdata.data /= (aq*johnson_factor)
        t = data.getaxis(newt2)
        mask = logical_and(t>mask_start,
            t<mask_stop)
        try:
            avg = plotdata.data[mask].mean() 
        except IndexError:
            raise CustomError('error trying to mask for the average because the mask is',mask,'of shape',shape(mask),'for shape(plotdata)=',shape(plotdata.data))
        retval = []
        if both or not smoothing:
            pval = plot(plotdata,'-',alpha=0.5,plottype = plottype)
            retval += ['%d: '%j+bruker.load_title(r'%s%d'%(path,j))+'$t_{dw}$ %0.1f $t_{dwov}$ %0.1f RG %d, DE %0.2f, mean %0.1f'%(dw*1e6,dwov*1e6,rg,de,avg)]
            axis('tight')
        if smoothing:
            # begin convolution
            originalt = plotdata.getaxis(newt2).copy()
            plotdata.ft(newt2,shift = True)
            sigma = smoothing
            siginv = 0.5*sigma**2 # here, sigma is given in the original units (i.e. what we're convolving)
            t = plotdata.getaxis(newt2)
            g = exp(-siginv*t.copy()**2) # we use unnormalized kernel (1 at 0), which is not what I thought!
            plotdata.data *= g
            plotdata.ift(newt2,shift = True)
            t = plotdata.getaxis(newt2).copy()
            t[:] = originalt
            # end convolution
            pval = plot(plotdata,'-',alpha=0.5,plottype = plottype)
            retval += ['%d: '%j+bruker.load_title(r'%s%d'%(path,j))+' $t_{dwov}$ %0.1f RG %d, DE %0.2f, mean %0.1f'%(dwov*1e6,rg,de,avg)]
            axis('tight')
        if retplot:
            return pval,retval
        else:
            return retval
    else:
        return []
#}}}
#}}}
