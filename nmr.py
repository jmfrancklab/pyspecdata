# just make this a library of all NMR reading software
from matlablike import *
from nmrfit import *
import re
import string
import struct
import os
import fornotebook
from scipy.io import loadmat

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
def auto_steps(filename,threshold = -35, upper_threshold = 5, t_minlength = 0.5*60,minstdev = 0.1,showplots = True, showdebug = False,t_start=0,t_stop=60*1000,tolerance = 2,t_maxlen = inf,return_lastspike = False,first_figure = None):
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
                        figlisterr(figurelist,basename = pdfstring)
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
#{{{ load the pulse sequence parameters
def load_acqu(filename,whichdim='',return_s = None):
    filename = dirformat(filename)
    if det_type(filename)[0] == 'bruker':
        if return_s is not None:
            return bruker_load_acqu(filename,whichdim=whichdim,return_s = return_s)
        else:
            return bruker_load_acqu(filename,whichdim=whichdim)
    elif det_type(filename)[0] == 'prospa':
        if det_type(filename)[1] == 't1_sub':
            filename = dirformat(filename)
            return prospa_load_acqu(filename+'../')
        else:
            return prospa_load_acqu(filename)
    else:
        raise CustomError(det_type(filename),'is not yet supported')
#}}}
#{{{ is this (bruker or prospa, 1d or nd)?
def det_type(filename):
    filetype = None
    #{{{ WinEPR
    if os.path.exists(filename+'.spc'):
        return ('winepr',True)
    #}}}
    else:
        filename = dirformat(filename)
        files_in_dir = os.listdir(filename)
        #{{{ Bruker 2D
        if os.path.exists(filename+'ser'):
            return ('bruker',True)
        #}}}
        #{{{ Prospa generic 2D
        elif os.path.exists(filename+'data.2d'):
            return ('prospa',True)
        #}}}
        #{{{ specific Prospa formats
        elif any(map((lambda x:'Delay' in x),files_in_dir)):
            return ('prospa','t1')
        elif os.path.exists(filename+'acqu.par'):
            return ('prospa',False)
        elif os.path.exists(filename+'../acqu.par'):
            return ('prospa','t1_sub')
        #}}}
        #{{{ Bruker 1D
        elif os.path.exists(filename+'acqus'):
            return ('bruker',False)
        #}}}
        else:
            raise CustomError('WARNING! unidentified file type '+filename)
#}}}
#{{{ load an nddata structure for a 2d set -- give the data needed to load
def format_listofexps(args):
    'aux function, just used to decode (filename,listofexpnos) vs. (listoffilenames) in an arbitrary way'
    if type(args[0]) is str: # even if it's just a string, make it a list
        args[0] = [args[0]]
    if len(args) > 1 and ((not isscalar(args[1])) and len(args[1]) == 0): args.pop(1)
    if len(args) > 1:
        if len(args) > 2: raise CustomError('wrong number of args!')
        if isscalar(args[1]): args[1] = [args[1]] # if the second argument is a single file number, make it into a list
        if len(args[0]) > 1: raise CustomError("you can't have both the filename and the expnos be longer than 1")
        filenames = [dirformat(args[0][0]) + '%d'%x for x in args[1]]
    else:
        filenames = args[0]
    return filenames
def load_file(*args,**kwargs):
    'load a file or series of files; load as load_file(filename) or load_file(dirname,expnos)'
    args = list(args)
    #{{{ manually do kwargs
    dimname = '' 
    if 'dimname' in kwargs.keys(): dimname = kwargs['dimname']
    calibration = 1.0
    if 'calibration' in kwargs.keys(): calibration = kwargs['calibration']
    add_sizes = []
    if 'add_sizes' in kwargs.keys(): add_sizes = kwargs['add_sizes']
    add_dims = []
    if 'add_dims' in kwargs.keys(): add_dims = kwargs['add_dims']
    #}}}
    filenames = format_listofexps(args)
    #{{{load all the data into a list
    data = [load_indiv_file(filenames[0],dimname=dimname,add_sizes = add_sizes,add_dims = add_dims)]
    for filename in filenames[1:]:
        data += [load_indiv_file(filename,dimname=dimname,add_sizes = add_sizes,add_dims = add_dims)]
    #}}}
    # for the following, I used to have a condition, but this is incompatible with the pop statement at the end
    newdata = concat(data,dimname) # allocate the size of the indirect array
    #print 'DEBUG concatenated list = ',data
    newdata_shape = ndshape(newdata)
    if all(map((lambda x:det_type(x)[0]=='prospa'),filenames)):
        if hasattr(data[0],'want_to_prospa_decim_correct'):
            if data[0].want_to_prospa_decim_correct is True:
                newdata = prospa_decim_correct(newdata)
    #print 'DEBUG concatenated list before pop = ',data
    if newdata_shape[dimname]==1:
        newdata.popdim(dimname)
    return newdata*calibration
def bruker_det_rg(a):
    '''determine the actual voltage correction from the value of rg for a bruker NMR file'''
    return a
def load_indiv_file(filename,dimname='',return_acq=False,add_sizes=[],add_dims=[]):
    filetype,twod = det_type(filename)
    #{{{ ESR spectra
    if filetype == 'winepr':
        fp = open(filename+'.spc','rb')
        data = fp.read()
        data = array(
                struct.unpack('<%df'%(len(data)/4),data),
                dtype='double')
        v = winepr_load_acqu(filename)
        xpoints = v['RES']
        rg = v['RRG']
        data /= rg
        modulation = v['RMA']
        #data /= modulation
        try:
            data /= v['JNS'] # divide by number of scans
        except:
            pass
        data /= v['MP'] # divide by power
        ypoints = len(data)/xpoints
        if ypoints>1:
            if ypoints != v['REY']:
                raise CustomError('I thought REY was the indirect dim, guess not')
            if dimname=='':
                dimname = v['JEY']
            data = nddata(data,[ypoints,xpoints],[dimname,b0])
        else:
            data = nddata(data,[xpoints],[b0])
        xlabels = linspace(v['HCF']-v['HSW']/2.,v['HCF']+v['HSW']/2.,xpoints)
        if len(data.dimlabels)>1:
            data.labels([dimname,b0],[linspace(0,1,ypoints),xlabels])
            data.reorder([b0,dimname])
        else:
            data.labels([b0],[xlabels])
        return data
    #}}}
    filename = dirformat(filename)
    #{{{ Bruker 2D
    if twod and filetype == 'bruker':
        v = bruker_load_acqu(filename)
        v2 = bruker_load_acqu(filename,whichdim='2')
        td2 = int(v['TD'])
        rg = bruker_det_rg(float(v['RG']))
        td1 = int(v2['TD'])
        td2_zf = int(ceil(td2/256.)*256) # round up to 256 points, which is how it's stored
        fp = open(filename+'ser','rb')
        data = fp.read()
        data = array(struct.unpack('>%di'%(len(data)/4),data),
                dtype='complex128')
        data = data[0::2]+1j*data[1::2]
        data /= rg
        mydimsizes = [td1,td2_zf/2]
        mydimnames = [dimname]+['t2']
        #print 'DEBUG: data going to nddata =',data
        try:
            data = nddata(data,mydimsizes,mydimnames)
        except:
            size_it_should_be = array(mydimsizes).prod()
            if size_it_should_be > len(data):
                zero_filled_data = zeros(size_it_should_be)
                zero_filled_data[0:len(data)] = data
                data = nddata(zero_filled_data,mydimsizes,mydimnames)
            else:
                new_guess = len(data)/(td2_zf/2)
                print lsafen("WARNING!, chopping the length of the data to fit the specified td1 of ",td1,"points!\n(specified ",zip(mydimnames,mydimsizes),' td2_zf=%d)'%td2_zf)
                #td2_zf_new = 2**ceil(log(td2)/log(2))
                #mydimsizes[1] = td2_zf_new
                #size_it_might_be = array(mydimsizes).prod()
                #print "maybe this works:",size_it_might_be == len(data)
                data = data[0:size_it_should_be]
                data = nddata(data,mydimsizes,mydimnames)
                #raise CustomError("found td1=",td1,"for",filename,"which I don't think is right, because the product of the dimensions",zip(mydimnames,mydimsizes),'=',size_it_should_be,'does not equal the length of the data',len(data),'I think that it should be',len(data)/(td2_zf/2))
        #print 'DEBUG: data straight from nddata =',data
        data = data['t2',0:td2/2] # now, chop out their zero filling
        t2axis = 1./v['SW_h']*r_[1:td2/2+1]
        t1axis = r_[0:td1]
        mylabels = [t1axis]+[t2axis]
        data.labels(mydimnames,mylabels)
        shiftpoints = int(bruker_det_phcorr(v)) # use the canned routine to calculate the first order phase shift
        data.circshift('t2',shiftpoints)
        data.set_units('t2','s')
        data.set_units('digital')
        data.other_info['title'] = bruker_load_title(filename)
        #print 'DEBUG 2: data from bruker file =',data
        #}}}
        #{{{ Prospa 2D
    elif twod and filetype == 'prospa':
        if twod == 't1_sub':
            v = prospa_load_acqu(filename+'../') # if it's subdirectory format, the file comes from one directory up
            indirect_dim_len = [1]
            indirect_dim_name = [dimname]
            dimshere = 1
        else:
            v = prospa_load_acqu(filename)
            indirect_dim_name = []
            indirect_dim_len = []
            dimshere = 2
        taxis = linspace(0,1,v['nrPnts'])*v['acqTime']/1e3 # this is the t2 dimension, and so is always true
        data = prospa_load_datafile(filename,dims=dimshere)/v['nrScans']#added this 2/20/13 to allow automatic signal averaging
        #{{{ Prospa CPMG
        if v['experiment'].find('cpmg') > -1:
            data = nddata(data,indirect_dim_len+[v['nrEchoes'],v['nrPnts']],indirect_dim_name+['echo','t2'])
            echotime = (r_[0:v['nrEchoes']]+0.5)*v['echoTime']/1e6
            data.labels(indirect_dim_name+['echo','t2'],indirect_dim_len+[echotime,taxis])
            data.want_to_prospa_decim_correct = False
        #}}}
        #{{{ Prospa where 1D subscan is not CPMG
        else:
            data = nddata(data,indirect_dim_len+[v['nrPnts']],indirect_dim_name+['t2'])
            data.labels([dimname,'t2'],[r_[1],taxis])
            data.want_to_prospa_decim_correct = True
        #}}}
        #}}}
        #{{{ bruker 1D
    else:
        if filetype == 'bruker':
            v = bruker_load_acqu(filename)
            td2 = int(v['TD'])
            td1 = 1
            td2_zf = int(ceil(td2/256.)*256) # round up to 256 points, which is how it's stored
            fp = open(filename+'fid','rb')
            data = fp.read()
            data = array(
                    struct.unpack('>%di'%(len(data)/4),data),
                    dtype='complex128')
            data = data[0::2]+1j*data[1::2]
            rg = bruker_det_rg(v['RG'])
            data /= rg
            data = nddata(data,[td1,td2_zf/2],[dimname,'t2'])
            data = data['t2',0:td2/2] # now, chop out their zero filling
            t2axis = 1./v['SW_h']*r_[1:td2/2+1]
            t1axis = r_[1]
            data.labels([dimname,'t2'],[t1axis,t2axis])
            shiftpoints = int(bruker_det_phcorr(v)) # use the canned routine to calculate the second order phase shift
            #print 'shiftpoints = ',shiftpoints
            data.circshift('t2',shiftpoints)
            # finally, I will probably need to add in the first order phase shift for the decimation --> just translate this
            data.other_info['title'] = bruker_load_title(filename)
        #}}}
        #{{{ prospa 1d
        elif filetype == 'prospa':
            v = prospa_load_acqu(filename)
            data = prospa_load_datafile(filename,dims=1)/v['nrScans']#added this 2/20/13 to allow automatic signal averaging
            data = nddata(data,[v['nrPnts']],['t2'])
            taxis = linspace(0,1,v['nrPnts'])*v['acqTime']/1e3
            data.labels(['t2'],[taxis])
        #}}}
        else:
            raise CustomError("can't load this file type $\\rightarrow$ \\verb+%s+"%filename)

    #{{{ return, and if necessary, reorganize
    if len(add_sizes)>0:
        data.labels([dimname],[[]]) # remove the axis, so we can reshape
        #print 'DEBUG: data before chunk = ',data
        data.chunkoff(dimname,add_dims,add_sizes)
        #print 'DEBUG: data after chunk = ',data
        data.labels(add_dims,
                [r_[0:x] for x in add_sizes])
    if return_acq:
        return (data,v,v2)
    else:
        return data
    #}}}
#}}}
#{{{ t1 axis
def load_t1_axis(file):
    if det_type(file)[0] == 'bruker':
        return bruker_load_t1_axis(file)
    elif det_type(file)[0] == 'prospa':
        return prospa_t1_info(file)[1]
    else:
        raise CustomError('Trying to load T1 axis on a file of unrecognized format!')
#}}}
#}}}
#{{{ lower level functions
#{{{ routines specific to prospa
#{{{ load acqu
def prospa_decim_correct(data):
    #{{{ get rid of the finite rise time    
    data_abs = abs(data)
    otherdims = ndshape(data)
    otherdims.pop('t2')
    for indirect_dim_name in otherdims.dimlabels:
        data_abs.mean(indirect_dim_name)
    data_abs = data_abs.run(argmax,'t2')
    top = int(data_abs.data)
    data.circshift('t2',top)
    #}}}
    print 'Applied prospa decimation correction'
    return data
def prospa_load_acqu(file):
    file = dirformat(file)
    fp = open(file+'acqu.par')
    lines = fp.readlines()
    line_re = re.compile(r'([^ \t]+) *= *(.+)')
    vars = {}
    for j in range(0,len(lines)):
        lines[j] = string.rstrip(lines[j])
        m = line_re.match(lines[j])
        if m:
            exec 'temp = %s'%m.groups()[1]
            vars.update({m.groups()[0]:temp})
        else:
            print "error, acqu.par line not parsed: ",lines[j]
    fp.close()
    return vars
#}}}
#{{{ load_datafile
def prospa_load_datafile(file,dims=1):
    r'''load a prospa datafile into a flat array as a 1D file
    use dims=2 if it's a 2D file'''
    file = dirformat(file)
    if dims == 1:
        fp = open(file+'data.1d','rb')
    elif dims == 2:
        fp = open(file+'data.2d','rb')
    else:
        print 'ERROR: wrong number of dims'
    data = fp.read()
    data = array(struct.unpack('%df'%(len(data)/4),data))
    data = data[7:]
    # the following is junk!!!
    #elif precision=='b':
    #   data = array(struct.unpack('%db'%(len(data)/1),data))
    #   data = data[7*4:]
    #else:
    #   print 'error, precision wrong'
    data = reshape(data,(-1,2))
    data = data[:,0]+1j*data[:,1]
    fp.close()
    return data
#}}}
#{{{ load the data from a t1 file based on file names
def prospa_t1_info(file):
    file = dirformat(file)
    if det_type(file) == ('prospa','t1_sub'):
        file += '../'
    elif not det_type(file) == ('prospa','t1'):
        raise CustomError("You're trying to get prospa T1 info from a file that's not a Prospa T1 file!")
    files = [x for x in os.listdir(file) if os.path.isdir(dirformat(file)+x)]
    file_re = re.compile(r'([0-9]+)msDelay$')
    datafiles = []
    wait_time = []
    print 'DEBUG: prospa is searching for times in the file list',files
    for j in range(0,len(files)):
        m = file_re.match(files[j])
        if m:
            datafiles += [file+files[j]]
            wait_time += [int(m.groups()[0])]
    return datafiles,array(wait_time)*1e-3
#}}}
#}}}
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
            wait_time += [bruker_load_vdlist(thisfile)]
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
#{{{ bruker_load_vdlist
def bruker_load_vdlist(file):
    fp = open(file+'vdlist')
    lines = fp.readlines()
    lines = map(string.rstrip,lines)
    lines = map((lambda x: x.replace('m','e-3')),lines)
    lines = map((lambda x: x.replace('s','')),lines)
    lines = map((lambda x: x.replace('u','e-6')),lines)
    lines = map(double,lines)
    return array(lines)
#}}}
#{{{ bruker_load_acqu
def load_title(file):
    if det_type(file)[0] == 'bruker':
        return bruker_load_title(file)
    else:
        return ''
def bruker_load_title(file):
    file = dirformat(file)
    fp = open(file+'pdata/1/title')
    lines = fp.readlines()
    emptystring = '\r\n'
    while emptystring in lines:
        lines.pop(lines.index(emptystring))
    emptystring = '\n'
    while emptystring in lines:
        lines.pop(lines.index(emptystring))
    return ''.join(lines)
def bruker_load_acqu(file,whichdim='',return_s = True):
    if return_s:
        fp = open(file+'acqu'+whichdim+'s')# this is what I am initially doing, and what works with the matched filtering, etc, as is, but it's actually wrong
    else:
        fp = open(file+'acqu'+whichdim)# this is actually right, but doesn't work with the matched filtering, etc.
    lines = fp.readlines()
    vars = {}
    number_re = re.compile(r'##\$([_A-Za-z0-9]+) *= *([0-9\-\.]+)')
    string_re = re.compile(r'##\$([_A-Za-z0-9]+) *= *<(.*)')
    array_re = re.compile(r'##\$([_A-Za-z0-9]+) *= *\(([0-9]+)\.\.([0-9]+)\)(.*)')
    lines = map(string.rstrip,lines)
    j=0
    retval =  bruker_match_line(lines[j],number_re,string_re,array_re)
    j = j+1
    retval2 =  bruker_match_line(lines[j],number_re,string_re,array_re) #always grab the second line
    while j < len(lines):
        isdata = False
        if retval[0]==1 or retval[0]==2:
            name = retval[1]
            thislen = retval[2]
            data = retval[3]
            while (retval2[0] == 3) and (j<len(lines)): # eat up the following lines
                data += ' '+retval2[1]
                j = j+1
                retval2 =  bruker_match_line(lines[j],number_re,string_re,array_re)
            isdata = True
        elif retval[0]==0:
            name = retval[1]
            data = retval[2]
            isdata = True
        #else:
        #   print 'not a data line:',retval[1]
        if(isdata):
            if retval[0]==2: #if it's an array
                data = data.split(' ')
                if len(data)>0:
                    while '' in data:
                        data.remove('')
                    data = map(double,data)
                    if len(data)-1!= thislen[1]:
                        print 'error:',len(data)-1,'!=',thislen[1]
            vars.update({name:data})
        # at this point, the string or array data is loaded into data and we have something in retval2 which is definitely a new line
        retval = retval2
        j = j+1
        if j<len(lines):
            retval2 =  bruker_match_line(lines[j],number_re,string_re,array_re)
    fp.close()
    return vars
#}}}
#{{{ winepr_load_acqu
def winepr_load_acqu(file):
    fp = open(file+'.par','rU') # the U automatically converts dos format
    lines = fp.readlines()
    vars = {}
    line_re = re.compile(r'([_A-Za-z0-9]+) +(.*)')
    lines = map(string.rstrip,lines)
    #lines = [j.rstrip('\n') for j in lines] # because it's just \n, even on windows
    v = {'DRS':4096,'RES':1024,'HSW':50}
    for line in lines:
        m = line_re.match(line)
        if m is None:
            print 'Warning:',lsafen(repr(line)),'does not appear to be a valid WinEPR format line, and I suspect this is a problem with the terminators!'
        else:
            name = m.groups()[0]
            value = m.groups()[1]
            try:
                value = int(value)
            except:
                try:
                    value = double(value)
                except:
                    pass
            v[name]=value
    jss = long(v['JSS'])
    parameters = [ 'DUAL', '2D', 'FT', 'MAN0', 'MAN1', 'PROT', 'VEPR', 'POW', 'ABS', 'FTX', 'FTY', 'POW2', 'ABS2']
    parameters = map((lambda x: 's_'+x),parameters)
    masks = [ 0x00000001L, 0x00000002L, 0x00000004L, 0x00000008L, 0x00000010L, 0x00000020L, 0x00000040L, 0x00000080L, 0x00000100L, 0x00000200L, 0x00000400L, 0x00000800L, 0x00001000L]
    values = map((lambda x: x&jss),masks)
    values = map(bool,values)
    values = map(bool,values)
    v.update(dict(zip(parameters,values)))
    return v
#}}}
#}}}
#}}}
#}}}
#{{{ higher level functions
#{{{ standard EPR processing
def standard_epr(dir = None,
        files = None,# this is a single string for the file, a
        #list of files, each of which can optionally be a tuple
        #with a sample length (to calculate a "concentration")
        normalize_field = False,
        find_maxslope = False,
        offset_spectra = False,
        subtract_first = False,
        normalize_peak = False,
        background = None,
        grid = False,
        figure_list = None):
    if type(files) is str or type(files) is tuple:
        files = [files]
    if background is not None:
        if type(background) is list:
            if len(background) > 1:
                raise CustomError("multiple backgrounds not yet supported")
            else:
                background = background[0]
        files = [background] + files
        subtract_first = True
    dir = dirformat(dir)
    mu_B = 9.27400915e-24
    h = 6.62606957e-34
    if subtract_first:
        firstdata = 2*load_indiv_file(dir+files.pop(0))
    legendtext = list(files)
    for j in range(0,len(legendtext)):
        if type(legendtext[j]) is tuple:
            legendtext[j] = legendtext[j][0]
    print lsafen("legendtext is",legendtext,"to start with")
    for index,file in enumerate(files):
        if type(file) is tuple:
            sample_len = file[1]
            file = file[0]
        else:
            sample_len = None
        try:
            data = load_indiv_file(dir+file)
        except:
            raise CustomError('Error loading file: type of dir',type(dir),'type of file',type(file),'file=',file)
        if subtract_first:
            data -= firstdata
        field = r'$B_0$'
        neworder = list(data.dimlabels)
        data.reorder([neworder.pop(neworder.index(field))]+neworder) # in case this is a saturation experiment
        data -= data.copy().run_nopop(mean,field)
        figure_list.next('epr')
        v = winepr_load_acqu(dir+file)
        if index == 0:
             fieldbar = data[field,lambda x: logical_and(x>x.mean(),x<x.mean()+10.)]
             fieldbar.data[:] = 0.5
             fieldbar.data[0] = 0.6
             fieldbar.data[-1] = 0.6
             fxaxis = fieldbar.getaxis(field)
        xaxis = data.getaxis(field)
        centerfield = None
        if normalize_field:
            xaxis /= v['MF']
            if index == 0:
                fxaxis /= v['MF']
            newname = r'$B_0/\nu_e$'
        else:
            deriv = data.copy()
            deriv.run_nopop(diff,field)
            deriv.data[abs(data.data) > abs(data.data).max()/10] = 0 # so it doesn't give a fast slope, non-zero area
            deriv = abs(deriv)
            deriv.argmax(field,raw_index = True)
            centerfield = mean(xaxis[int32(deriv.data)])
            if find_maxslope:
                xaxis -= centerfield
                if index == 0:
                    fxaxis -= centerfield
                newname = r'$\Delta B_0$'
            else:
                newname = field
        data.rename(field,newname)
        if index == 0:
            fieldbar.rename(field,newname)
        mask = data.getaxis(newname)
        mask = mask > mask[int32(len(mask)-len(mask)/20)]
        snr = abs(data.data).max()/std(data.data[mask])
        integral = data.copy()
        integral.data -= integral.data.mean() # baseline correct it
        integral.integrate(newname)
        figure_list.next('epr_int')
        pc = plot_color_counter()
        plot_color_counter(pc)
        plot(integral,alpha=0.5,linewidth=0.3)
        integral.integrate(newname)
        figure_list.next('epr')
        if normalize_peak:
           normalization = abs(data).run_nopop(max,newname)
           data /= normalization
        ax = gca()
        if offset_spectra:
            myoffset = array(ax.get_ylim()).min()
        else:
            myoffset = 0.
        plot_color_counter(pc)
        plot(data+myoffset,alpha=0.5,linewidth=0.3,label=file)
        minval = abs(data.getaxis(newname)-centerfield).argmin()
        centerpoint = data[newname,minval]
        #plot_color_counter(pc)
        #plot(centerfield,centerpoint.data+myoffset,'o',markersize = 5,alpha=0.3)
        axis('tight')
        if index == 0:
            fieldbar *= array(ax.get_ylim()).max()
        legendtext[index] += '\n'
        if centerfield != None:
            legendtext[index] += ', %0.03f $G$'%centerfield
            gfactor = h*v['MF']*1e9/mu_B/(centerfield*1e-4)
            legendtext[index] += ', g=%0.03f'%gfactor
        if sample_len is None:
            legendtext[index] += r', SNR %0.2g $\int\int$ %0.3g'%(snr,integral[newname,-1].data[-1])
        else:
            legendtext[index] += r', SNR %0.2g $\int\int$ %0.3g / cm'%(snr,integral[newname,-1].data[-1]/double(sample_len))
    #xtl = ax.get_xticklabels()
    #at.xaxis.tick_top()
    #map( (lambda x: x.set_visible(False)), xtl)
    plot_color_counter(pc)
    plot(data.getaxis(newname)[mask],zeros(shape(data.getaxis(newname)[mask])),'k',alpha=0.2,linewidth=10)
    figure_list.next('epr')
    plot(fieldbar,'k',linewidth = 2.0)
    if grid:
        gridandtick(gca())
    #autolegend()
    axis('tight')
    figure_list.next('epr_int')
    autolegend(legendtext)
    axis('tight')
    return figure_list
#}}}
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
            plot_newdata.integrate('t2') # this is apparently a function to do integral with all the correct bells and whistles
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
        v = bruker_load_acqu(dirformat(filename))
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
    v = bruker_load_acqu(r'%s%d/'%(path,j))
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
            retval += ['%d: '%j+bruker_load_title(r'%s%d'%(path,j))+' $t_{dwov}$ %0.1f RG %d, DE %0.2f, mean %0.1f'%(dwov*1e6,rg,de,avg)]
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
            retval += ['%d: '%j+bruker_load_title(r'%s%d'%(path,j))+' $t_{dwov}$ %0.1f RG %d, DE %0.2f, mean %0.1f'%(dwov*1e6,rg,de,avg)]
            axis('tight')
        if retplot:
            return pval,retval
        else:
            return retval
    else:
        return []
#}}}
#}}}
