# just make this a library of all NMR reading software
from matlablike import *
import re
import string
import struct
import os
import fornotebook
from scipy.io import loadmat

def OUTPUT_notebook():
    return True
#{{{ general, non file-format specific functions
def dbm_to_power(dbm):
    return 1e-3*10.0**((dbm+40.0)/10.0) #20 db for each atten
def power_to_dbm(power):
    return 10.0*log(power/1e-3)/log(10.0)-40 #20 db for each atten
#{{{ auto_steps
def auto_steps(filename,threshold = -35, upper_threshold = 0, t_minlength = 0.5*60,minstdev = 0.1,showplots = True, showdebug = False,t_start=0,t_stop=60*1000,tolerance = 2,t_maxlen = inf):
    r'Plot the raw power output in figure 1, and chop into different powers in plot 2'
    v = loadmat(filename)
    p_ini = v['powerlist']
    t_ini = v['timelist']
    #{{{ plot the raw power log, then just pull out the times we're interested in
    if showplots:
        figure(1)
        plot(t_ini/60,p_ini)
        title(filename)
    mask = t_ini < t_stop
    minsteps = sum(t_ini<t_minlength)
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
    while nextpos < len(t)-1:
        figure(2)
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
        maxdiff_index = nextpos
        maxdiff = 0
        while (nextpos < len(t)-1) and (abs(p[nextpos]-curmean)<tolerance*stdev or ~isfinite(p[nextpos])):
            if abs(p[nextpos]-curmean) > maxdiff:
                maxdiff_index = nextpos
                maxdiff = abs(p[nextpos]-curmean)
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
            if t[nextpos]-t[blockstart] > t_maxlen:
                nextpos = maxdiff_index+1
                break
        track_curmean[nextpos] = curmean
        #uptohere = flattened[0:nextpos]
        #uptohere[uptohere==0] = cursum/curnum
        #flattened[nextpos-1] = cursum/curnum
        subset = flattened[0:nextpos] # subset of all the flattened stuff up to here
        subset[isnan(subset)] = curmean
        spikes[nextpos] = maxp # show a spike to see clearly where the data is divided
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
    return array(powerlist)
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
def load_acqu(filename):
    filename = dirformat(filename)
    if det_type(filename)[0] == 'bruker':
        return bruker_load_acqu(filename)
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
def load_file(filenames,dimname='',calibration=1.0,add_sizes = [], add_dims = []):
    if type(filenames) is not list:
        filenames = [filenames]
    #{{{load all the data into a list
    data = [load_indiv_file(filenames[0],dimname=dimname,add_sizes = add_sizes,add_dims = add_dims)]
    for filename in filenames[1:]:
        data += [load_indiv_file(filename,dimname=dimname,add_sizes = add_sizes,add_dims = add_dims)]
    #}}}
    # for the following, I used to have a condition, but this is incompatible with the pop statement at the end
    newdata = concat(data,dimname) # allocate the size of the indirect array
    newdata_shape = ndshape(newdata)
    if all(map((lambda x:det_type(x)[0]=='prospa'),filenames)):
        if hasattr(data[0],'want_to_prospa_decim_correct'):
            if data[0].want_to_prospa_decim_correct is True:
                newdata = prospa_decim_correct(newdata)
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
        if len(add_sizes)>0:
            mydimsizes = [td1/prod(add_sizes)]+add_sizes+[td2_zf/2]
        else:
            mydimsizes = [td1,td2_zf/2]
        mydimnames = [dimname]+add_dims+['t2']
        data = nddata(data,mydimsizes,mydimnames)
        data = data['t2',0:td2/2] # now, chop out their zero filling
        t2axis = 1./v['SW_h']*r_[1:td2/2+1]
        t1axis = r_[0:td1]
        add_labels = []
        for j in range(0,len(add_dims)):
            add_labels += [r_[0:add_sizes[j]]]
        mylabels = [t1axis]+add_labels+[t2axis]
        data.labels(mydimnames,mylabels)
        shiftpoints = int(bruker_det_phcorr(v)) # use the canned routine to calculate the first order phase shift
        data.circshift('t2',shiftpoints)
        if return_acq:
            return (data,v,v2)
        else:
            return data
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
        data = prospa_load_datafile(filename,dims=dimshere)
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
        return data
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
            if return_acq:
                return (data,v,v2)
            else:
                return data
        #}}}
        #{{{ prospa 1d
        elif filetype == 'prospa':
            v = prospa_load_acqu(filename)
            data = prospa_load_datafile(filename,dims=1)
            data = nddata(data,[v['nrPnts']],['t2'])
            taxis = linspace(0,1,v['nrPnts'])*v['acqTime']/1e3
            data.labels(['t2'],[taxis])
            return data
        #}}}
        else:
            raise CustomError("can't load this file type")
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
        decimarray=array([2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024])
        dspfvs = v['DSPFVS']
        decim = v['DECIM']
        return gdparray[dspfvs,where(decimarray==decim)[0]]/2/decim;
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
def bruker_load_acqu(file,whichdim=''):
    fp = open(file+'acqu'+whichdim+'s')
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
    v = {'DRS':4096,'RES':1024,'HSW':50}
    for line in lines:
        m = line_re.match(line)
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
#{{{ optimize the first-order phase the right way
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
def rg_check(file,expno,number_of_samples = 75,dynamic_range = 19,first_figure = 1,show_complex = True):
    r'''show the fid, plotted vs. the max possible value to check the max of the receiver gain in figure 1
    in figure 2, plot in the complex plain to check digitization limit'''
    data = load_emax(file,expno,printinfo = False) # load the data
    if len(expno)>0:
        params = load_acqu('%s/%d'%(dirformat(file),expno[0])) # load the parameters
    else:
        params = load_acqu(dirformat(file[0])) # load the parameters
    maxfloat = 3.4028234e38
    maxint = 2147483648
    max_for_this_dynamic_range= 2**(dynamic_range-1) # the -1 is the bit used for the sign
    data.reorder(['t2','power'])
    if type(file) is list:
        spec_type = det_type(file[0])[0]
    else:
        spec_type = det_type(file)[0]
    if spec_type == 'bruker':
        data.data *= params['RG'] # scale back up by RG, so we get raw numbers
    data.data /= max_for_this_dynamic_range # scale down by the max int number, so we test digitization
    figure(first_figure)
    plot(data['t2',0:number_of_samples],'k')
    plot(data['t2',0:number_of_samples]*(-1j),'b',alpha=0.3)
    title('receiver gain upper value check')
    axis('tight')
    if show_complex:
        figure(first_figure+1)
        OLDplot(real(data['t2',0:number_of_samples].data),imag(data['t2',0:number_of_samples].data),'.',markersize=0.6)
        #polar(angle(data['t2',0:number_of_samples].data),abs(data['t2',0:number_of_samples].data),'.',markersize=0.1)
        title('receiver gain minimum value check')
        axis('tight')
    return
#}}}
#{{{ integrate_emax --> old integration function
def integrate_emax(file,expno,
        integration_width=1e3,
        intpoints=None,
        showimage=False,
        usephase=True,
        filteredge = 10,
        center_peak=False,
        usebaseline=False,
        plotcheckbaseline=False,
        filter_direct = False,
        return_noise=False,
        show_integral = False,
        indiv_phase = False,
        abs_image = False,
        add_sizes = [],
        add_dims = [],
        add_dims_cycle = []):
    r''' old (obsolete) integration function'''
    #print 'diagnose: before emax'
    data = load_emax(file,expno,add_sizes = add_sizes,add_dims = add_dims) # load the data
    for j in range(0,len(add_dims)):
        if bool(add_dims[j].find('ph')==0):# starts with ph
            data.ft(add_dims[j])
            data = data[add_dims[j],add_dims_cycle[j]] # pull the correct phase cycle --> awesome!
    see_fid = False
    if see_fid:
        #{{{ for debug
        data.reorder(['t2','power'])
        plot(abs(data['t2',0:500]))
        return r_[0,1,2],r_[0,1,2]
        #}}}
    debug_ini_data_size = ndshape(data)
    data.ft('t2',shiftornot=True) # ft along t2
    #{{{ apply matched filter, which I've decided should be obsolete
    if filter_direct:
        data.ift('t2',shiftornot=True) # ft along t2
        filter = matched_filter(data,'t2')
        data *= filter
        data.ft('t2',shiftornot=True) # ft along t2
    #}}}
    data_shape = ndshape(data)
    #{{{ make the absolute value, so we can find the top of the peak
    data_abs = data.copy()
    #{{{ apply the matched filter to maximize our SNR while picking the peak
    data_abs.ift('t2',shiftornot=True) # ft along t2
    filter = matched_filter(data_abs,'t2')
    data_abs *= filter
    data_abs.ft('t2',shiftornot=True) # ft along t2
    data_abs = abs(data_abs)
    #}}}
    #}}}
    #{{{ find the top of the peak, or alternatively, put an array into data_abs which we can use for matched filtering
    if not(center_peak):
        data_abs.mean('power')
    data_topvals = data_abs.copy()
    data_topvals.run(amax,'t2') # amax returns the max along the appropriate axis
    data_top = data_abs.copy() # since data_abs is currently not used, but I want to use it to do matched filtered integration, really need to make a separate variable here
    data_top.run(argmax,'t2') # put the index of the top peak there
    top = int32(data_top.data)
    topvals = data_topvals.data
    if not(center_peak):
        data_abs.data = ones(data_shape['power'])*data_abs.data
    #}}}
    #{{{ if integration points unspec, pull out the integration width
    if intpoints==None:
        df = data.getaxis('t2').copy()
        df = df[1]-df[0]
        intpoints = floor(integration_width/(df))
    #}}}
    #{{{actually pull the integration points, + and - intpoints from the peak, putting them into newdata and the noise into newnoise
    data_shape['t2'] = intpoints*2+1
    newdata = []
    newnoise = []
    top[top<intpoints] = intpoints # prevent a bug where the integration range exceeds the spectrum
    #{{{baseline correction
    if usebaseline:
        for j in range(0,ndshape(data)['power']):
            if plotcheckbaseline:
                hold(True)
                baseline_data = baseline_spectrum_peakpick(data['power',j].copy(),showplots = True)
            else:
                baseline_data = baseline_spectrum_peakpick(data['power',j].copy())
                data['power',j].data[:] -= baseline_data.data.flatten()
                if any(isnan(data['power',j].data[:])):
                    print 'isnan!!'
                if any(isinf(data['power',j].data[:])):
                    print 'isinf!!'
        if plotcheckbaseline:
            error_plot('Wanted to check baseline')
    #}}}
    for j in range(0,data_shape['power']):
        newdata += [data['power',j,'t2',top[j]-intpoints:top[j]+intpoints+1]]
        if return_noise:
            #newnoise += [data['power',j,'t2',10:10+intpoints]]
            newnoise += [data['power',j,'t2',top[j]+intpoints+1:top[j]+2*intpoints+1]] # pull noise from an integration range to the right of the current one --> this probably needs work
    debug_top = top
    debug_slice = [top[j]-intpoints,top[j]+intpoints+1]
    debug_newdata_len = len(newdata)
    newdata = concat(newdata,'power')
    if return_noise:
        newnoise = concat(newnoise,'power')
    debug_newdata_shape = ndshape(newdata)
    debug_newdata_data_shape = newdata.data.shape
    #}}}
    #{{{ if we're using absolute value, calculate the noiselevel --> just ignore this for now, since I don't really forsee using this
    if not usephase:
        noiselevel = data['t2',-(2*intpoints+1)-filteredge:-filteredge]
        noiselevel = abs(noiselevel).mean('t2').mean_nopop('power')
        newdata = abs(newdata)
    #}}}
    #{{{ show what we're integrating
    if showimage:
        clf()
        figure(1)
        try:
            #image(newdata)
            center_for_image = int32(sum(top*topvals)/sum(topvals))
            if abs_image:
                image(abs(data['t2',center_for_image-3*intpoints:center_for_image+3*intpoints+1]))
                title('2D plot of spectra')
            else:
                image(data['t2',center_for_image-3*intpoints:center_for_image+3*intpoints+1])
                title('2D plot of spectra')
        except:
            print 'top = ',top
            print 'topvals = ',topvals
            if(any(isnan(newdata.data))):
                raise CustomError('image isnan')
            elif(any(isinf(newdata.data))):
                raise CustomError('image isinf')
            else:
                clf()
                plot(debug_data_afterbaseline.T)
                error_plot()
                raise CustomError('cant make image from type',type(newdata.data),newdata.data.dtype,'data_shape:',data_shape,'debug_newdata_len:',debug_newdata_len,'debug_newdata_shape:',debug_newdata_shape,'debug_newdata_data_shape:',debug_newdata_data_shape,'debug_top:',debug_top,'debug_slice:',debug_slice,'debug_ini_data_size:',debug_ini_data_size,'debug_data_afterbaseline(',shape(debug_data_afterbaseline),'):',debug_data_afterbaseline)
    newdatacopy = newdata.copy()
    #}}}
    newdata.mean('t2') # integrate
    #{{{ either autophase, or subtract out the abs baseline 
    if usephase:
        if not indiv_phase:
            #newdata.data *= exp(-1j*angle(newdata.data.sum()))
            phaseoptval = phaseopt(newdata.data)
            newdata.data *= phaseoptval
            newdatacopy.data *= phaseoptval
            if return_noise:
                newnoise.data *= phaseoptval
        else:
            for j in range(0,len(newdata.data)):
                phcorr =  newdata['power',j]
                phcorr /= abs(phcorr)
                try:
                    newdatacopy['power',j] *= phcorr
                    newdata['power',j] *= phcorr
                except:
                    print 'shape of newdatacopy',ndshape(newdatacopy)
                    print 'shape of newdata',ndshape(newdata)
    else:
        newdata -= noiselevel
        newdatacopy.data -= noiselevel
    if showimage:
        figure(2)
        clf()
        plot((newdatacopy).reorder(['t2','power']))
        title('Peaks, zoomed in to integration region')
        if show_integral:
            #{{{this does work to plot the integral
            newdatacopy.integrate('t2') #newdatacopy.run_nopop(cumsum,'t2')
            gridandtick(gca())
            ax = gca()
            myxlim = ax.get_xlim()
            twinx()
            plot(newdatacopy,alpha=0.5)
            axes(ax)
            ax.set_xlim(myxlim)
            #}}}
        figure(1)
    #}}}
    if return_noise:
        newnoise.data = abs(newnoise.data)
        newnoise.data **= 2
        newnoise.mean('t2') # which should give us the sigma --> note that if there is no baseline correction, this will not be a sigma
        return real(newdata.data.flatten()),sqrt(newnoise.data.flatten()*sqrt(2.*intpoints+1.)) # we have to account for the signal averaging from the integral
    else:
        return real(newdata.data.flatten())
#}}}
#{{{ integrate --> new integration function
def integrate(file,expno,
        integration_width=1e3,
        intpoints=None,
        show_image=True,
        filter_edge = 10,
        center_peak=True,
        use_baseline=True,
        plot_check_baseline=False,
        filter_direct = False,
        return_noise=False,
        show_integral = False,
        indiv_phase = False,
        peak_within = 20e3,
        abs_image = False,
        max_drift = 1e3,
        first_figure = 1):
    r'''new integration function, which replaces integrate_emax, and is used to integrate data, as for Emax and T1 curves'''
    if type(plot_check_baseline) is bool:
        if plot_check_baseline:
            plot_check_baseline = 0 
        else:
            plot_check_baseline = -1
    data = load_emax(file,expno) # load the data
    # see_fid obsolete by rg_check
    # also remove all debug statements
    #data = data['power',r_[0,1,-2,-1]];print 'WARNING, DEBUG --> pulling only 3 powers!'
    data.ft('t2',shiftornot=True) # ft along t2
    data_shape = ndshape(data) # this is used to shape the output
    #{{{ abs w/ max SNR, so we can pick the peak
    data_abs = data.copy()
    t2temp = data_abs.getaxis('t2')
    data_abs['t2',abs(t2temp)>peak_within].data *= 0
    #{{{ apply the matched filter to maximize our SNR while picking the peak
    data_abs.ift('t2',shiftornot=True) # ft along t2
    filter = matched_filter(data_abs,'t2')
    data_abs *= filter
    data_abs.ft('t2',shiftornot=True) # ft along t2
    #}}}
    data_abs = abs(data_abs)
    #}}}
    #{{{ generate topavg --> the index at the top of the average
    data_mean = data_abs.copy()
    data_mean.mean('power') # note that we are taking the mean over the abs here, which would not be great for the noise, but the center should still be in the right place
    data_mean.run(argmax,'t2') # put the index of the top peak there
    topavg = int32(data_mean.data)
    #}}}
    #{{{ generate center --> an array with the value at the center of each scan
    data_center = data_abs.copy() # since data_abs is currently not used, but I want to use it to do matched filtered integration, really need to make a separate variable here
    f = data_center.getaxis('t2')
    data_center['t2',abs(f-f[topavg])>max_drift] = 0# we need to keep the indeces in place, but don't want to pick anything too far out of the way
    test_drift_limit = False
    if test_drift_limit:
        plot(data_center.reorder(['t2','power']))
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
    #{{{ show what we're integrating
    if show_image:
        clf()
        figure(first_figure)
        try:
            #image(newdata)
            center_for_image = topavg
            if abs_image:
                image(abs(data['t2',center_for_image-3*intpoints:center_for_image+3*intpoints+1]))
                title('2D plot of spectra')
            else:
                image(data['t2',center_for_image-3*intpoints:center_for_image+3*intpoints+1])
                title('2D plot of spectra')
        except:
            print 'center_for_image = ',center_for_image
            print 'topavg = ',topavg
            print 'intpoints = ',intpoints
            if(any(isnan(newdata.data))):
                raise CustomError('image isnan')
            elif(any(isinf(newdata.data))):
                raise CustomError('image isinf')
            else:
                clf()
                plot(debug_data_afterbaseline.T)
                error_plot()
                raise CustomError('cant make image from type',type(newdata.data),newdata.data.dtype,'data_shape:',data_shape,'debug_newdata_len:',debug_newdata_len,'debug_newdata_shape:',debug_newdata_shape,'debug_newdata_data_shape:',debug_newdata_data_shape,'debug_center:',debug_center,'debug_slice:',debug_slice,'debug_ini_data_size:',debug_ini_data_size,'debug_data_afterbaseline(',shape(debug_data_afterbaseline),'):',debug_data_afterbaseline)
    #}}}
    #{{{baseline correction
    if use_baseline:
        #data.data += 0.1 #for debug baseline
        for j in range(0,ndshape(data)['power']):
            if plot_check_baseline == j: # if I've passed this index, show this baseline
                baseline_data = baseline_spectrum(data['power',j].copy(),center[j],intpoints,showplots = True) # call baseline_spectrum with center and intpoints, which should already be defined
                error_plot('Wanted to check baseline on scan ',j)
            else:
                baseline_data = baseline_spectrum(data['power',j].copy(),center[j],intpoints) # call baseline_spectrum with center and intpoints, which should already be defined
                data['power',j].data[:] -= baseline_data.data.flatten()
                if any(isnan(data['power',j].data[:])):
                    print 'isnan!!'
                if any(isinf(data['power',j].data[:])):
                    print 'isinf!!'
    #}}}
    #{{{ actually pull the integral points and the standard deviation of the noise
    #print 'DEBUG: center points are',center
    plot_noise = [] # array where we can put the noise for plotting
    for j in range(0,data_shape['power']):
        newdata += [data['power',j,'t2',center[j]-intpoints:center[j]+intpoints+1]]
        if return_noise:
            #newnoise += [data['power',j,'t2',10:10+intpoints]]
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
            temp = data['power',j,'t2',list_of_noise_indeces]
            if show_image:
                plot_noise += [temp.copy()]
            temp.data = abs(temp.data-temp.data.mean()) # need to explicitly do the abs, since the data is complex
            temp.data **= 2
            temp.mean('t2') # I DO NOT take the sqrt here, because it's taken at the very end
            #}}}
            newnoise += [temp.copy()] # find the standard deviation of the noise which we have pulled --> it should be independent of the number of points that we're using
    newdata = concat(newdata,'power')
    if return_noise:
        newnoise = concat(newnoise,'power')
        if show_image:
            plot_noise = concat(plot_noise,'power')
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
            phcorr =  newdata['power',j]
            phcorr /= abs(phcorr)
            try:
                if show_image:
                    plot_newdata['power',j] *= phcorr
                newdata['power',j] *= phcorr
                if return_noise:
                    #newnoise['power',j] *= phcorr # again, don't rotate the real noise
                    if show_image:
                        plot_noise['power',j] *= phcorr
            except:
                print 'shape of newdatacopy',ndshape(newdatacopy)
                print 'shape of newdata',ndshape(newdata)
    #}}}
    #{{{ show what we're integrating
    if show_image:
        figure(first_figure+1)
        clf()
        plot_newdata.reorder(['t2','power'])
        plot(plot_newdata,alpha=0.5)
        title('Peaks, zoomed in to integration region')
        if return_noise:
            plot_noise.reorder(['t2','power'])
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
        figure(first_figure)
    #}}}
    #{{{ return integral and the noise of the 
    if return_noise:
        number_of_integral_points = 2.*intpoints+1. # the integral is a sum over this many points
        newdata.set_error(sqrt(newnoise.data.flatten()*number_of_integral_points))
    return newdata
    #}}}
#}}}
#}}}
#{{{ load the data from a emax series based on input array
def print_info(filename):
    filetype,twod = det_type(filename)
    if filetype == 'bruker':
        v = bruker_load_acqu(dirformat(filename))
        if OUTPUT_notebook():
            print r'\fn{%s}: sfo1:%0.5f aq:%0.3f swh:%0.3f ns:%d ds: %d rg:%0.1f d1:%0.1f p1:%0.2f pl1:%0.1f'%(filename,v['SFO1'],v['TD']/v['SW_h']/2.0,v['SW_h'],v['NS'],v['DS'],v['RG'],v['D'][1],v['P'][1],v['PL'][1])
            data = fornotebook.save_data()
            pptvalue = v['SFO1']/data['current_frequency']
            if abs(pptvalue-data['current_ppt'])>1e-4:
                print ('WARNING! ppt value is off!! (',pptvalue,' ppt)')
            else:
                print '%0.4f $ppt$'%pptvalue
            print '\n\n'
def load_emax(file,datafiles,printinfo=True,add_sizes = [],add_dims = []):
    datafiles = list(datafiles)
    if len(datafiles)==0:
        data = load_file(file,dimname='power',add_sizes = add_sizes,add_dims = add_dims) # got rid of a [0] here --> don't know why it was there --> load_file loads lists of files
        if printinfo:
            print_info(file[0])
    else:
        file = dirformat(file)
        datafiles = map(str,datafiles)
        spectype,dim = det_type(file+datafiles[0])
        if spectype == 'bruker':
            #{{{ concatenate all the files into one 
            data = []
            for j in range(0,len(datafiles)):
                data += [load_file(file+datafiles[j],dimname='power',add_sizes = add_sizes,add_dims = add_dims)]
            data = concat(data,'power')
            #}}}
            if printinfo:
                print_info(file+datafiles[0])
                print_info(file+datafiles[-1])
        elif spectype == 'prospa':
            #{{{ kea version
            datafiles = map((lambda x: file+str(x)), datafiles)
            vars = prospa_load_acqu(datafiles[0])
            for j in range(0,len(datafiles)):
                vars_check = prospa_load_acqu(datafiles[j])
                if (vars_check['acqTime']!=vars['acqTime']) or (vars_check['nrPnts']!=vars['nrPnts']):
                    raise CustomError('Problem! doesn\'t match first file:\n')
                    show_acqu(vars)
                    show_acqu(vars_check)
                    return
            npoints = vars['nrPnts']
            taxis = linspace(0,1,npoints)*vars['acqTime']/1e3
            #npoints = load_datafile(datafiles[0]).shape[0] #hack
            data_shape = ndshape([npoints,len(datafiles)],['t2','power'])
            data = data_shape.alloc()
            data.labels(['power','t2'],[range(0,len(datafiles)),taxis])
            for j in range(0,data_shape['power']):
                    data['power',j] = prospa_load_datafile(datafiles[j])
            data = prospa_decim_correct(data)
            #}}}
    return data
#}}}
#{{{ process_t1
def process_t1(file,expno,usebaseline = None,showimage = None,plotcheckbaseline = None,saturation = False,first_figure = 3,**kwargs):
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
    integral = integrate(file,expno,first_figure = first_figure,**kwargs)
    t1name = r'$t_1$'
    integral.rename('power',t1name)
    integral = t1curve(integral,fit_axis = t1name) # make this into an integral class, which fits along the dimension t1
    integral.labels([t1name],[wait_time]) # before, I had to sort them manually, but now, I don't
    #print 'DEBUG wait times:',integral.getaxis(t1name)
    integral.sort(t1name)
    #{{{ finally, show the fit  
    figure(first_figure+2)
    taxis = wait_time
    integral.data *= phaseopt(integral.data)
    plot(integral.runcopy(real),'ko')
    plot(integral.runcopy(imag),'yo')
    integral.makereal() # otherwise, it won't fit
    integral.fit()
    plot(integral.eval(taxis)) # evaluate the fit function on the axis taxis
    #{{{ for now, do not plot the modified versions
    #plot(taxis,t1_fitfunc(r_[p[0:2],p[2]*1.2],taxis),'y')
    #plot(taxis,t1_fitfunc(r_[p[0:2],p[2]*0.8],taxis),'y')
    #}}}
    ax = gca()
    text(0.5,0.5,integral.latex(),transform = ax.transAxes,size = 'x-large', horizontalalignment = 'center',color = 'g')
    title(titlestr)
    #}}}
    #{{{ and the straight line plot
    figure(first_figure+3)
    #print 'linear data:',integral.linear().data
    plot(integral.linear(),'o')
    plot(integral.linear(taxis))
    #print ndshape(integral.linear())
    #}}}
    return integral # there is never a return_fit, since the fit is stored in the class itsself
#}}}
#{{{ process cpmg 
def process_cpmg(file,dimname=''):
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
    data = load_file(r'%s%d'%(path,j),calibration=calibration)
    k_B = 1.3806504e-23
    data.ft('t2',shift = True)
    newt2 = r'F2 / $Hz$'
    data.rename('t2',newt2)
    v = bruker_load_acqu(r'%s%d/'%(path,j))
    dw = 1/v['SW_h']
    dwov = dw/v['DECIM']
    rg = v['RG']
    aq = v['TD']*dw
    if rg>rgmin:
        plotdata = abs(data)
        plotdata.data **= 2
        johnson_factor = 4.0*k_B*T
        plotdata.data /= (aq*johnson_factor)
        t = data.getaxis(newt2)
        mask = logical_and(t>mask_start,
            t<mask_stop)
        avg = plotdata.data[mask].mean() 
        retval = []
        if both or not smoothing:
            pval = plot(plotdata,'-',alpha=0.5,plottype = plottype)
            retval += ['%d: '%j+bruker_load_title(r'%s%d'%(path,j))+' $t_{dwov}$ %0.1f RG %d, mean %0.1f'%(dwov*1e6,rg,avg)]
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
            retval += ['%d: '%j+bruker_load_title(r'%s%d'%(path,j))+' $t_{dwov}$ %0.1f RG %d, mean %0.1f'%(dwov*1e6,rg,avg)]
            axis('tight')
        if retplot:
            return pval,retval
        else:
            return retval
    else:
        return []
#}}}
#{{{ different types of fit classes
class t1curve(fitdata):
    def guess(self):
        r'''provide the guess for our parameters, which is specific to the type of function'''
        x = self.getaxis(self.fit_axis)
        y = self.data
        testpoint = len(y)/3
        initial_slope = (y[testpoint]-y[0])/(x[testpoint]-x[0])
        A = y[-1]
        B = y[testpoint]-x[testpoint]*initial_slope
        C = (A-B)/initial_slope
        if (C < 0):
            print repr(CustomError(maprep('Negative T1!!! A-B=',A-B,'initial_slope=',initial_slope,x,y)))
        #else:
        #   print 'C is ',C
        return r_[A,B,C]
    def fitfunc_raw(self,p,x):
        '''just the actual fit function to return the array y as a function of p and x'''
        return p[0]+(p[1]-p[0])*exp(-x/p[2])
    def linfunc(self,x,y,xerr = None,yerr = None):
        '''just the actual fit function to return the pair of arrays x',y' that should be linear
        it accepts as inputs x and y, and it uses the output from the fit, where necessary
        also optionally propagates the error based on yerr and xerr, which can be passed in to it
        For the case of T1, we want to return ln(y-M(\infty)) = ln(M(0)-M(\infty)) - t/T_1
        '''
        #print 'DEBUG: y is',y
        #print 'DEBUG: M(\infty) is',self.output(r'M(\infty)')
        temp = self.output(r'M(\infty)')-y # the argument for log
        #print 'DEBUG: temp is',temp
        # note that there is some error associated with m(\infty) that I'm just taking for granted
        rety = log(temp)
        if yerr != None:
            reterr = yerr/abs(temp)
        mask = isfinite(rety)
        retx = x # for instance, in emax, this is not just x
        xname = self.fit_axis # same as the fit axis
        yname = r'$ln(M(\infty)-M(t))$'
        #{{{ this should be pretty standardized
        retval = nddata(rety,
                [size(rety),1],
                [xname,yname])
        retval.labels([self.fit_axis],
                [retx.copy()])
        if yerr != None:
            retval.set_error(reterr)
        #}}}
        return retval
    def linerror(self,x,y):
        '''propagate the error for linfunc
        '''
        rety = log(y-self.output(r'M(\infty)'))
        mask = isfinite(rety)
        x_axis_of_linear_plot = x # for instance, in emax, this is not just x
        retval = nddata(rety,
                [size(rety),1],
                [self.fit_axis,r'$ln(M(t)-M(\infty))$'])
        retval.labels([self.fit_axis],
                [x_axis_of_linear_plot.copy()])
        return retval
    def __init__(self,*args,**kwargs):
        '''here, we give the particular latex representation and list of symbols for this particular child class'''
        fitdata.__init__(self,*args,**kwargs)
        self.function_string = r'$M(t)=M(\infty)+(M(0)-M(\infty))e^{-t/T_1}$'
        self.symbol_list = [r'M(\infty)',r'M(0)',r'T_1'] # note that it must notbe possible to find part of one of the later strings by searching for one of the earlier strings
        return
class emax(fitdata):
    def guess(self):
        r'''provide the guess for our parameters, which is specific to the type of function'''
        newdata = self.copy()
        newdata.sort(self.fit_axis)
        power = newdata.getaxis(self.fit_axis)
        integral = newdata.data
        largest = len(power)
        initial_slope = (integral[largest/4]-integral[0])/(power[largest/4]-power[0])
        approx_emax = integral[-1]
        print 'DEBUG: guessed initial slope',initial_slope,'approx emax',approx_emax
        return [1,-initial_slope,-initial_slope/(1-approx_emax)]
    def fitfunc_raw(self,p,x):
        '''just the actual fit function to return the array y as a function of p and x'''
        return (p[0]-(p[1]*x/(1.+p[2]*x)))
    def linfunc(self,x,y,xerr = None,yerr = None):
        '''just the actual fit function to return the pair of arrays x',y' that should be linear
        it accepts as inputs x and y, and it uses the output from the fit, where necessary
        also optionally propagates the error based on yerr and xerr, which can be passed in to it
        For the case of E_max, we want 1/(1-E) = 
        '''
        # note that there is some error associated with m(\infty) that I'm just taking for granted
        #print "linfunc passed x=",x,"and y=",y
        rety = 1./(1.-y)
        if yerr != None:
            reterr = yerr/((1.-y)**2)
        mask = isfinite(rety)
        retx = 1./x # for instance, in emax, this is not just x
        xname = r'1 / '+self.fit_axis # same as the fit axis
        yname = r'$\frac{1}{1-E(p)}$'
        #{{{ this should be pretty standardized
        retval = nddata(rety,
                [size(rety),1],
                [xname,yname])
        retval.labels([xname],
                [retx.copy()])
        if yerr != None:
            retval.set_error(reterr)
        #}}}
        return retval
    def __init__(self,*args,**kwargs):
        '''here, we give the particular latex representation and list of symbols for this particular child class'''
        fitdata.__init__(self,*args,**kwargs)
        self.function_string = r'$E(t)=c_0-Ap/(1+Bp)$'
        self.symbol_list = [r'c_0',r'A',r'B'] # note that it must notbe possible to find part of one of the later strings by searching for one of the earlier strings
        return
#}}}
