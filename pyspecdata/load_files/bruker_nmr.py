from ..core import *
from ..datadir import dirformat
from .open_subpath import open_subpath
import os.path
import re
import struct
import numpy as np
from numpy import nan

bruker_data = nddata # should work by inheritance but doesn't

def det_phcorr(v):
    if v['DIGMOD']==1:
        logger.debug('DIGMOD is 1')
        # table from Matlab program from C. Hilty
        gdparray=np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[179,201,533,709,1097,1449,2225,2929,4481,5889,8993,11809,18017,23649,36065,47329,72161,94689,144353,189409,288737],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[184,219,384,602,852,1668,2292,3368,4616,6768,9264,13568,18560,27392,36992,55040,73856,110336,147584,220928,295040]])
        decimarray=np.array([2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024]) # the -1 is because this is an index, and copied from matlab code!!!
        dspfvs = v['DSPFVS']
        decim = v['DECIM']
        if 'GRPDLY' not in list(v.keys()):
            grpdly = -1
        grpdly = v['GRPDLY'] # later versions of topspin
        if grpdly == -1:
            try:
                retval = gdparray[dspfvs,where(decimarray==decim)[0]]//2/decim
            except:
                if len(where(decimarray==decim)[0]) == 0:
                    raise CustomError("Not able to find decim",decim,"in decimarray")
                raise CustomError('Problem returning',dspfvs,where(decimarray==decim)[0],'elements of gdparray from gdparray of size',shape(gdparray),'because decimarray is of size',shape(decimarray))
                retval = 0
            return retval
        else:
            return grpdly
    else:
        logger.debug('DIGMOD is %d'%v['DIGMOD'])
        if v['DIGMOD']==3 and 'GRPDLY' in list(v.keys()):
            logger.debug('GRPDLY is %f'%v['GRPDLY'])
            return v['GRPDLY']
        return np.array([0])
def det_rg(a):
    '''determine the actual voltage correction from the value of rg for a bruker NMR file'''
    return a
def match_line(line,number_re,string_re,array_re):
    m = number_re.match(line)
    if m:
        retval = (0,m.groups()[0],np.double(m.groups()[1]))
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
                name = m.groups()[0]
                thislen = (np.double(m.groups()[1]),np.double(m.groups()[2]))
                thisdata = m.groups()[3]
                retval = (2,name,thislen,thisdata)
            else:
                retval = (3,line)
    return retval
def series(file_reference, *subpath, **kwargs):
    """For opening Bruker ser files.  Note that the expno is included as part of the subpath.
    
    Parameters
    ----------
    filename:
        see :func:`open_subpath`
    subpath:
        the path within the directory or zip file that's one level up
        from the ser file storing the raw data
        (this is typically a numbered directory).
    """
    dimname = process_kwargs([('dimname',''),
        ],kwargs)
    #{{{ Bruker 2D
    v = load_acqu(file_reference, *subpath)
    v2 = load_acqu(file_reference, *subpath, whichdim='2')
    td2 = int(v['TD'])
    rg = det_rg(float(v['RG']))
    td1 = int(v2['TD'])
    td2_zf = int(np.ceil(td2/256.)*256) # round up to 256 points, which is how it's stored
    fp = open_subpath(file_reference,*(subpath+('ser',)), mode='rb')
    data = fp.read()
    fp.close()
    if int(v['BYTORDA']) == 1:
        data = np.fromstring(data, dtype=np.dtype('>i4'), count=(len(data)//4))
    else:
        data = np.fromstring(data, dtype=np.dtype('<i4'), count=(len(data)//4))
    data = np.complex128(data)
    data = data[0::2]+1j*data[1::2]
    data /= rg
    mydimsizes = [td1,td2_zf//2]
    mydimnames = [dimname]+['t2']
    try:
        data = bruker_data(data,mydimsizes,mydimnames)
    except:
        size_it_should_be = np.array(mydimsizes).prod()
        if size_it_should_be > len(data):
            zero_filled_data = np.zeros(size_it_should_be)
            zero_filled_data[0:len(data)] = data
            data = bruker_data(zero_filled_data,mydimsizes,mydimnames)
        else:
            new_guess = len(data)/(td2_zf//2)
            print(lsafen("WARNING!, chopping the length of the data to fit the specified td1 of ",td1,"points!\n(specified ",list(zip(mydimnames,mydimsizes)),' td2_zf=%d)'%td2_zf))
            data = data[0:size_it_should_be]
            data = bruker_data(data,mydimsizes,mydimnames)
    logger.debug(strm('data straight from nddata =',data))
    data = data['t2',0:td2//2] # now, chop out their zero filling
    t2axis = 1./v['SW_h']*r_[1:td2//2+1]
    t1axis = r_[0:td1]
    mylabels = [t1axis]+[t2axis]
    data.labels(mydimnames,mylabels)
    shiftpoints = int(det_phcorr(v)) # use the canned routine to calculate the first order phase shift
    data.setaxis('t2',lambda x: x-shiftpoints/v['SW_h'])
    data.set_units('t2','s')
    data.set_units('digital')
    logger.debug("I'm about to try to load the title from: "+strm(file_reference,*subpath))
    data.set_prop('title',
            load_title(file_reference, *subpath))
    SFO1 = v['SFO1']
    BF1 = v['BF1']
    O1 = v['O1']
    v['TD2'] = int(v['TD'])
    del v['TD']
    v.update(v2)
    v['TD1'] = int(v['TD'])
    del v['TD']
    if v['SFO1'] != SFO1:
        # for, e.g. 2H experiments, a bad SFO1 (1H) is stored in acqu2, which we don't want
        print("warning: ignoring second dimension SFO1, since it's probably wrong")
        v['SFO1'] = SFO1
        v['O1'] = O1
        v['BF1'] = BF1
    with open_subpath(file_reference, *(subpath+('pulseprogram',)),mode='r') as fp:
        ppg = fp.read()
        if type(ppg) is bytes:
            ppg = ppg
        data.set_prop('pulprog',ppg)
    data.set_prop('acq',
            v)
    if isinstance(file_reference,str):
        data.set_prop('file_reference',
                file_reference)
    else:
        # if it's a zip file
        data.set_prop('file_reference',
                file_reference[1])
    if open_subpath(file_reference, *(subpath +
        ('pdata','1','procs')), test_only=True):
        data.set_prop('proc',
                load_jcamp(file_reference, *(subpath +
                    ('pdata','1','procs'))))
    if open_subpath(file_reference, *(subpath+('vdlist',)), test_only=True):
        data.set_prop('vd',
                load_vdlist(file_reference, *subpath))
    else:
        logger.info(strm("vdlist doesn't exist",file_reference,'vdlist'))
    if open_subpath(file_reference, *(subpath+('gradient_calib',)), test_only=True):
        logger.debug(strm("I found a gradient_calib file"))
        fp = open_subpath(file_reference, *(subpath+('gradient_calib',)),mode='r')
        data.set_prop('gradient_calib',
                fp.read())
        fp.close()
    if open_subpath(file_reference, *(subpath+('difflist',)), test_only=True):
        data.set_prop('diff',
                load_vdlist(file_reference, *subpath,
                    **dict(name='difflist')))
    else:
        logger.debug(strm("difflist doesn't exist",file_reference,'difflist'))
    logger.debug(strm('data from bruker file =',data))
    #}}}
    return data
def load_1D(file_reference, *subpath, **kwargs):
    """Load 1D bruker data into a file.  Load acquisition parameters into
    property 'acq' and processing parameters *from procno 1 only* into
    'proc'
    
    Note that is uses the 'procs' file, which appears to contain the correct data
    """
    dimname = process_kwargs([('dimname','')], kwargs)
    v = load_acqu(file_reference, *subpath)
    td2 = int(v['TD'])
    td1 = 1
    td2_zf = int(np.ceil(td2/256.)*256) # round up to 256 points, which is how it's stored
    fp = open_subpath(file_reference, *(subpath+('fid',)),mode='rb')
    data = fp.read()
    if int(v['BYTORDA']) == 1:
        data = np.fromstring(data, dtype=np.dtype('>i4'), count=(len(data)//4))
    else:
        data = np.fromstring(data, dtype=np.dtype('<i4'), count=(len(data)//4))
    data = np.complex128(data)
    data = data[0::2]+1j*data[1::2]
    rg = det_rg(v['RG'])
    data /= rg
    data = bruker_data(data,[td1,td2_zf//2],[dimname,'t2'])
    data = data['t2',0:td2//2] # now, chop out their zero filling
    t2axis = 1./v['SW_h']*r_[1:td2//2+1]
    t1axis = r_[1]
    data.labels([dimname,'t2'],[t1axis,t2axis])
    shiftpoints = int(det_phcorr(v)) # use the canned routine to calculate the second order phase shift
    #print 'shiftpoints = ',shiftpoints
    data.setaxis('t2',lambda x: x-shiftpoints/v['SW_h'])
    logger.debug('yes, I called with %d shiftpoints'%shiftpoints)
    # finally, I will probably need to add in the first order phase shift for the decimation --> just translate this
    data.set_prop('title',
            load_title(file_reference,*subpath))
    data.set_prop('acq',
            v)
    with open_subpath(file_reference, *(subpath+('pulseprogram',)),mode='r') as fp:
        ppg = fp.read()
        data.set_prop('pulprog',ppg)
    if type(file_reference) is tuple:
        data.set_prop('filename',
                file_reference[1])
    else:
        data.set_prop('filename',
                file_reference)
    if open_subpath(file_reference,
            *(subpath+('pdata','1','procs')),
            test_only=True):
        data.set_prop('proc',
                load_jcamp(file_reference,
                    *(subpath+('pdata','1','procs'))))
    return data
def load_vdlist(file_reference, *subpath, **kwargs):
    name = process_kwargs([('name','vdlist')], kwargs)
    subpath += (name,)
    print("subpath is",subpath)
    fp = open_subpath(file_reference,*subpath)
    lines = fp.readlines()
    if isinstance(lines[0],bytes):
        lines = map(lambda x: x, lines)
    lines = list(map(lambda x: x.rstrip(),lines))
    lines = list(map((lambda x: x.replace('m','e-3')),lines))
    lines = list(map((lambda x: x.replace('s','')),lines))
    lines = list(map((lambda x: x.replace('u','e-6')),lines))
    lines = list(map(np.double,lines))
    fp.close()
    return np.array(lines)
def load_acqu(file_reference,*subpath,**kwargs):
    """based on file_reference, determine the jcamp file that stores the acquisition info, and load it

    Parameters
    ----------
    file_reference:
        the file reference -- see open_subpath
    subpath:
        the subpath -- see open_subpath
    whichdim: optional string
        Default null string -- for multi-dimensional data, there is a separate acqu file for each dimension
    return_s: bool
        Default True -- whether to return the parameters of the saved data, or those manipulated since then.
    """
    whichdim, return_s = process_kwargs([('whichdim',''),
        ('return_s',True)],kwargs)
    if return_s:
        subpath += ('acqu'+whichdim+'s',) # this is what I am initially doing, and what works with the matched filtering, etc, as is, but it's actually wrong
    else:
        subpath += ('acqu'+whichdim,) # this is actually right, but doesn't work with the matched filtering, etc.
    return load_jcamp(file_reference,*subpath)
def load_jcamp(file_reference,*subpath):
    "return a dictionary with information for a jcamp file"
    def convert_to_num(val):
        if val == '<>':
            return nan
        elif val[0] == '<' and val[-1] == '>':
            return val[1:-1]
        elif '.' in val:
            return np.double(val)
        elif 'e-' in val.lower():
            return np.double(val)
        else:
            return int(val)
    fp = open_subpath(file_reference,*subpath)
    lines = fp.readlines()
    if isinstance(lines[0],bytes):
        lines = map(lambda x: x, lines)
    vars = {}
    number_re = re.compile(r'##\$([_A-Za-z0-9]+) *= *([0-9\-\.]+)')
    string_re = re.compile(r'##\$([_A-Za-z0-9]+) *= *<(.*)')
    array_re = re.compile(r'##\$([_A-Za-z0-9]+) *= *\(([0-9]+)\.\.([0-9]+)\)(.*)')
    lines = list(map(lambda x: x.rstrip(),lines))
    j=0
    retval =  match_line(lines[j],number_re,string_re,array_re)
    j = j+1
    retval2 =  match_line(lines[j],number_re,string_re,array_re) #always grab the second line
    while j < len(lines):
        isdata = False
        if retval[0]==1 or retval[0]==2:
            name = retval[1]
            thislen = retval[2]
            data = retval[3]
            while (retval2[0] == 3) and (j<len(lines)): # eat up the following lines
                data += ' '+retval2[1]
                j = j+1
                retval2 =  match_line(lines[j],number_re,string_re,array_re)
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
                    data = list(map(convert_to_num,data))
                    if len(data)-1!= thislen[1]:
                        print('error:',len(data)-1,'!=',thislen[1])
            vars.update({name:data})
        # at this point, the string or array data is loaded into data and we have something in retval2 which is definitely a new line
        retval = retval2
        j = j+1
        if j<len(lines):
            try:
                retval2 =  match_line(lines[j],number_re,string_re,array_re)
            except:
                raise RuntimeError(strm('matching line "',lines[j],'"',"in file",file_reference))
    fp.close()
    return vars
def load_title(file_reference,*subpath):
    logger.debug("I'm attempting to load the subpath with the following arguments: "+strm(
        file_reference,subpath,'pdata','1','title'))
    if not open_subpath(file_reference,*(subpath + ('pdata','1','title')), test_only=True):
        return None
    else:
        fp = open_subpath(file_reference,*(subpath + ('pdata','1','title')))
        lines = fp.readlines()
        logger.debug("I get %d lines"%len(lines))
        if len(lines) == 0:
            lines = []
            logger.warning("You do not have a title set -- this is highly unusual!!")
        else:
            lines = [j for j in lines if j not in ['\r\n','\n']]
        fp.close()
        return ''.join(lines)
