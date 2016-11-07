from ..core import *
from ..datadir import dirformat
import os.path
import re
import string
import struct
def det_phcorr(v):
    if v['DIGMOD']==1:
        gdparray=array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[179,201,533,709,1097,1449,2225,2929,4481,5889,8993,11809,18017,23649,36065,47329,72161,94689,144353,189409,288737],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[184,219,384,602,852,1668,2292,3368,4616,6768,9264,13568,18560,27392,36992,55040,73856,110336,147584,220928,295040]])
        decimarray=array([2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024]) # the -1 is because this is an index, and copied from matlab code!!!
        dspfvs = v['DSPFVS']
        decim = v['DECIM']
        if 'GRPDLY' not in v.keys():
            grpdly = -1
        grpdly = v['GRPDLY'] # later versions of topspin
        if grpdly == -1:
            try:
                retval = gdparray[dspfvs,where(decimarray==decim)[0]]/2/decim
            except:
                if len(where(decimarray==decim)[0]) == 0:
                    raise CustomError("Not able to find decim",decim,"in decimarray")
                raise CustomError('Problem returning',dspfvs,where(decimarray==decim)[0],'elements of gdparray from gdparray of size',shape(gdparray),'because decimarray is of size',shape(decimarray))
                retval = 0
            return retval
        else:
            return grpdly
    else:
        return array([0])
def det_rg(a):
    '''determine the actual voltage correction from the value of rg for a bruker NMR file'''
    return a
def match_line(line,number_re,string_re,array_re):
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
                name = m.groups()[0]
                thislen = (double(m.groups()[1]),double(m.groups()[2]))
                thisdata = m.groups()[3]
                retval = (2,name,thislen,thisdata)
            else:
                retval = (3,line)
    return retval
def dirformat(filename):
    "adds a path separator at the end of the file as needed"
    if os.path.isdir(filename):
        thislist = list(os.path.split(filename))
        thislist = filter(lambda x: len(x)>0,thislist)
        thislist = thislist + ['']
    else:
        raise ValueError(strm(filename,"doesn't appear to be a directory!!"))
    return os.path.sep.join(thislist)
def series(filename, dimname=''):
    "For opening Bruker ser files"
    #{{{ Bruker 2D
    filename = dirformat(filename)
    v = load_acqu(filename)
    v2 = load_acqu(filename,whichdim='2')
    td2 = int(v['TD'])
    rg = det_rg(float(v['RG']))
    td1 = int(v2['TD'])
    td2_zf = int(ceil(td2/256.)*256) # round up to 256 points, which is how it's stored
    fp = open(filename+'ser','rb')
    data = fp.read()
    if int(v['BYTORDA']) == 1:
        data = fromstring(data, dtype=dtype('>i4'), count=(len(data)/4))
    else:
        data = fromstring(data, dtype=dtype('<i4'), count=(len(data)/4))
    data = complex128(data)
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
    shiftpoints = int(det_phcorr(v)) # use the canned routine to calculate the first order phase shift
    data.circshift('t2',shiftpoints)
    data.set_units('t2','s')
    data.set_units('digital')
    data.set_prop('title',
            load_title(filename))
    #print 'DEBUG 2: data from bruker file =',data
    #}}}
    return data
def load_1D(filename, dimname=''):
    filename = dirformat(filename)
    v = load_acqu(filename)
    td2 = int(v['TD'])
    td1 = 1
    td2_zf = int(ceil(td2/256.)*256) # round up to 256 points, which is how it's stored
    fp = open(filename+'fid','rb')
    data = fp.read()
    if int(v['BYTORDA']) == 1:
        data = fromstring(data, dtype=dtype('>i4'), count=(len(data)/4))
    else:
        data = fromstring(data, dtype=dtype('<i4'), count=(len(data)/4))
    data = complex128(data)
    data = data[0::2]+1j*data[1::2]
    rg = det_rg(v['RG'])
    data /= rg
    data = nddata(data,[td1,td2_zf/2],[dimname,'t2'])
    data = data['t2',0:td2/2] # now, chop out their zero filling
    t2axis = 1./v['SW_h']*r_[1:td2/2+1]
    t1axis = r_[1]
    data.labels([dimname,'t2'],[t1axis,t2axis])
    shiftpoints = int(det_phcorr(v)) # use the canned routine to calculate the second order phase shift
    #print 'shiftpoints = ',shiftpoints
    data.circshift('t2',shiftpoints)
    # finally, I will probably need to add in the first order phase shift for the decimation --> just translate this
    data.set_prop('title',
            load_title(filename))
    return data
def load_vdlist(file):
    fp = open(file+'vdlist')
    lines = fp.readlines()
    lines = map(string.rstrip,lines)
    lines = map((lambda x: x.replace('m','e-3')),lines)
    lines = map((lambda x: x.replace('s','')),lines)
    lines = map((lambda x: x.replace('u','e-6')),lines)
    lines = map(double,lines)
    return array(lines)
def load_acqu(file,whichdim='',return_s = True):
    file = dirformat(file)
    def convert_to_num(val):
        if val == '<>':
            return NaN
        elif val[0] == '<' and val[-1] == '>':
            return val[1:-1]
        else:
            return double(val)
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
                    data = map(convert_to_num,data)
                    if len(data)-1!= thislen[1]:
                        print 'error:',len(data)-1,'!=',thislen[1]
            vars.update({name:data})
        # at this point, the string or array data is loaded into data and we have something in retval2 which is definitely a new line
        retval = retval2
        j = j+1
        if j<len(lines):
            retval2 =  match_line(lines[j],number_re,string_re,array_re)
    fp.close()
    return vars
def load_title(file):
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
