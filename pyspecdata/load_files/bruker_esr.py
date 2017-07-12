from ..core import *
from ..general_functions import strm
from numpy import fromstring
import re, string
from StringIO import StringIO
b0_texstr = r'$B_0$'
def xepr(filename, dimname=''):
    """For opening Xepr files.
    
    Parameters
    ----------
    filename : str
        The filename that ends with either ``.dsc`` or ``.dta``.
    """
    # {{{ determine the pair of filenames that we need
    filename = filename[:-4]+filename[-4:].upper()# case insensitive extension
    if filename[-4:] == '.DTA':
        filename_spc,filename_par = filename,filename.replace('.DTA','.DSC')
    elif filename[-4:] == '.DSC':
        filename_spc,filename_par = filename.replace('.DSC','.DTA'),filename
    else:
        raise ValueError(strm("When guessing that the filename is a"
                " WinEPR file, the extension must be either .SPC or"
                " .PAR\n"
                "This one is called",repr(filename)))
    # {{{ check if the extension is upper or lowercase
    if not os.path.exists(filename_spc):
        filename_spc = filename_spc[:-4] + filename_spc[-4:].lower()
        filename_par = filename_par[:-4] + filename_par[-4:].lower()
    # }}}
    # }}}
    # {{{ load the data
    with open(filename_spc,'rb') as fp:
        data = fromstring(fp.read(),'>f8')
    # }}}
    # load the parameters
    v = xepr_load_acqu(filename_par)
    # {{{ flatten the dictionary (remove the uppermost/block
    #     level)
    new_v = {}
    for k_a,v_a in v.iteritems():
        new_v.update(v_a)
    v = new_v
    # }}}
    # {{{ use the parameters to determine the axes
    #     pop parameters that are just part of the axes
    x_points = v.pop('XPTS')
    print (x_points)
    x_axis = r_[0:x_points]
    # the following is NOT the same as *=, which preserves the
    # type (=bad!)!!!
    x_axis = x_axis*v.pop('XWID')/x_axis[-1] # wouldn't be a huge
    #         difference, but I'm not sure if this is correct
    #         (I think so)
    x_axis += v.pop('XMIN')
    harmonics = array([[False] * 5]*2) # outer dimension for the 90 degree phase
    for j,jval in enumerate(['1st','2nd','3rd','4th','5th']):
        for k,kval in enumerate(['','90']):
            thiskey = 'Enable'+jval+'Harm'+kval
            if thiskey in v.keys() and v[thiskey]:
                harmonics[k,j] = True
    n_harmonics = sum(harmonics)
    y_points = len(data)/x_points/n_harmonics
    if 'YPTS' in v.keys():
        raise ValueError(strm("I looks like this is a 2D file YPTS=",v['YPTS'],", which I"
            " just haven't bothered to program yet"))
    dimname_list = [b0_texstr]
    dimsize_list = [x_points]
    if n_harmonics > 1:
        dimname_list = dimname_list + ['harmonic']
        dimsize_list = dimsize_list + [n_harmonics]
        # {{{ generate a grid of labels and mask out the ones we want
        harmonic_axes = array([[(1,0),(2,0),(3,0),(4,0),(5,0)],
            [(1,90),(2,90),(3,90),(4,90),(5,90)]],
            dtype=[('harmonic','int'),('phase','int')])
        harmonic_axes = harmonic_axes[harmonics]
        print "I found multiple harmonics, and am loading them into the 'harmonics' axis.  This is experimental.  You most likely will want to select the 0th element of the harmonic axis."
        # }}}
    if y_points>1:
        raise ValueError(strm("I looks like this is a 2D file len(data)/x_points/n_harmonics=",y_points,", which I"
            " just haven't bothered to program yet -- the"
            " code here is just copied from the WinEPR"
            " reader"))
        if y_points != v['REY']:
            raise CustomError('I thought REY was the indirect dim, guess not')
        if dimname=='':
            dimname = v['JEY']
        dimname_list = [dimname] + dimname_list
        dimsize_list = [y_points] + dimsize_list
    data = nddata(data,dimsize_list,dimname_list)
    if len(set(data.dimlabels)^{'harmonic'})>1:
        raise ValueError("Code below is just copied from winepr"
                " -- need to update")
        #yaxis = r_[0:v['REY']]
        #if dimname == 'mw-power-sweep':
        #    yaxis *= v['MPS']
        #    yaxis += v['XYLB'] # the starting attenuation
        #    yaxis = 10**(-yaxis/10.) # convert to linear power
        #    yaxis *= v['MP']/yaxis[0] # the initial power
        #    yaxis *= 1e-3 # convert from mW to W
        #    data.rename('mw-power-sweep','power')
        #    dimname = 'power'
        #data.labels([dimname,b0_texstr],[yaxis,x_axis])
        #data.reorder([b0_texstr,dimname])
    else:
        data.labels([b0_texstr],[x_axis])
    # }}}
    # {{{ use the parameters to rescale the data
    logger.info("There is a parameter called DModGain as well as"
            " Gain -- not sure what that is")
    rg = v.pop('Gain')
    if isscalar(rg) or rg[1] != 'dB':
        raise ValueError(strm("The gain from the file is not given in"
                " dB -- not sure what's up.  I get",rg,"for gain"))
    #data /= 10**(rg[0]/10.0)
    #data /= modulation
    # here, for winepr, I divided by the number of scans, but I'm
    # fairly sure I don't wan to do that
    # }}}
    data.other_info.update(v)
    data.reorder(b0_texstr)
    return data
def winepr(filename, dimname=''):
    """For opening WinEPR files.
    
    Parameters
    ----------
    filename : str
        The filename that ends with either ``.par`` or ``.spc``.
    """
    # {{{ determine the pair of filenames that we need
    filename = filename[:-4]+filename[-4:].upper()# case insensitive extension
    if filename[-4:] == '.SPC':
        filename_spc,filename_par = filename,filename.replace('.SPC','.PAR')
    elif filename[-4:] == '.PAR':
        filename_spc,filename_par = filename.replace('.PAR','.SPC'),filename
    else:
        raise ValueError(strm("When guessing that the filename is a"
                " WinEPR file, the extension must be either .SPC or"
                " .PAR\n"
                "This one is called",repr(filename)))
    # {{{ check if the extension is upper or lowercase
    if not os.path.exists(filename_spc):
        filename_spc = filename_spc[:-4] + filename_spc[-4:].lower()
        filename_par = filename_par[:-4] + filename_par[-4:].lower()
    # }}}
    # }}}
    # {{{ load the data
    with open(filename_spc,'rb') as fp:
        data = fp.read()
    data = fromstring(data,'<f4')
    # }}}
    # load the parameters
    v = winepr_load_acqu(filename_par)
    # {{{ use the parameters to determine the axes
    xpoints = v['RES']
    ypoints = len(data)/xpoints
    if ypoints>1:
        if ypoints != v['REY']:
            raise CustomError('I thought REY was the indirect dim, guess not')
        if dimname=='':
            dimname = v['JEY']
        data = nddata(data,[ypoints,xpoints],[dimname,b0_texstr])
    else:
        data = nddata(data,[xpoints],[b0_texstr])
    xlabels = linspace(v['HCF']-v['HSW']/2.,v['HCF']+v['HSW']/2.,xpoints)
    if len(data.dimlabels)>1:
        yaxis = r_[0:v['REY']]
        if dimname == 'mw-power-sweep':
            yaxis *= v['MPS']
            yaxis += v['XYLB'] # the starting attenuation
            yaxis = 10**(-yaxis/10.) # convert to linear power
            yaxis *= v['MP']/yaxis[0] # the initial power
            yaxis *= 1e-3 # convert from mW to W
            data.rename('mw-power-sweep','power')
            dimname = 'power'
        data.labels([dimname,b0_texstr],[yaxis,xlabels])
        data.reorder([b0_texstr,dimname])
    else:
        data.labels([b0_texstr],[xlabels])
    # }}}
    # {{{ use the parameters to rescale the data
    rg = v['RRG']
    data /= rg
    modulation = v['RMA']
    #data /= modulation
    try:
        data /= v['JNS'] # divide by number of scans
    except:
        pass
    #data /= v['MP'] # divide by power <-- weird, don't do this!
    # }}}
    data.other_info.update(v)
    return data
def winepr_load_acqu(filename):
    "Load the parameters for the winepr filename"
    with open(filename,'rU') as fp:# the U automatically converts dos format
        lines = fp.readlines()
    vars = {}
    line_re = re.compile(r'([_A-Za-z0-9]+) +(.*)')
    lines = map(string.rstrip,lines)
    #lines = [j.rstrip('\n') for j in lines] # because it's just \n, even on windows
    v = {'DRS':4096,'RES':1024,'HSW':50}
    for line in lines:
        m = line_re.match(line)
        if m is None:
            raise RuntimeError('Warning:',lsafen(repr(line)),'does not appear to be a valid WinEPR format line, and I suspect this is a problem with the terminators!')
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
    masks = [ 0x00000001L, 0x00000002L, 0x00000004L, 0x00000008L, 0x00000010L,
            0x00000020L, 0x00000040L, 0x00000080L, 0x00000100L, 0x00000200L,
            0x00000400L, 0x00000800L, 0x00001000L]
    values = map((lambda x: x&jss),masks)
    values = map(bool,values)
    values = map(bool,values)
    v.update(dict(zip(parameters,values)))
    return v
def xepr_load_acqu(filename):
    '''Load the Xepr acquisition parameter file, which should be a .dsc extension.

    Returns
    -------
    A dictionary of the relevant results.
    Because of the format of the .dsc files, this is a dictionary of
    dictionaries, where the top-level keys are the hash-block (*i.e.*
    ``#DESC``, *etc.*).
    '''
    def auto_string_convert(x):
        '''genfromtxt is from numpy -- with dtype=None, it does
        automatic type conversion -- note that strings with
        spaces will be returned as a record array it appears to
        need this StringIO function rather than a string because
        it's designed to read directly from a file.  The tolist
        converts the record array to a list.'''
        if len(x):
            try:
                return genfromtxt(StringIO(x),dtype=None).tolist()
            except:
                raise ValueError("genfromtxt chokes on "+repr(x))
        else:
            return None
    which_block = None
    block_re = re.compile(r'^ *#(\w+)')
    comment_re = re.compile(r'^ *\*')
    variable_re = re.compile(r'^ *([^\s]*)\s+(.*?) *$')
    comma_re = re.compile(r'\s*,\s*')
    with open(filename,'r') as fp:
        blocks = {}
        # {{{ read lines and assign to the appropriate block
        for line in fp:
            m = comment_re.search(line)
            if m:
                pass
            else:
                m = block_re.search(line)
                if m:
                    if which_block is not None:
                        blocks.update({which_block:dict(block_list)})
                    which_block = m.groups()[0]
                    block_list = []
                else:
                    if which_block is None:
                        raise ValueError("Appears to be stuff outside the first hashed block which, as far as I know, should not be allowed.  The first non-comment line I see is: "+repr(line))
                    else:
                        m = variable_re.search(line)
                        if m:
                            if ',' in m.groups()[1]:
                                # {{{ break into lists
                                block_list.append((m.groups()[0],
                                        map(auto_string_convert,
                                            comma_re.split(
                                                m.groups()[1]))))
                                # }}}
                            else:
                                block_list.append((m.groups()[0],
                                        auto_string_convert(
                                            m.groups()[1])))
                        else:
                            raise ValueError("I don't know what to do with the line:\n"+line)
        blocks.update({which_block:dict(block_list)})
        # }}}
    return blocks
