from ..core import *
from numpy import fromstring
def winepr(filename, dimname=''):
    """For opening WinEPR files.
    
    Parameters
    ----------
    filename : str
        The filename that ends with either ``.par`` or ``.spc``.
    """
    # {{{ determine the pair of filenames that we need
    filename = filename[:-4]+filename[-4:].upper()# case insensitive extension
    if filename[:-4] == '.SPC':
        filename_spc,filename_par = filename,filename.replace('.SPC','.PAR')
    elif filename[:-4] == '.PAR':
        filename_spc,filename_par = filename.replace('.PAR','.SPC'),filename
    else:
        raise ValueError("When guessing that the filename is a WinEPR file, the"
                " extension must be either .SPC or .PAR")
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
    # {{{ use the parameters to rescale the data and determine the axes
    xpoints = v['RES']
    rg = v['RRG']
    data /= rg
    modulation = v['RMA']
    #data /= modulation
    try:
        data /= v['JNS'] # divide by number of scans
    except:
        pass
    #data /= v['MP'] # divide by power <-- weird, don't do this!
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
        yaxis = r_[0:v['REY']]
        if dimname == 'mw-power-sweep':
            yaxis *= v['MPS']
            yaxis += v['XYLB'] # the starting attenuation
            yaxis = 10**(-yaxis/10.) # convert to linear power
            yaxis *= v['MP']/yaxis[0] # the initial power
            yaxis *= 1e-3 # convert from mW to W
            data.rename('mw-power-sweep','power')
            dimname = 'power'
        data.labels([dimname,b0],[yaxis,xlabels])
        data.reorder([b0,dimname])
    else:
        data.labels([b0],[xlabels])
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
