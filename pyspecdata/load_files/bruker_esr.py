from ..core import *
from ..general_functions import strm
from ..datadir import rclone_search
import numpy as np
import re
from io import StringIO
import os
logger = logging.getLogger('pyspecdata.load_files.bruker_esr')
b0_texstr = r'$B_0$'
def xepr(filename, exp_type=None, dimname='', verbose=False):
    """For opening Xepr files.
    
    Parameters
    ----------
    filename : str
        The filename that ends with ``.DSC``, ``.DTA``, or ``.YGF``.
    """
    # {{{ determine the pair of filenames that we need
    filename = filename[:-4]+filename[-4:].upper()# case insensitive extension
    if filename[-4:] == '.DTA':
        filename_spc,filename_par = filename,filename.replace('.DTA','.DSC')
    elif filename[-4:] == '.DSC':
        filename_spc,filename_par = filename.replace('.DSC','.DTA'),filename
    elif filename[-4:] == '.YGF':
        filename_spc,filename_par = filename.replace('.YGF','.DSC'),filename.replace('.YGF','.DTA')
    else:
        raise ValueError(strm("When guessing that the filename is a"
                " WinEPR file, the extension must be either .SPC or"
                " .PAR\n"
                "This one is called",repr(filename)))
    # {{{ check if the extension is upper or lowercase
    orig_spc = filename_spc
    if not os.path.exists(filename_spc):
        logging.debug(f"{filename_spc} doesn't exist -- trying lowercase")
        filename_spc = filename_spc[:-4] + filename_spc[-4:].lower()
    if not os.path.exists(filename_spc):
        filename_spc = orig_spc
    orig_par = filename_par
    if not os.path.exists(filename_par):
        logging.debug(f"{filename_par} doesn't exist -- trying lowercase")
        filename_par = filename_par[:-4] + filename_par[-4:].lower()
    if not os.path.exists(filename_par):
        filename_par = orig_par
    # }}}
    # }}}
    # {{{ load the parameters
    v = xepr_load_acqu(filename_par)
    # {{{ flatten the dictionary (remove the uppermost/block
    #     level)
    new_v = {}
    for k_a,v_a in v.items():
        new_v.update(v_a)
    v = new_v
    # }}}
    ikkf = v['IKKF']
    if type(ikkf) == str: ikkf = [ikkf]
        
    # }}}
    # {{{ load the data
    if not os.path.exists(filename_spc):
        # because the spc isn't part of the original search, we
        # need to log the fact that it's missing manually
        if exp_type is None:
            raise ValueError("I could probably find "
            "your file remotely, but you called with "
            "exp_type None!")
        rclone_search(
                os.path.split(filename_spc)[-1],
                exp_type,
                os.path.split(filename_spc)[0])
    with open(filename_spc,'rb') as fp:
        if all([j == 'REAL' for j in ikkf]):
            data = np.frombuffer(fp.read(),'>f8')
        elif all([j == 'CPLX' for j in ikkf]):
            data = np.frombuffer(fp.read(),'>c16')
        else:
            raise ValueError('the data type (IKKF) is givn as '
                    +' '.join(ikkf)
                    +" and I don't support mixed types!")
    # }}}
    # {{{ use the parameters to determine the axes
    #     pop parameters that are just part of the axes
    x_points = v.pop('XPTS')
    x_axis = r_[0:x_points]
    # the following is NOT the same as *=, which preserves the
    # type (=bad!)!!!
    x_axis = x_axis*v.pop('XWID')/x_axis[-1] # wouldn't be a huge
    #         difference, but I'm not sure if this is correct
    #         (I think so)
    x_axis += v.pop('XMIN') # actually using pop is better than calling these, so that we don't have redundant information
    harmonics = np.array([[False] * 2]*5) # inner dimension for the 90 degree phase
    for j,jval in enumerate(['1st','2nd','3rd','4th','5th']):
        for k,kval in enumerate(['','90']):
            thiskey = 'Enable'+jval+'Harm'+kval
            if thiskey in v.keys() and v[thiskey]:
                harmonics[j,k] = True
    n_harmonics = np.sum(harmonics)
    if n_harmonics == 0: n_harmonics = 1
    logger.debug('there are %d harmonics'%n_harmonics)
    if n_harmonics > 1:
        logger.debug('there are %d harmonics, first is of type %s'%(n_harmonics,ikkf[0]))
    # }}}
    # {{{ check that calculated axes match dimensions
    y_points_calcd = len(data)//x_points//n_harmonics
    if 'XNAM' in list(v.keys()):
        x_dimname = v['XNAM']
    else:
        x_dimname = b0_texstr
    dimname_list = [x_dimname]
    dimsize_list = [x_points]
    dims_accounted_for = {x_dimname}
    dims_to_label = {x_dimname:x_axis}
    def interpret_units(un_key):
        retval = v.pop(un_key)
        if retval[0] == "'": retval = retval.replace("'","")
        return retval
    if 'XUNI' in list(v.keys()):
        dim_units = {x_dimname:interpret_units('XUNI')}
    else:
        dim_units = {}
    if n_harmonics > 1:
        dimname_list = dimname_list + ['harmonic']
        dimsize_list = dimsize_list + [n_harmonics]
        # {{{ generate a grid of labels and mask out the ones we want
        #harmonic_axes = np.array([[(1,0),(2,0),(3,0),(4,0),(5,0)],
        #    [(1,90),(2,90),(3,90),(4,90),(5,90)]],
        #    dtype=[('harmonic','int'),('phase','int')])
        harmonic_axes = np.array([(1,0),(2,0),(3,0),(4,0),
            (5,0),(1,90),(2,90),(3,90),(4,90),(5,90)],
            dtype=[('harmonic','int'),('phase','int')])
        harmonic_axes = harmonic_axes.reshape(2,5).T
        harmonic_axes = harmonic_axes[harmonics]
        dims_to_label.update({'harmonic':harmonic_axes})
        logger.debug("I found multiple harmonics, and am loading them into the 'harmonics' axis.  This is experimental.  You most likely will want to select the 0th element of the harmonic axis.")
        dims_accounted_for |= {'harmonic'}
        # }}}
    y_dim_name = None
    if y_points_calcd>1:
        if 'YPTS' in list(v.keys()):
            assert v['YPTS']==y_points_calcd, ("y points left over after"
                    " reshaping according to harmonics (y_points_calcd="
                    + str(y_points_calcd)
                    + ") doesn't match the number of data points")
            assert 'YTYP' in list(v.keys()), ("No parameter YTYP -- how do you expect me to know the type of 2D dataset??")
            if v['YTYP'] == 'IGD':
                logger.debug('Found YTYP=IGD, assuming this is a power series')
                assert 'YNAM' in list(v.keys()), ("No parameter YNAM -- how do you expect me to know the name of the second dimension??")
                y_dim_name = v.pop('YNAM')
                if isinstance(y_dim_name, list):
                    y_dim_name = ' '.join(y_dim_name) # it gets split into a list, which for XEpr files shouldn't be happening, but fix later
                    if y_dim_name[0] == "'": y_dim_name = y_dim_name.replace("'","")
                if y_dim_name.startswith("'") and y_dim_name.endswith("'"):
                    y_dim_name = y_dim_name[1:-1]
                assert 'YUNI' in list(v.keys()), ("No parameter YUNI -- how do you expect me to know the units of the second dimension??")
                dim_units.update({y_dim_name:interpret_units('YUNI')})
                filename_ygf = filename_par[:-4] + '.YGF'
                for j in range(2):
                    if not os.path.exists(filename_ygf) and j>0:
                        filename_ygf = filename_ygf[:-4] + filename_ygf[-4:].lower()
                    if not os.path.exists(filename_ygf):
                        # because the spc isn't part of the original search, we
                        # need to log the fact that it's missing manually
                        if exp_type is None:
                            raise ValueError("I could probably find "
                            "your file remotely, but you called with "
                            "exp_type None!")
                        rclone_search(
                                os.path.split(filename_ygf)[-1],
                                exp_type,
                                os.path.split(filename_ygf)[0])
                with open(filename_ygf,'rb') as fp:
                    y_axis = fp.read()
                y_axis = np.fromstring(y_axis,'>f8')
                assert len(y_axis)==y_points_calcd, "Length of the power axis doesn't seem to match!"
            else:
                raise ValueError(strm("found YTYP=",v['YTYP']," which is not currently programmed"))
        else:
            logger.debug("y_points_calcd is greater than 1, but YPTS is not set, so assuming WinEPR style")
            assert y_points_calcd == v['REY'], 'Trying WinEPR style and I thought REY was the indirect dim, guess not'
            if dimname=='':
                y_dim_name = v['JEY']
            dimsize_list = [y_points_calcd] + dimsize_list
            y_axis = r_[0:v['REY']]
            y_axis *= v['MPS']
            y_axis += v['XYLB'] # the starting attenuation
            y_axis = 10**(-yaxis/10.) # convert to linear power
            y_axis *= v['MP']/yaxis[0] # the initial power
            y_axis *= 1e-3 # convert from mW to W
        dimsize_list = [y_points_calcd] + dimsize_list
        dimname_list = [y_dim_name] + dimname_list
        dims_accounted_for |= {y_dim_name}
        dims_to_label.update({y_dim_name:y_axis})
    data = nddata(data,dimsize_list,dimname_list)
    extra_dims = set(data.dimlabels) - dims_accounted_for
    assert len(extra_dims)==0, strm("You seem to have one or more extra dimension(s)"
            "called",str(extra_dims),
            "(Code below is just copied from winepr"
            " -- need to update)")
    data.labels(dims_to_label)
    for k,val in dim_units.items():
        data.set_units(k,val)
    # }}}
    data.other_info.update(v)
    data.reorder(x_dimname)
    if 'Microwave Power' in data.dimlabels:
        # convert from W to dBm
        data.setaxis('Microwave Power',lambda x: 10*np.log10(x)).set_units('Microwave Power','dBm')
    if data.get_units(x_dimname) == 's': data.rename(x_dimname,'t')
    # {{{ this is a hack -- should be replaced with pint, etc, but needs
    # more thought for a general solution!
    for thisdim in data.dimlabels:
        if data.get_units(thisdim) == 'ns':
            data[thisdim] /= 1e9
            data.set_units(thisdim,'s')
    if 'Field' in data.dimlabels:
        data.rename('Field',b0_texstr)
    # }}}
    return data
def winepr(filename, dimname='', exp_type=None):
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
    if not os.path.exists(filename_spc):
        if exp_type is None:
            raise ValueError("I could probably find "
            "your file remotely, but you called with "
            "exp_type None!")
        rclone_search(
                os.path.split(filename_spc)[0],
                exp_type,
                os.path.split(filename_spc)[-1])
    with open(filename_spc,'rb') as fp:
        data = fp.read()
    data = np.fromstring(data,'<f4')
    # }}}
    # load the parameters
    v = winepr_load_acqu(filename_par)
    # {{{ use the parameters to determine the axes
    xpoints = v['RES']
    ypoints = len(data)/xpoints
    if ypoints>1:
        if ypoints != v['REY']:
            raise RuntimeError('I thought REY was the indirect dim, guess not')
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
    lines = list(map(lambda x: x.rstrip(),lines))
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
    jss = int(v['JSS'])
    parameters = [ 'DUAL', '2D', 'FT', 'MAN0', 'MAN1', 'PROT', 'VEPR', 'POW', 'ABS', 'FTX', 'FTY', 'POW2', 'ABS2']
    parameters = list(map((lambda x: 's_'+x),parameters))
    masks = [ 0x00000001, 0x00000002, 0x00000004, 0x00000008, 0x00000010,
            0x00000020, 0x00000040, 0x00000080, 0x00000100, 0x00000200,
            0x00000400, 0x00000800, 0x00001000]
    values = list(map((lambda x: x&jss),masks))
    values = list(map(bool,values))
    values = list(map(bool,values))
    v.update(dict(list(zip(parameters,values))))
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
                result = np.genfromtxt(StringIO(x),dtype=None,encoding='utf-8').tolist()
            except:
                raise ValueError("genfromtxt chokes on "+repr(x))
            if type(result) is str and result[0] == "'" and result[-1] == "'":
                result = result[1:-1]
            return result
        else:
            return None
    which_block = None
    block_re = re.compile(r'^ *#(\w+)')
    comment_re = re.compile(r'^ *\*')
    variable_re = re.compile(r'^ *([^\s]*)\s+(.*?) *$')
    comma_re = re.compile(r'\s*,\s*')
    with open(filename,'r',encoding='utf-8') as fp:
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
                                        list(map(auto_string_convert,
                                            comma_re.split(
                                                m.groups()[1])))))
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
