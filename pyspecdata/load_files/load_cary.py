""" written based on igor code
(under GPL license)
taken from `here <https://github.com/Yohko/importtool/blob/master/Igor%20Procedures/import%20tool/import_VARIAN_CARY_BIN.ipf>`_.

"""
from numpy import *
from ..general_functions import strm
from ..core import nddata
from ..placeholder import nddata_placeholder
import os, logging

def load_bindata(fp,param,filename):
    x_mode_rep = "nm;Å;cm-1;°".split(';')
    y_mode_rep = "Abs;%T;Absorptivity;%R;4?;5?;Log(Abs);Absolute %R;F(R);Log(F(R));Log(1/R)".split(';')
    factoreV_list = [1239.8424,1239.8424,12398.424,8065.54429,0]
    x_unit = x_mode_rep[int(param['x_mode'])]
    factoreV = factoreV_list[int(param['x_mode'])]
    y_unit = y_mode_rep[int(param['y_mode'])]
    # I'm not creating the spec_counter variable that is is
    #tmps = str(param['spec_counter'])+"_"+param['spectrum_name']+"_"+x_unit
    logging.debug(strm("about to read",param['points'],"double points, with x units",x_unit,"and y unit",y_unit))
    current_pos = fp.tell()
    logging.debug(strm("inside bindata, position is",fp.tell()))
    num_points = param['points']
    #logging.debug(retval)
    def retrieve_data():
        with open(filename,'rb') as innerfp: # b/c file will be closed by the time we try to access the data
            innerfp.seek(current_pos,0)
            data_struct_array = fromfile(innerfp, dtype=[('x','<f4'),('y','<f4')], count=num_points)
        def followup(self):
            self.setaxis('wavelength', data_struct_array['x'])
        return followup,data_struct_array['y']
    retval = nddata_placeholder(retrieve_data)
    retval.dimlabels = ['wavelength']
    retval.set_units('wavelength',x_unit).set_units(y_unit)
    return retval
def load_header(fp, param):
    retval = fromfile(fp, dtype=[
        ("softversion","<u4"),
        ("tmp1","<u4"),
        ("end_x","<i4"),
        ("start_x","<i4"),
        ("V_min","<i4"),
        ("V_max","<i4"),
        ("points","<u4"),
        ("tmp2","<i4"),
        ("tmp3","<u2"),
        ("tmp4","<u2"),
        ("tmp5","<i4"),
        ("spectrum_number","<i4"),
        ("tmp6","<f4"),
        ("tmp7","<i4"),
        ("tmp8","<f4"),
        ("tmp9","<i4"),
        ("x_mode","<f4"),
        ("y_mode","<f4")], count=1)
    if param['Tstore_type'] == "TContinuumStore":
        fp.seek(789+param['blockoffset']+4+2,0)
    else:
        raise ValueError(strm("not yet set up for Tstore_type",param['Tstore_type']))
    param.update(dict(zip(retval.dtype.fields,retval.item())))
    spectrum_name = fromfile(fp,dtype='a256', count=1).item()
    spectrum_name = spectrum_name[0:spectrum_name.find(b'\x00')].decode('ascii')
    param['spectrum_name'] = spectrum_name
    return param

def load_cary(filename):
    filesize = os.stat(filename).st_size
    logging.debug(strm("filesize",filesize))
    param = {}
    with open(filename,'rb') as fp:
        # line 92
        thislen = fromfile(fp,dtype='<u1', count=1).item()
        logging.debug(strm("thislen is",thislen))
        magic = fromfile(fp,dtype=f'a{thislen}', count=1).item()
        logging.debug(strm("magic is",magic))
        _ = fromfile(fp,dtype=f'a{61-thislen}', count=1)
        marker = fp.tell()
        fp.seek(0,2) # seek to end
        file_end = fp.tell()
        fp.seek(marker,0)
        alldata = []
        while fp.tell() < file_end-4:
            param['blockoffset'] = fp.tell()
            logging.debug(strm("blockoffset",param['blockoffset']))
            # line 100
            thislen = fromfile(fp,dtype='<u4', count=1).item()
            param['Tstore_type'] = fromfile(fp,dtype=f'a{thislen}', count=1).item().decode('ascii')
            logging.debug("Tstore_type \"%s\""%param['Tstore_type'])
            # line 102
            param['blocklen'] = fromfile(fp,dtype='<u4', count=1).item()
            if param['Tstore_type'] == "TContinuumStore":
                # line 106
                param = load_header(fp, param)
                logging.debug(strm("param:",param))
                data = load_bindata(fp, param, filename)
                data.name(param['spectrum_name'])
                logging.debug(strm("at this point, we are at",
                        fp.tell(),"relative to end of block",
                        param['blockoffset']+param['blocklen'],
                        "and end of file",file_end))
                fp.seek(param['blockoffset']+param['blocklen'], 0)
            alldata.append(data)
    retval = {}
    for j in alldata:
        while j.name() in retval.keys():
            logging.warn("You have a duplicate spectrum name!!! -- renaming it from %s to %s_rep"%(j.name(),j.name()))
            j.name(j.name()+'_rep')
        retval[j.name()] = j
    return retval
