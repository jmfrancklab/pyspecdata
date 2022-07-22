""" written based on igor code
(under GPL license)
taken from `here <https://github.com/Yohko/importtool/blob/master/Igor%20Procedures/import%20tool/import_VARIAN_CARY_BIN.ipf>`_.

"""
from numpy import *
from ..general_functions import strm
from ..core import nddata
import os, logging

def load_bindata(fp,param):
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
    #for j in range(500):
    #    fp.seek(current_pos+j,0)
    #    result = fromfile(fp, dtype='<f4', count=1)
    #    if result>1 and result<1000:
    #        logging.debug(strm("try",j,"offset",result))
    retval = fromfile(fp, dtype=[('x','<f4'),('y','<f4')], count=param['points'])
    #logging.debug(retval)
    retval = nddata(retval['y'],'$\lambda$' ).setaxis('$\lambda$',retval['x']
            ).set_units('$\lambda$',x_unit)#.set_units(y_unit)
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
        names = []
        while fp.tell() < file_end-4:
            param['blockoffset'] = fp.tell()
            logging.debug(strm("blockoffset",param['blockoffset']))
            # line 100
            thislen = fromfile(fp,dtype='<u4', count=1).item()
            temp = fromfile(fp,dtype=f'a{thislen}', count=1).item()
            try:
                param['Tstore_type'] = temp.decode('ascii')
            except:
                raise IOError(strm("problem decoding or assigning Tstore type",repr(temp)))
            logging.debug("Tstore_type \"%s\""%param['Tstore_type'])
            # line 102
            param['blocklen'] = fromfile(fp,dtype='<u4', count=1).item()
            if param['Tstore_type'] == "TContinuumStore":
                # line 106
                param = load_header(fp, param)
                logging.debug(strm("param:",param))
                data = load_bindata(fp, param)
                #data.name(param['spectrum_name'])
                names.append(param['spectrum_name'])
                alldata.append(data)

                logging.debug(strm("at this point, we are at",
                        fp.tell(),"relative to end of block",
                        param['blockoffset']+param['blocklen'],
                        "and end of file",file_end))
                fp.seek(param['blockoffset']+param['blocklen'], 0)
            else:
                logging.info(strm("skipping block of type",param['Tstore_type']))
                fp.seek(param['blockoffset']+param['blocklen'], 0)
                #raise ValueError(strm("not yet set up for Tstore_type",param['Tstore_type']))
            #alldata.append(data)
    retval = {}
    for n,j in enumerate(alldata):
        repcounter = 0
        #orig_name = j.name()
        if n < len(names):
            orig_name = names[n]
        else:
            orig_name = names[len(names)-1]

        logging.debug("orig name for %d is %s"%(n,orig_name))
        new_name = orig_name
        while new_name in retval.keys():
            repcounter += 1
            new_name = orig_name + '_rep%03d'%repcounter
        #if repcounter > 0:
            #logging.warn("You have a duplicate spectrum name!!! -- renamed it from %s to %s"%(orig_name,new_name))
        retval[new_name] = j
        #retval[new_name].name(new_name) # this messes with the keys -- not sure why?
    return retval
