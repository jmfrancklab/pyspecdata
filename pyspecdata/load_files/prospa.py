"routines specific to loading information from prospa files"
from ..core import *
def decim_correct(data):
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
    print('Applied prospa decimation correction')
    return data
def load_acqu(file):
    file = dirformat(file)
    fp = open(file+'acqu.par')
    lines = fp.readlines()
    line_re = re.compile(r'([^ \t]+) *= *(.+)')
    vars = {}
    for j in range(0,len(lines)):
        lines[j] = string.rstrip(lines[j])
        m = line_re.match(lines[j])
        if m:
            exec('temp = %s'%m.groups()[1])
            vars.update({m.groups()[0]:temp})
        else:
            print("error, acqu.par line not parsed: ",lines[j])
    fp.close()
    return vars
def t1_info(file):
    file = dirformat(file)
    if det_type(file) == ('prospa','t1_sub'):
        file += '../'
    elif not det_type(file) == ('prospa','t1'):
        raise CustomError("You're trying to get prospa T1 info from a file that's not a Prospa T1 file!")
    files = [x for x in os.listdir(file) if os.path.isdir(dirformat(file)+x)]
    file_re = re.compile(r'([0-9]+)msDelay$')
    datafiles = []
    wait_time = []
    print('DEBUG: prospa is searching for times in the file list',files)
    for j in range(0,len(files)):
        m = file_re.match(files[j])
        if m:
            datafiles += [file+files[j]]
            wait_time += [int(m.groups()[0])]
    return datafiles,array(wait_time)*1e-3
def load_datafile(file,dims=1):
    r'''load a prospa datafile into a flat array as a 1D file
    use dims=2 if it's a 2D file'''
    file = dirformat(file)
    if dims == 1:
        fp = open(file+'data.1d','rb')
    elif dims == 2:
        fp = open(file+'data.2d','rb')
    else:
        print('ERROR: wrong number of dims')
    data = fp.read()
    data = array(struct.unpack('%df'%(len(data)//4),data))
    data = data[7:]
    # the following is junk!!!
    #elif precision=='b':
    #   data = array(struct.unpack('%db'%(len(data)//1),data))
    #   data = data[7*4:]
    #else:
    #   print 'error, precision wrong'
    data = reshape(data,(-1,2))
    data = data[:,0]+1j*data[:,1]
    fp.close()
    return data
def load_1D(filename):
    v = prospa_load_acqu(filename)
    data = prospa_load_datafile(filename,dims=1)/v['nrScans']#added this 2/20/13 to allow automatic signal averaging
    data = nddata(data,[v['nrPnts']],['t2'])
    taxis = linspace(0,1,v['nrPnts'])*v['acqTime']/1e3
    data.labels(['t2'],[taxis])
    return data
def load_2D(filename, dimname=''):
    #{{{ Prospa 2D
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
    return data
