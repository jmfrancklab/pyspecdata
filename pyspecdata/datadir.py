import os
def get_notebook_dir(*arg):
    r'returns the notebook directory -- make this settable later'
    if len(arg) == 0:
        return os.path.expanduser('~') + os.path.sep + 'notebook' + os.path.sep
    elif len(arg) > 0:
        return os.path.expanduser('~') + os.path.sep + 'notebook' + os.path.sep + os.path.sep.join(arg)
def dirformat(file):
        #{{{ format strings
        if file[-1]!=os.sep:
            file += os.sep
        #}}}
        return file
def grab_data_directory():
    'this checks if the environment variable for the python data directory, PYTHONDATADIR, is set, and that .matplotlib directory has been created\n if needed, it sets and creates them, respectively\n\t\treturns: grabbed_datadir_from_file --> it sucessfully pulled PYTHONDATADIR from the .datadir file\n\t\t\t datadir_error --> it couldn\'t figure out what the data directory was, and returns an appropriate error message'
    grabbed_datadir_from_file = False
    datadir_error = False
    if not os.path.exists(os.getcwd()+'/.matplotlib'):
        os.mkdir(os.getcwd()+'/.matplotlib')
    if 'PYTHONDATADIR' in os.environ.keys():
        pass
    elif os.path.exists('.datadir'):
        fp_datadir = open('.datadir','r')
        mydatadir = fp_datadir.read().strip()
        if mydatadir[-1] not in ['/','\\']:
            if os.name == 'posix':
                mydatadir = mydatadir+'/'
            else:
                mydatadir = mydatadir+'\\'
        os.environ['PYTHONDATADIR'] = mydatadir
        fp_datadir.close()
        grabbed_datadir_from_file = True
    else:
        mydatadir = os.path.expanduser('~') + os.path.sep + 'exp_data'
        os.environ['PYTHONDATADIR'] = mydatadir
    return grabbed_datadir_from_file,datadir_error
def getDATADIR(*args):
    r'''The data directory is set through the PYTHONDATADIR environment variable.
    This assumes that the arguments are a series of subdirectories under that directory, and returns the resulting path, with trailing backslash/slash as appropriate.'''
    if 'PYTHONDATADIR' in os.environ.keys():
        DATADIR = os.environ['PYTHONDATADIR']
    else:
        grabbed_datadir_from_file,datadir_error = grab_data_directory()
        if datadir_error is not False:
            raise RuntimeError(datadir_error)
        DATADIR = os.environ['PYTHONDATADIR']
    if len(args)>0:
        return dirformat(dirformat(DATADIR)+os.sep.join(args))
    else:
        return dirformat(DATADIR)
