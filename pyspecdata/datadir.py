import os
def dirformat(file):
        #{{{ format strings
        if file[-1]!='/':
            file += '/'
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
        datadir_error = ('''Since this is your first time running the notebook code, you need to create a file called .datadir in the directory

                %s

                (your base notebook directory) that tells where all your data is stored
                
                '''%os.getcwd()).replace('\\','\\textbackslash ').replace('_','\\_')
    return grabbed_datadir_from_file,datadir_error
def getDATADIR():
    if 'PYTHONDATADIR' in os.environ.keys():
        DATADIR = os.environ['PYTHONDATADIR']
    else:
        grabbed_datadir_from_file,datadir_error = grab_data_directory()
        if datadir_error is not False:
            raise RuntimeError(datadir_error)
        DATADIR = os.environ['PYTHONDATADIR']
    return dirformat(DATADIR)
