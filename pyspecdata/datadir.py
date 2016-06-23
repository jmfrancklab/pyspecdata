import os
import ConfigParser
import platform
from .general_functions import process_kwargs
# {{{ http://stackoverflow.com/questions/2533120/show-default-value-for-editing-on-python-input-possible
if platform.platform().startswith('Windows'):
    import win32console
    _stdin = win32console.GetStdHandle(win32console.STD_INPUT_HANDLE)
    def rlinput(prompt, default=''):
        keys = []
        for c in unicode(default):
            evt = win32console.PyINPUT_RECORDType(win32console.KEY_EVENT)
            evt.Char = c
            evt.RepeatCount = 1
            evt.KeyDown = True
            keys.append(evt)

        _stdin.WriteConsoleInput(keys)
        return raw_input(prompt)
else:
    def rlinput(prompt, prefill=''):
        import readline
        readline.set_startup_hook(lambda: readline.insert_text(prefill))
        try:
            return raw_input(prompt)
        finally:
            readline.set_startup_hook()
# }}}
class MyConfig(object):
    def __init__(self):
        self._config_parser = None
        if platform.platform().startswith('Windows'):
            self.hide_start = '_' # the default hidden/config starter for vim, mingw applications, etc
            "This filename prefix denotes a configuration file on the OS."
        else:
            self.hide_start = '.'
        self.config_location = os.path.join(os.path.expanduser('~'),self.hide_start+'pyspecdata')
        self.config_vars = {}
        "The dictionary that stores the current settings -- keep these in a dictionary, which should be faster than reading from environ, or from a file."
        return
    def set_setting(self,this_section,this_key,this_value):
        "set `this_key` to `this_value` inside section `this_section`, creating it if necessary"
        if self._config_parser is None:
            self._config_parser = ConfigParser.ConfigParser()
            read_cfg = self._config_parser.read(self.config_location)
            if not read_cfg:
                print "\nWarning!! There was no file at",self.config_location,"so I'm creating one"
        if not self._config_parser.has_section(section):
            self._config_parser.add_section(section)
        self._config_parser.set(this_section,this_key,this_value)
    def __exit__(self):
        self.__del__()
    def __del__(self):
        with open(self.config_location,'w') as fp:
            self._config_parser.write(fp)
    def get_setting(self,this_key,environ = None,default = None,section = 'General'):
        """Get a settings from the "General" group.
        If the file does not exist, or the option is not set, then set the option, creating the file as needed.
        The option is determined in the following order:

            * The value in the `config_vars` dictionary.
            * The value of the environment variable named `environ`.
            * The value stored in the configuration file at ``~/.pyspecdata`` (``~/_pyspecdata`` on Windows).
            * The value given by the `default` argument.  If `default` is ``None``, then return ``None``.

        Parameters
        ----------
        this_key : str
            The name of the settings that we want to retrieve.
        environ : str
            If the value corresponding to `this_key` is not present in the
            `self.config_vars` dictionary, look for an environment variable
            called `environ` that stores the value of `this_key`.  If this is
            not set, it's set from the config file or the default argument.
        default : str
            If the value for `this_key` is not in
            `self.config_vars`, not in the environment variable
            called `environ`, and not in the config file, set it
            to this value, and ask for user response to confirm.
            Then, set the config file value, the environment, and the
            `self.config_vars` value to the result.

            For platform compatibility, leading period characters are converted
            to `self.hide_start`, and '/' are converted to `os.path.sep`.
        section : str
            Look in this section of the config file.
            
        Returns
        -------
        The value corresponding to `this_key`.
        """
        if this_key in self.config_vars.keys():
            return self.config_vars[this_key]
        if environ is not None and environ in os.environ.keys():
                retval = os.environ[environ]
        else:
            if self._config_parser is None:
                self._config_parser = ConfigParser.ConfigParser()
                read_cfg = self._config_parser.read(self.config_location)
                if not read_cfg:
                    print "\nWarning!! There was no file at",self.config_location,"so I'm creating one"
            if self._config_parser.has_section(section):
                try:
                    retval = self._config_parser.get(section,this_key)
                except NoOptionError:
                    retval = None
            else:
                self._config_parser.add_section(section)
                retval = None
            if retval is None:# it wasn't found from the config file
                if default is None:
                    return None
                default = default.split('/')
                for j in range(len(default)):
                    if default[j][0] == '.':
                        default[j] = self.hide_start + default[j][1:]
                default = os.path.sep.join(default)
                retval = rlinput("\nI didn't find the value corresponding to "+this_key
                        +" in the environment variable "+repr(environ)
                        +" or in the config file, so confirm what value you would like here:\t",default)
                if '~' in retval:
                    retval = os.path.expanduser(retval)
                    retval = rlinput("\nI detected a tilde (~) character -- check that this looks like what you meant and edit it or (more likely) just press enter:\t",retval)
                self._config_parser.set(section,this_key,retval)
            if environ is not None:
                os.environ[environ] = retval
        self.config_vars[this_key] = retval
        return retval
_my_config = MyConfig()
def get_notebook_dir(*arg):
    r'''Returns the notebook directory.  If arguments are passed, it returns the directory underneath the notebook directory, ending in a trailing (back)slash
    
    It is determined by a call to `MyConfig.get_setting` with the environment variable set to ``PYTHON_NOTEBOOK_DIR`` and default ``~/notebook``.
    '''
    base_notebook_dir = _my_config.get_setting('notebook_directory',environ = 'PYTHON_NOTEBOOK_DIR',default = '~/notebook')
    retval = (base_notebook_dir,) + arg
    if len(retval[-1]) != 0:
        retval = retval + ('',)
    return os.path.join(*retval)
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
def getDATADIR(*args,**kwargs):
    r'''Returns the base directory where you put all your data.  If arguments
    are passed, it returns the directory underneath the data directory, ending
    in a trailing (back)slash
    
    It is determined by a call to `MyConfig.get_setting` with the environment
    variable set to ``PYTHON_DATA_DIR`` and default ``~/experimental_data``.

    Parameters
    ----------
    exp_type : str
        A string identifying the experiment type.
        It searches for `exp_type` in this order:

        * Look in the ``ExpTypes`` section of the config file.
        * use `os.walk` to search for a directory with this name,
            excluding things that start with '.', '_' or
            containing '.hfssresults', always choosing the
            thing that's highest up in the tree.
    '''
    exp_type = process_kwargs(('exp_type',None),kwargs)
    base_data_dir = _my_config.get_setting('data_directory',environ = 'PYTHON_DATA_DIR',default = '~/experimental_data')
    if exp_type is not None:
        # {{{ determine the experiment subdirectory
        exp_directory = _my_config.get_setting(exp_type,section = 'ExpTypes')
        if exp_directory is None:
            equal_matches = []
            containing_matches = []
            for d,s,_ in os.walk(base_data_dir):
                s[:] = [j for j in s if j[0]
                        not in ['.','_']
                        and '.hfssresults' not in j]# prune the walk
                equal_matches.extend([os.path.join(d,j) for j in s if j == exp_type])
                containing_matches.extend([os.path.join(d,j) for j in s if exp_type.lower() in j.lower()])
            def grab_smallest(matches):
                if len(matches) > 0:
                    if len(matches) == 1:
                        return matches[0]
                    else:
                        min_length_match = min(map(len,matches))
                        matches = filter(lambda x: len(x) == min_length_match,matches)
                        if len(matches) != 1:
                            raise ValueError("I found multiple equivalent matches when searching for exp_type: "+repr(matches))
                        return matches[0]
            if len(equal_matches) > 0:
                exp_directory = grab_smallest(equal_matches)
            elif len(containing_matches) > 0:
                exp_directory = grab_smallest(containing_matches)
            else:
                raise ValueError("I found no directory matches for exp_type "+exp_type)
            _my_config.set_setting('ExpTypes',exp_type,exp_directory)
        retval = (exp_directory,) + arg
        # }}}
    else:
        retval = (base_data_dir,) + arg
    if len(retval[-1]) != 0:
        retval = retval + ('',)
    return os.path.join(*retval)
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
