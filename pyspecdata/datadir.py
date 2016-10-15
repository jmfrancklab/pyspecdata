import os
import ConfigParser
import platform
from .general_functions import process_kwargs
# {{{ http://stackoverflow.com/questions/2533120/show-default-value-for-editing-on-python-input-possible
if platform.platform().startswith('Windows'):
    import win32console
    _stdin = win32console.GetStdHandle(win32console.STD_INPUT_HANDLE)
    def rlinput_forwindows(prompt, default=''):
        keys = []
        for c in unicode(default):
            evt = win32console.PyINPUT_RECORDType(win32console.KEY_EVENT)
            evt.Char = c
            evt.RepeatCount = 1
            evt.KeyDown = True
            keys.append(evt)

        _stdin.WriteConsoleInput(keys)
        return raw_input(prompt)
    def rlinput(prompt, default=''):
        try:
            rlinput_forwindows(prompt, default='')# doesn't work from inside the newer git window
        except:
            print "My guess:",default
            result = raw_input(prompt+"\nMy guess: "+default+'\n')
            if result.strip() == '':
                result = default
            return default
else:
    def rlinput(prompt, prefill=''):
        import readline
        readline.set_startup_hook(lambda: readline.insert_text(prefill))
        try:
            return raw_input('\n'+prompt)
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
            self._config_parser = ConfigParser.SafeConfigParser()
            read_cfg = self._config_parser.read(self.config_location)
            if not read_cfg:
                print "\nWarning!! There was no file at",self.config_location,"so I'm creating one"
        if not self._config_parser.has_section(this_section):
            self._config_parser.add_section(this_section)
        self._config_parser.set(this_section,this_key,this_value)
        return
    def __exit__(self):
        self.__del__()
    def __del__(self):
        if self._config_parser is not None:
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
                self._config_parser = ConfigParser.SafeConfigParser()
                read_cfg = self._config_parser.read(self.config_location)
                if not read_cfg:
                    print "\nWarning!! There was no file at",self.config_location,"so I'm creating one"
            if self._config_parser.has_section(section):
                try:
                    retval = self._config_parser.get(section,this_key)
                except ConfigParser.NoOptionError:
                    retval = None
            else:
                self._config_parser.add_section(section)
                retval = None
            if retval in [None,'']:# it wasn't found from the config file
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
def get_notebook_dir(*args):
    r'''Returns the notebook directory.  If arguments are passed, it returns the directory underneath the notebook directory, ending in a trailing (back)slash
    
    It is determined by a call to `MyConfig.get_setting` with the environment variable set to ``PYTHON_NOTEBOOK_DIR`` and default ``~/notebook``.
    '''
    base_notebook_dir = _my_config.get_setting('notebook_directory',environ = 'PYTHON_NOTEBOOK_DIR',default = '~/notebook')
    retval = (base_notebook_dir,) + args
    if len(retval[-1]) != 0:
        retval = retval + ('',)
    return os.path.join(*retval)
def dirformat(file):
    raise ValueError("dirformat should now be obsolete, and replaced by the ability to pass a list of subdirectories to getDATADIR or get_notebook_dir")
def grab_data_directory():
    raise RuntimeError("This used to return two arguments,"
            +"\n\t(1) if it successfully grabbed the data directory from"
            +"a file and set the PYTHONDATADIR environment variable"
            +"\n\t(2) set to true if it couldn't figure out where the"
            +"data directory was"
            +"\nAll the functionality of this function should now be"
            +"replaced by getDATADIR")
def getDATADIR(*args,**kwargs):
    r'''Returns the base directory where you put all your data.  If arguments
    are passed, it returns the directory underneath the data directory, ending
    in a trailing (back)slash
    
    It is determined by a call to `MyConfig.get_setting` with the setting name
    `experimental_data` and the environment variable set to ``PYTHON_DATA_DIR``
    and default ``~/experimental_data``.

    Parameters
    ----------
    exp_type : str
        A string identifying the experiment type.
        It searches for `exp_type` in this order:

        * Look in the ``ExpTypes`` section of the config file.
            * Note that by using this, you can store data in locations other
                than your `experimental_data` directory.
                For example, consider the following section of the
                ``~/.pyspecdata`` config file:
                ```
                [ExpTypes]
                alternate_base = /opt/other_data
                alternate_type_one = %(alternate_base)s/type_one
                ```
                which would find data with `exp_type` ``alternate_type_one`` in
                ``/opt/other_data/type_one``.
        * use `os.walk` to search for a directory with this name,
            starting from the data identified by `experimental_data`.
            excluding things that start with '.', '_' or
            containing '.hfssresults', always choosing the
            thing that's highest up in the tree.
    '''
    exp_type, = process_kwargs([('exp_type',None)],kwargs)
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
            # {{{ I would like to do something like the following, but it's not allowed in either ConfigParser or SafeConfigParser
            #base_dir = _my_config.get_setting('ExpTypes','base')
            #if base_dir is None: _my_config.set_setting('ExpTypes','base',base_dir)
            #if base_dir in exp_directory: exp_directory = [exp_directory.replace(base_dir,'%(base)s')]
            # }}}
            _my_config.set_setting('ExpTypes',exp_type,exp_directory)
        retval = (exp_directory,) + args
        # }}}
    else:
        retval = (base_data_dir,) + args
    if len(retval[-1]) != 0:
        retval = retval + ('',)
    return os.path.join(*retval)
