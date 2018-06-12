r'''Allows the user to run the same code on different machines,
even though the location of the raw spectral data might change.

This is controlled by the ``~/.pyspecdata`` or ``~/_pyspecdata`` config file.
'''
import os
import ConfigParser
import platform
from .general_functions import process_kwargs, strm
import logging
logger = logging.getLogger('pyspecdata.datadir')
class MyConfig(object):
    r'''Provides an easy interface to the pyspecdata configuration file.
    Only one instance _my_config should be created -- this instance is used by
    the other functions in this module.
    '''
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
                raise RuntimeError("\nI didn't find the value corresponding to "+this_key
                        +" in the environment variable "+repr(environ)+'\n'+
                        "--> You probably want to run the command-line tool pyspecdata_dataconfig to set up a configuration file")
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
    if not os.path.exists(base_notebook_dir):
        base_notebook_dir = os.path.expanduser(base_notebook_dir)
        if not os.path.exists(base_notebook_dir):
            raise ValueError("It seems that your notebook directory (the main directory containing your latex files) isn't either (1) called \"notebook\" and immediately underneath your home directory or (2) registered in the [General] block of your "+_my_config.config_location+"file.\nThis probably means that you want to add a line like\nnotebook_directory = [path to your main notebook directory here]\nTo the [General] block of "+_my_config.config_location)
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
    `data_directory` and the environment variable set to ``PYTHON_DATA_DIR``.

    Parameters
    ----------
    exp_type : str
        A string identifying the experiment type.
        It searches for `exp_type` in this order:

        * Look in the ``ExpTypes`` section of the config file.
            * Note that by using this, you can store data in locations other
                than your main data directory.
                For example, consider the following section of the
                ``~/.pyspecdata`` config file:
                ```
                [ExpTypes]
                alternate_base = /opt/other_data
                alternate_type_one = %(alternate_base)s/type_one
                ```
                which would find data with `exp_type` ``alternate_type_one`` in
                ``/opt/other_data/type_one``.
        * use `os.walk` to search for a directory with this name
            inside the directory identified by `experimental_data`.
            excluding things that start with '.', '_' or
            containing '.hfssresults', always choosing the
            thing that's highest up in the tree.
            If it doesn't find a directory inside `experimental_data`, it will
            search inside all the directories already listed in `ExpTypes`.
            Currently, in both attempts, it will only walk 2 levels deep (since NMR directories
            can be rather complex, and otherwise it would take forever).
    '''
    # the following is from https://stackoverflow.com/questions/229186/os-walk-without-digging-into-directories-below
    def walklevel(some_dir, level=1):
        some_dir = some_dir.rstrip(os.path.sep)
        assert os.path.isdir(some_dir),strm(some_dir,"is not a directory (probably an invalid entry in your pyspecdata config file)")
        num_sep = some_dir.count(os.path.sep)
        for root, dirs, files in os.walk(some_dir):
            yield root, dirs, files
            num_sep_this = root.count(os.path.sep)
            if num_sep + level <= num_sep_this:
                del dirs[:]
    def walk_and_grab_best_match(walking_top_dir):
        logger.info(strm("Walking inside",walking_top_dir,"to find",exp_type,"will only walk 2 directories deep!"))
        equal_matches = []
        containing_matches = []
        for d,s,_ in walklevel(walking_top_dir, level=2):
            logger.debug(strm("walking: ",d,s))
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
            return None
        # {{{ I would like to do something like the following, but it's not allowed in either ConfigParser or SafeConfigParser
        #base_dir = _my_config.get_setting('ExpTypes','base')
        #if base_dir is None: _my_config.set_setting('ExpTypes','base',base_dir)
        #if base_dir in exp_directory: exp_directory = [exp_directory.replace(base_dir,'%(base)s')]
        # }}}
        _my_config.set_setting('ExpTypes',exp_type,exp_directory)
        return exp_directory
    exp_type = process_kwargs([('exp_type',None)],kwargs)
    base_data_dir = _my_config.get_setting('data_directory',environ = 'PYTHON_DATA_DIR',default = '~/experimental_data')
    if exp_type is not None:
        # {{{ determine the experiment subdirectory
        exp_directory = _my_config.get_setting(exp_type, section='ExpTypes')
        if exp_directory is None:
            exp_directory = walk_and_grab_best_match(base_data_dir)
            if exp_directory is None:
                logger.info(strm("I found no directory matches for exp_type "+exp_type+", so now I want to look inside all the known exptypes"))
                for t,d in dict(_my_config._config_parser.items('ExpTypes')).iteritems():
                    exp_directory = walk_and_grab_best_match(d)
                    if exp_directory is not None:
                        break
                if exp_directory is None:
                    raise ValueError(strm("I found no directory matches for exp_type "+exp_type+", even after searching inside all the known exptypes"))
        retval = (exp_directory,) + args
        # }}}
    else:
        retval = (base_data_dir,) + args
    if len(retval[-1]) != 0:
        retval = retval + ('',)
    return os.path.join(*retval)
