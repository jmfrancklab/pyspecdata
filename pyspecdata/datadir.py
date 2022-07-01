r'''Allows the user to run the same code on different machines,
even though the location of the raw spectral data might change.

This is controlled by the ``~/.pyspecdata`` or ``~/_pyspecdata`` config file.
'''
import os, sys
import configparser
import platform
from .general_functions import process_kwargs, strm
import logging
import atexit
from subprocess import Popen,PIPE,check_output
from subprocess import call as subprocess_call
logger = logging.getLogger('pyspecdata.datadir')
class MyConfig(object):
    r'''Provides an easy interface to the pyspecdata configuration file.
    Only one instance pyspec_config should be created -- this instance is used by
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
        atexit.register(self.__exit__,None,None,None)
        return
    def set_setting(self,this_section,this_key,this_value):
        "set `this_key` to `this_value` inside section `this_section`, creating it if necessary"
        if self._config_parser is None:
            self._config_parser = configparser.SafeConfigParser()
            read_cfg = self._config_parser.read(self.config_location)
            if not read_cfg:
                print("\nWarning!! There was no file at",self.config_location,"so I'm creating one")
        if not self._config_parser.has_section(this_section):
            self._config_parser.add_section(this_section)
        self._config_parser.set(this_section,this_key,this_value)
        return
    def __exit__(self, exception_type, exception_value, traceback):
        if self._config_parser is not None:
            # {{{ reset to standard figures
            self.set_setting('mode','figures','standard')
            # }}}
            with open(self.config_location,'w') as fp:
                self._config_parser.write(fp)
        return
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
            logger.debug(strm("I pulled",this_key,"from the config_vars"))
            return self.config_vars[this_key]
        if environ is not None and environ in os.environ.keys():
            logger.debug(strm("about to look for environment variable",environ))
            retval = os.environ[environ]
            logger.debug(strm("I pulled",environ,"from the environment variables -- it is",retval))
        else:
            if self._config_parser is None:
                self._config_parser = configparser.SafeConfigParser()
                read_cfg = self._config_parser.read(self.config_location)
                if not read_cfg:
                    logger.debug("\nWarning!! There was no file at",self.config_location,"so I'm creating one")
            if self._config_parser.has_section(section):
                try:
                    retval = self._config_parser.get(section,this_key)
                except configparser.NoOptionError:
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
                retval = default
                if retval is None:
                    raise RuntimeError("\nI didn't find the value corresponding to "+this_key
                            +" in the environment variable "+repr(environ)+'\n'+
                            "--> You probably want to run the command-line tool pyspecdata_dataconfig to set up a configuration file")
                self._config_parser.set(section,this_key,retval)
            if environ is not None:
                os.environ[environ] = retval
            logger.debug(strm("I pulled",this_key,"from the configuration file -- it is",retval))
        self.config_vars[this_key] = retval
        return retval
pyspec_config = MyConfig()
def get_notebook_dir(*args):
    r'''Returns the notebook directory.  If arguments are passed, it returns the directory underneath the notebook directory, ending in a trailing (back)slash
    
    It is determined by a call to `MyConfig.get_setting` with the environment variable set to ``PYTHON_NOTEBOOK_DIR`` and default ``~/notebook``.
    '''
    base_notebook_dir = pyspec_config.get_setting('notebook_directory',environ = 'PYTHON_NOTEBOOK_DIR',default = '~/notebook')
    if not os.path.exists(base_notebook_dir):
        base_notebook_dir = os.path.expanduser(base_notebook_dir)
        if not os.path.exists(base_notebook_dir):
            raise ValueError("It seems that your notebook directory (the main directory containing your latex files) isn't either (1) called \"notebook\" and immediately underneath your home directory or (2) registered in the [General] block of your "+pyspec_config.config_location+"file.\nThis probably means that you want to add a line like\nnotebook_directory = [path to your main notebook directory here]\nTo the [General] block of "+pyspec_config.config_location)
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
    r'''Used to find a directory containing data in a way that works
    seamlessly across different computers (and operating systems).
    Supports the case where data is processed both on a laboratory
    computer and (*e.g.* after transferring via ssh or a syncing client) on a
    user's laptop.
    While it will return a default directory without any arguments, it is
    typically used with the keyword argument `exp_type`, described below.

    Note that **the most common way** to use this mechanism is to set up your directories using
    the pyspecdata_register_dir shell command -- see
    :func:`~pyspecdata.datadir.register_directory`.
    

    It returns the directory ending in a trailing (back)slash.
    
    It is determined by a call to `MyConfig.get_setting` with the setting name
    `data_directory` and the environment variable set to ``PYTHON_DATA_DIR``.

    Parameters
    ----------
    exp_type : str
        A string identifying the name of a subdirectory where the data is stored.
        It can contain slashes.
        Typically, this gives the path relative to a google drive, rclone,
        dropbox, etc, repository.
        To make code portable, `exp_type` should **not** contain a full path or
        or portions of the path that are specific to the computer/user.

        If the directory has note been used before, all the directories listed
        in the user's `_pyspecdata` or `.pyspecdata` config file will be
        searched recursively up to 2 levels deep.

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
    exp_type = process_kwargs([('exp_type',None)],kwargs)
    base_data_dir = pyspec_config.get_setting('data_directory',environ = 'PYTHON_DATA_DIR',default = '~/experimental_data')
    if exp_type is not None and '/' in exp_type:
        exp_type = os.path.normpath(exp_type)
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
        logger.debug(strm("Walking inside",walking_top_dir,"to find",exp_type,"will only walk 2 directories deep!"))
        equal_matches = []
        containing_matches = []
        for d,s,_ in walklevel(walking_top_dir, level=2):
            logger.debug(strm("walking: ",d,s))
            s[:] = [j for j in s if j[0]
                    not in ['.','_']
                    and '.hfssresults' not in j]# prune the walk
            equal_matches.extend([os.path.join(d,j) for j in s if
                os.path.join(d,j).lower().endswith(exp_type.lower())])
            containing_matches.extend([os.path.join(d,j) for j in s if
                exp_type.lower() in os.path.join(d,j).lower()])
        def grab_smallest(matches):
            if len(matches) > 0:
                if len(matches) == 1:
                    return matches[0]
                else:
                    min_length_match = min(list(map(len,matches)))
                    matches = [x for x in matches if len(x) == min_length_match]
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
        #base_dir = pyspec_config.get_setting('ExpTypes','base')
        #if base_dir is None: pyspec_config.set_setting('ExpTypes','base',base_dir)
        #if base_dir in exp_directory: exp_directory = [exp_directory.replace(base_dir,'%(base)s')]
        # }}}
        pyspec_config.set_setting('ExpTypes',exp_type,exp_directory)
        return exp_directory
    if exp_type is not None:
        # {{{ determine the experiment subdirectory
        exp_directory = pyspec_config.get_setting(exp_type, section='ExpTypes')
        if exp_directory is None:
            logger.debug(strm("I found no directory matches for exp_type "+exp_type+", so now I want to look inside all the known exptypes"))
            for t,d in dict(pyspec_config._config_parser.items('ExpTypes')).items():
                exp_directory = walk_and_grab_best_match(d)
                if exp_directory is not None:
                    break
            if exp_directory is None:
                logger.debug(strm("I found no directory matches for exp_type "+exp_type+", even after searching inside all the known exptypes"))
                d = dict(pyspec_config._config_parser.items('General'))['data_directory']
                exp_directory = walk_and_grab_best_match(d)
                if exp_directory is None:
                    raise ValueError("even after walking the whole data directory, I can't find a match for "+exp_type+".  If that's really a directory that exists on a remote server, etc, then you should add the empty directory to your local file structure, somewhere where it's findable by pyspecdata (listed in your pyspecdata config file)")
        if exp_directory is None:
            logger.debug(strm("I found no directory matches for exp_type "+exp_type+", after walking the known exptypes, so I'm going to walk data_directory"))
            exp_directory = walk_and_grab_best_match(base_data_dir)
        retval = (exp_directory,) + args
        # }}}
    else:
        retval = (base_data_dir,) + args
    if len(retval[-1]) != 0:
        retval = retval + ('',)
    return os.path.join(*retval)
def rclone_search(fname,exp_type,dirname):
    logger.debug(strm("rclone search called with",fname,exp_type))
    remotelocation = pyspec_config.get_setting('RcloneRemotes',exp_type)
    if remotelocation is not None:
        logger.debug("remote location previously stored")
        cmd = strm(
                "rclone copy -v --include '%s' %s %s"%(fname,
                    remotelocation,
                    # dirname below needs to be replaced with path relative to current directory
                    os.path.normpath(os.path.join(dirname)).replace('\\','\\\\')))
        retval += '\nBased on previous searches for this exp_type, you should be able to retrieve this file with:\n'+cmd
        return retval
    retval = "I can't find %s in %s, so I'm going to search for t in your rclone remotes"%(fname,exp_type)
    rclone_remotes = []
    try:
        with Popen('rclone', stdout=PIPE, stderr=PIPE) as p:
            pass
    except OSError as e:
        raise RuntimeError("I can't find the rclone program in your"
                "path -- please install from http://rclone.org")
    with Popen(['rclone','listremotes'],
            stdout=PIPE, stderr=PIPE, encoding='utf-8') as proc:
        for j in proc.stdout:
            rclone_remotes.append(j.rstrip())
    for thisremote in rclone_remotes:
        retval += '\n'+strm("checking remote",thisremote)
        # do NOT quote the filename -- quotes are typically stripped off by the
        # shell -- they would be literal here
        cmd = ['rclone','--include',f'**{exp_type}**/*{fname}*', 'ls',thisremote]
        logger.info("trying to find a file, and running:\n\t"+strm(*cmd))
        with Popen(cmd, stdout=PIPE, stderr=PIPE, encoding='utf-8') as proc:
            for j in proc.stdout:
                foundpath = j.split()
                foundpath = [' '.join(foundpath[1:])] # leave out size, just filename
                lastlen = 0
                while len(foundpath) > lastlen:
                    lastlen = len(foundpath)
                    foundpath = [j for k in foundpath
                            for j in os.path.split(k)
                            if len(j) > 0]
                    logging.debug(strm("foundpath is",foundpath))
                remotelocation = thisremote+'/'.join(foundpath[:-1])
                logging.debug("about to write to RcloneRemotes")
                pyspec_config.set_setting('RcloneRemotes',exp_type,remotelocation)
                cmd = strm(
                        "rclone copy -v --include '%s' %s %s"%(foundpath[-1],
                            remotelocation,
                            # dirname below needs to be replaced with path relative to current directory
                            os.path.normpath(os.path.join(dirname)).replace('\\','\\\\')))
                retval += '\nYou should be able to retrieve this file with:\n'+cmd
        logger.debug("Popen has finally terminated!")
    return retval
def log_fname(logname,fname,dirname,exp_type,err=False):
    r"""logs the file name either used or missing.

    Also, by setting the `err` flag to True, you can generate an error message
    that will guide you on how to selectively copy down this data from a remote source
    (google drive, etc.), *e.g.*:

    ``Traceback (most recent call last):
      File "proc_square_refl.py", line 21, in <module>
        directory=getDATADIR(exp_type='test_equip'))
      File "c:\users\johnf\notebook\pyspecdata\pyspecdata\core.py", line 6630, in __init__
        check_only=True, directory=directory)
      File "c:\users\johnf\notebook\pyspecdata\pyspecdata\core.py", line 1041, in h5nodebypath
        +errmsg)
    AttributeError: You're checking for a node in a file (200110_pulse_2.h5) that does not exist
    I can't find 200110_pulse_2.h5 in C:\Users\johnf\exp_data\test_equip\, so I'm going to search for t in your rclone remotes
    checking remote g_syr:
    You should be able to retrieve this file with:
    rclone copy -v --include '200110_pulse_2.h5' g_syr:exp_data/test_equip C:\\Users\\johnf\\exp_data\\test_equip``
    """
    if err:
        logger.debug(strm("about to call rclone search on fname:",fname,"dirname:",exp_type))
        rclone_suggest = rclone_search(fname,exp_type,dirname)# eventually lump into the error message
        logger.debug("rclone search done")
    with open(logname+'.log','a+',encoding='utf-8') as fp:
        already_listed = False
        fp.seek(0,0)
        for j in fp:
            j = j.replace(r'\ ','LITERALSPACE')
            try:
                f, e, d = j.split()
            except:
                raise RuntimeError(strm("there seems to be something wrong with your",logname+'.log',"file (in the current directory).  It should consist of one line per file, with each file containing a file and directory name.  Instead, I find a line with the following elements",j.split(),'\n',"You might try deleting the",logname+'.log',"file"))
            f = f.replace('LITERALSPACE',' ')
            d = d.replace('LITERALSPACE',' ')
            if f == fname and d == dirname:
                already_listed = True
                break
        if not already_listed:
            fp.seek(0,os.SEEK_END)# make sure at end of file
            fp.write('%-70s%-50s%-50s\n'%(fname.replace(' ','\\ '),exp_type.replace(' ','\\ '),dirname.replace(' ','\\ ')))
    if err:
        return rclone_suggest
def register_directory():
    r"""The shell command `pyspecdata_register_dir WHICHDIR` will register the
    directory WHICHDIR (substitute with the name of a directory on your
    computer) so that it can be automatically discovered by
    :func:`~pyspecdata.find_file` or
    :func:`~pyspecdata.search_filename`
    after executing this shell command
    you can use the `exp_type` argument of those commands where you only give
    the lowest level subdirectory (or the final couple subdirectories) that
    contains your data.

    Key to the way this mechanism works

    .. note::
        this feature was installed on 9/24/20: you need to re-run
        `setup.py` in order to get this command to work for the first time if
        you installed pyspecdata before that date.
    """
    assert len(sys.argv) == 2,"Only give one argument -- the directory!"
    exp_directory = sys.argv[1]
    exp_directory = os.path.normpath(os.path.expanduser(exp_directory))
    _,exp_type = os.path.split(exp_directory)
    logger.debug(strm("trying to register directory",exp_directory,"as",exp_type))
    pyspec_config.set_setting('ExpTypes',exp_type,exp_directory)
    return
