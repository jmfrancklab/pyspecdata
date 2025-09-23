r"""Allows the user to run the same code on different machines,
even though the location of the raw spectral data might change.

This is controlled by the ``~/.pyspecdata`` or ``~/_pyspecdata`` config file.
"""

import os, sys, csv
import configparser
import platform
from difflib import get_close_matches
from .general_functions import process_kwargs, strm
import logging
import atexit
import collections
from subprocess import Popen, PIPE
from pathlib import Path, PureWindowsPath

logger = logging.getLogger("pyspecdata.datadir")
unknown_exp_type_name = "XXX--unknown--XXX"


class MyConfig(object):
    r"""Provides an easy interface to the pyspecdata configuration file.
    Only one instance pyspec_config should be created -- this instance is used
    by the other functions in this module.
    """

    def __init__(self):
        self._config_parser = None
        if platform.platform().startswith("Windows"):
            self.hide_start = "_"  # the default hidden/config starter for
            #                        vim, mingw applications, etc
            "This filename prefix denotes a configuration file on the OS."
        else:
            self.hide_start = "."
        self.config_location = os.path.join(
            os.path.expanduser("~"), self.hide_start + "pyspecdata"
        )
        self.config_vars = collections.defaultdict(dict)
        (
            "The dictionary that stores the current settings -- keep these in"
            " a dictionary, which should be faster than reading from environ,"
            " or from a file."
        )
        atexit.register(self.__exit__, None, None, None)
        return

    def set_setting(self, this_section, this_key, this_value):
        """set `this_key` to `this_value` inside section `this_section`,
        creating it if necessary"""
        if any(
            list(
                type(j) is not str
                for j in [this_section, this_key, this_value]
            )
        ):
            raise ValueError(
                "one of section, key, value are not a string! They are: "
                + str(
                    list(type(j) for j in [this_section, this_key, this_value])
                )
            )
        if self._config_parser is None:
            self._config_parser = configparser.ConfigParser()
            read_cfg = self._config_parser.read(self.config_location)
            if not read_cfg:
                print(
                    "\nWarning!! There was no file at",
                    self.config_location,
                    "so I'm creating one",
                )
        if not self._config_parser.has_section(this_section):
            self._config_parser.add_section(this_section)
        self._config_parser.set(this_section, this_key, this_value)
        return

    def __exit__(self, exception_type, exception_value, traceback):
        if self._config_parser is not None:
            # {{{ reset to standard figures
            self.set_setting("mode", "figures", "standard")
            # }}}
            with open(self.config_location, "w") as fp:
                self._config_parser.write(fp)
        return

    def get_setting(
        self, this_key, environ=None, default=None, section="General"
    ):
        """Get a settings from the "General" group.
        If the file does not exist, or the option is not set, then set the
        option, creating the file as needed.
        The option is determined in the following order:

            * The value in the `config_vars` dictionary.
            * The value of the environment variable named `environ`.
            * The value stored in the configuration file at ``~/.pyspecdata``
              (``~/_pyspecdata`` on Windows).
            * The value given by the `default` argument.  If `default` is
              ``None``, then return ``None``.

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

            For platform compatibility, leading period characters are
            converted to `self.hide_start`, and '/' are converted to
            `os.path.sep`.
        section : str
            Look in this section of the config file.

        Returns
        -------
        The value corresponding to `this_key`.
        """
        if this_key in self.config_vars[section].keys():
            logger.debug(
                strm(
                    "I pulled",
                    this_key,
                    "in",
                    section,
                    "from the",
                    "config_vars",
                )
            )
            return self.config_vars[section][this_key]
        if environ is not None and environ in os.environ.keys():
            logger.debug(
                strm("about to look for environment variable", environ)
            )
            retval = os.environ[environ]
            logger.debug(
                strm(
                    "I pulled",
                    environ,
                    "from the environment variables -- it is",
                    retval,
                )
            )
        else:
            if self._config_parser is None:
                self._config_parser = configparser.ConfigParser()
                read_cfg = self._config_parser.read(self.config_location)
                if not read_cfg:
                    logger.debug(
                        "\nWarning!! There was no file at",
                        self.config_location,
                        "so I'm creating one",
                    )
            if self._config_parser.has_section(section):
                try:
                    retval = self._config_parser.get(section, this_key)
                except configparser.NoOptionError:
                    retval = None
            else:
                self._config_parser.add_section(section)
                retval = None
            if retval in [None, ""]:  # it wasn't found from the config file
                if default is None:
                    return None
                if type(default) is str and "/" in default:
                    default = default.split("/")
                    for j in range(len(default)):
                        if default[j][0] == ".":
                            default[j] = self.hide_start + default[j][1:]
                    default = os.path.sep.join(default)
                retval = default
                if retval is None:
                    raise RuntimeError(
                        "\nI didn't find the value corresponding to "
                        + this_key
                        + " in the environment variable "
                        + repr(environ)
                        + "\n"
                        + "--> You probably want to run the command-line tool"
                        + " pyspecdata_dataconfig to set up a configuration"
                        + " file"
                    )
                if type(retval) in [int, float]:
                    retval = str(retval)
                self._config_parser.set(section, this_key, retval)
            if environ is not None:
                os.environ[environ] = retval
            logger.debug(
                strm(
                    "I pulled",
                    this_key,
                    "from the configuration file -- it is",
                    retval,
                )
            )
        self.config_vars[section][this_key] = retval
        return retval


pyspec_config = MyConfig()


def get_notebook_dir(*args):
    r"""Returns the notebook directory.  If arguments are passed, it returns
    the directory underneath the notebook directory, ending in a trailing
    (back)slash

    It is determined by a call to `MyConfig.get_setting` with the environment
    variable set to ``PYTHON_NOTEBOOK_DIR`` and default ``~/notebook``.
    """
    base_notebook_dir = pyspec_config.get_setting(
        "notebook_directory",
        environ="PYTHON_NOTEBOOK_DIR",
        default="~/notebook",
    )
    if not os.path.exists(base_notebook_dir):
        base_notebook_dir = os.path.expanduser(base_notebook_dir)
        if not os.path.exists(base_notebook_dir):
            raise ValueError(
                "It seems that your notebook directory (the main"
                "directory containing your latex files) isn't either (1)"
                'called "notebook" and immediately underneath your home'
                "directory or (2) registered in the [General] block of"
                " your"
                + pyspec_config.config_location
                + "file.\nThis probably means that you want to add a line"
                " like \nnotebook_directory = [path to your main notebook"
                " directory here]\nTo the [General] block of "
                + pyspec_config.config_location
            )
    retval = (base_notebook_dir,) + args
    if len(retval[-1]) != 0:
        retval = retval + ("",)
    return os.path.join(*retval)


def dirformat(file):
    raise ValueError(
        "dirformat should now be obsolete, and replaced by the ability to"
        " pass a list of subdirectories to getDATADIR or get_notebook_dir"
    )


def grab_data_directory():
    raise RuntimeError(
        "This used to return two arguments,"
        + "\n\t(1) if it successfully grabbed the data directory from"
        + "a file and set the PYTHONDATADIR environment variable"
        + "\n\t(2) set to true if it couldn't figure out where the"
        + "data directory was"
        + "\nAll the functionality of this function should now be"
        + "replaced by getDATADIR"
    )


def proc_data_target_dir(exp_type):
    """A convenience function for getting a data directory you're going to
    write data to.

    Provides the full path to the directory corresponding to exp_type.

    If this is the first time you're trying to write to that exp_type, it will
    create the directory
    """
    data_target = os.path.normpath(getDATADIR("AG_processed_data"))
    if not os.path.exists(data_target):
        os.mkdir(data_target)
    return data_target


def getDATADIR(*args, **kwargs):
    r"""Used to find a directory containing data in a way that works
    seamlessly across different computers (and operating systems).

    **This is not intended as a user-level function** use
    :func:`~pyspecdata.find_file`
    or
    :func:`~pyspecdata.search_filename` (especially with the `unique`
    parameter set to true) instead!

    Supports the case where data is processed both on a laboratory
    computer and (*e.g.* after transferring via ssh or a syncing client) on a
    user's laptop.
    While it will return a default directory without any arguments, it is
    typically used with the keyword argument `exp_type`, described below.

    Note that **the most common way** to use this mechanism is to set up your
    directories using
    the pyspecdata_register_dir shell command -- see
    :func:`~pyspecdata.datadir.register_directory`.


    It returns the directory ending in a trailing (back)slash.

    It is determined by a call to `MyConfig.get_setting` with the setting name
    `data_directory` and the environment variable set to ``PYTHON_DATA_DIR``.

    Parameters
    ----------
    exp_type : str
        A string identifying the name of a subdirectory where the data is
        stored.
        It can contain slashes.
        Typically, this gives the path relative to a google drive, rclone,
        dropbox, etc, repository.
        To make code portable, `exp_type` should **not** contain a full path
        or or portions of the path that are specific to the computer/user.

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
                which would find data with `exp_type` ``alternate_type_one``
                in ``/opt/other_data/type_one``.
        * use `os.walk` to search for a directory with this name
            inside the directory identified by `experimental_data`.
            excluding things that start with '.', '_' or
            containing '.hfssresults', always choosing the
            thing that's highest up in the tree.
            If it doesn't find a directory inside `experimental_data`, it will
            search inside all the directories already listed in `ExpTypes`.
            Currently, in both attempts, it will only walk 2 levels deep
            (since NMR directories can be rather complex, and otherwise it
            would take forever).
    """
    exp_type = process_kwargs([("exp_type", None)], kwargs)
    if exp_type is not None:
        exp_type_path = PureWindowsPath(exp_type)
        exp_type_key = exp_type_path.as_posix()
        exp_type_match_value = exp_type_path.as_posix().casefold()
    else:
        exp_type_path = None
        exp_type_key = None
        exp_type_match_value = None
    base_data_dir = pyspec_config.get_setting(
        "data_directory",
        environ="PYTHON_DATA_DIR",
        default="~/experimental_data",
    )
    if base_data_dir is not None:
        base_data_dir = Path(base_data_dir).expanduser()

    # the following is from
    # https://stackoverflow.com/questions/229186
    # /os-walk-without-digging-into-directories-below
    def iter_candidate_dirs(some_dir, max_depth=1):
        some_dir = Path(some_dir)
        assert some_dir.is_dir(), strm(
            some_dir,
            "is not a directory (probably an invalid entry in your"
            "pyspecdata config file)",
        )
        queue = collections.deque([(some_dir, 0)])
        while queue:
            current_dir, depth = queue.popleft()
            if depth >= max_depth:
                continue
            for child in current_dir.iterdir():
                if not child.is_dir():
                    continue
                name = child.name
                if name[0] in [".", "_"] or ".hfssresults" in name:
                    continue
                yield child
                queue.append((child, depth + 1))

    def walk_and_grab_best_match(walking_top_dir):
        logger.debug(
            strm(
                "Walking inside",
                walking_top_dir,
                "to find",
                exp_type,
                "will only walk 2 directories deep!",
            )
        )
        equal_matches = []
        containing_matches = []
        pathlevel = len(exp_type_path.parts) or 1
        for candidate in iter_candidate_dirs(
            walking_top_dir, max_depth=pathlevel
        ):
            logger.debug(strm("walking: ", candidate))
            candidate_key = candidate.as_posix().casefold()
            if candidate_key.endswith(exp_type_match_value):
                equal_matches.append(candidate)
            if exp_type_match_value in candidate_key:
                containing_matches.append(candidate)

        def grab_smallest(matches):
            if len(matches) > 0:
                if len(matches) == 1:
                    return matches[0]
                else:
                    min_length_match = min(len(x.parts) for x in matches)
                    matches = [
                        x for x in matches if len(x.parts) == min_length_match
                    ]
                    if len(matches) != 1:
                        raise ValueError(
                            "I found multiple equivalent matches when"
                            + " searching for exp_type: "
                            + repr(matches)
                        )
                    return matches[0]

        if len(equal_matches) > 0:
            exp_directory = grab_smallest(equal_matches)
        elif len(containing_matches) > 0:
            exp_directory = grab_smallest(containing_matches)
        else:
            return None
        # {{{ I would like to do something like the following, but it's not
        # allowed in either ConfigParser or SafeConfigParser base_dir =
        # pyspec_config.get_setting('ExpTypes','base') if base_dir is None:
        # pyspec_config.set_setting('ExpTypes','base',base_dir) if base_dir in
        # exp_directory: exp_directory =
        # [exp_directory.replace(base_dir,'%(base)s')]
        # }}}
        pyspec_config.set_setting("ExpTypes", exp_type_key, str(exp_directory))
        return exp_directory

    if exp_type is not None:
        # {{{ determine the experiment subdirectory
        exp_directory = pyspec_config.get_setting(
            exp_type_key, section="ExpTypes"
        )
        if exp_directory is None and exp_type != exp_type_key:
            exp_directory = pyspec_config.get_setting(
                exp_type, section="ExpTypes"
            )
        if (
            exp_directory is None
            and pyspec_config._config_parser is not None
            and pyspec_config._config_parser.has_section("ExpTypes")
        ):
            for stored_key, stored_value in pyspec_config._config_parser.items(
                "ExpTypes"
            ):
                normalized_key = (
                    PureWindowsPath(stored_key).as_posix().casefold()
                )
                if normalized_key == exp_type_match_value:
                    exp_directory = stored_value
                    pyspec_config.set_setting(
                        "ExpTypes", exp_type_key, stored_value
                    )
                    break
        if exp_directory is not None:
            exp_directory = Path(exp_directory).expanduser()
        if exp_directory is None:
            logger.debug(
                strm(
                    "I found no directory matches for exp_type "
                    + exp_type
                    + ", so now I want to look inside all the known exptypes"
                )
            )
            for t, d in dict(
                pyspec_config._config_parser.items("ExpTypes")
            ).items():
                exp_directory = walk_and_grab_best_match(d)
                if exp_directory is not None:
                    break
            if exp_directory is None:
                logger.debug(
                    strm(
                        "I found no directory matches for exp_type "
                        + exp_type
                        + ", even after searching inside all the known"
                        + " exptypes"
                    )
                )
                d = dict(pyspec_config._config_parser.items("General"))[
                    "data_directory"
                ]
                d = Path(d).expanduser()

                exp_directory = walk_and_grab_best_match(d)
                if exp_directory is None:
                    partial_match_value, partial_match_prefix_len = (
                        find_registered_prefix_match(
                            exp_type_path,
                            section="ExpTypes",
                            prefer_longest=True,
                        )
                    )
                    if partial_match_value is not None:
                        partial_match_directory = Path(partial_match_value).expanduser()
                        unmatched_parts = exp_type_path.parts[
                            partial_match_prefix_len:
                        ]
                        target_dir = partial_match_directory.joinpath(
                            *unmatched_parts
                        )
                        mkdir_hint = " (from partial path match)"
                    else:
                        target_dir = d / Path(*exp_type_path.parts)
                        mkdir_hint = " (under datadir)"
                    message = (
                        "even after walking the whole data directory, I can't"
                        + " find a match for "
                        + exp_type
                        + ".  If that's really a directory that exists on a"
                        + " remote server, etc, then you should add the empty"
                        + " directory to your local file structure, somewhere"
                        + " where it's findable by pyspecdata (listed in your"
                        + " pyspecdata config file) -- i.e. mkdir -p "
                        + str(target_dir)
                        + mkdir_hint
                    )
                    print(message)
                    prompt = input(
                        "Would you like pyspecdata to create this directory"
                        " for you? [y/N]: "
                    ).strip()
                    if prompt.lower() in {"y", "yes"}:
                        target_dir.mkdir(parents=True, exist_ok=True)
                        pyspec_config.set_setting(
                            "ExpTypes", exp_type_key, str(target_dir)
                        )
                        exp_directory = target_dir
                    else:
                        raise ValueError(message)
        if exp_directory is None:
            logger.debug(
                strm(
                    "I found no directory matches for exp_type "
                    + exp_type
                    + ", after walking the known exptypes, so I'm going to"
                    + " walk data_directory"
                )
            )
            exp_directory = walk_and_grab_best_match(base_data_dir)
        exp_directory = Path(exp_directory)
        retval = (str(exp_directory),) + args
        # }}}
    else:
        retval = (str(base_data_dir),) + args
    if len(retval[-1]) != 0:
        retval = retval + ("",)
    return os.path.join(*retval)


def find_registered_prefix_match(
    exp_type_path,
    *,
    section,
    allow_full_match=True,
    prefer_longest=True,
):
    """Return the stored value corresponding to the best prefix match."""

    if exp_type_path is None:
        return None, 0
    if (
        pyspec_config._config_parser is None
        or not pyspec_config._config_parser.has_section(section)
    ):
        return None, 0
    exp_dir_parts = exp_type_path.parts
    if allow_full_match:
        max_prefix_len = len(exp_dir_parts)
    else:
        max_prefix_len = len(exp_dir_parts) - 1
    if max_prefix_len <= 0:
        return None, 0
    matches = []
    for prefix_len in range(max_prefix_len, 0, -1):
        partial_key = PureWindowsPath(*exp_dir_parts[:prefix_len]).as_posix()
        partial_key_casefold = partial_key.casefold()
        for stored_key, stored_value in pyspec_config._config_parser.items(section):
            normalized_key = PureWindowsPath(stored_key).as_posix().casefold()
            if normalized_key == partial_key_casefold:
                if prefer_longest:
                    return stored_value, prefix_len
                matches.append((prefix_len, stored_value))
                break
    if not prefer_longest and matches:
        return matches[-1][1], matches[-1][0]
    return None, 0


class cached_searcher(object):
    def __init__(self):
        self.has_run = False
        self.dirlist = []

    def grab_dirlist(self, specific_remote=None):
        if specific_remote is None:
            rclone_remotes = []
            try:
                with Popen("rclone", stdout=PIPE, stderr=PIPE) as _:
                    pass
            except OSError:
                raise RuntimeError(
                    "I can't find the rclone program in your"
                    "path -- please install from http://rclone.org"
                )
            with Popen(
                ["rclone", "listremotes"],
                stdout=PIPE,
                stderr=PIPE,
                encoding="utf-8",
            ) as proc:
                for j in proc.stdout:
                    rclone_remotes.append(j.rstrip())
        else:
            rclone_remotes = [specific_remote]
        for thisremote in rclone_remotes:
            print(
                "checking remote", thisremote, "... this might take a minute"
            )
            cmd = ["rclone", "lsf", "-R", "--dirs-only", thisremote]
            logger.info(
                "grabbing all dir info with "
                + strm(*cmd)
                + " ... this might take a minute"
            )
            with Popen(
                cmd, stdout=PIPE, stderr=PIPE, encoding="utf-8"
            ) as proc:
                for j in proc.stdout:
                    child = PureWindowsPath(j.strip())
                    remote_root = PureWindowsPath(thisremote)
                    child_fragment = child.as_posix()
                    if child_fragment == ".":
                        combined_str = remote_root.as_posix()
                    elif remote_root.as_posix().endswith(":"):
                        combined_str = remote_root.as_posix() + child_fragment
                    else:
                        combined_str = (remote_root / child).as_posix()
                    combined = PureWindowsPath(combined_str)
                    self.dirlist.append(combined)
            logger.debug(f"done checking that remote -- got {self.dirlist}")
        self.has_run = True
        return

    def search_for(
        self,
        exp_type,
        specific_remote=None,
        *,
        require_match=False,
        suggest_limit=3,
    ):
        exp_type_path = PureWindowsPath(exp_type)
        exp_type_casefold = exp_type_path.as_posix().casefold()

        def _raise_no_match(search_space):
            if not require_match:
                return []
            if not search_space:
                search_space = self.dirlist
            search_space_strings = [
                (
                    candidate.as_posix()
                    if hasattr(candidate, "as_posix")
                    else str(candidate)
                )
                for candidate in search_space
            ]
            lower_to_original = {
                candidate.casefold(): candidate
                for candidate in search_space_strings
            }
            close_matches = get_close_matches(
                exp_type_casefold,
                list(lower_to_original.keys()),
                n=suggest_limit,
            )
            if close_matches:
                suggestions = [
                    lower_to_original[match] for match in close_matches
                ]
            else:
                suggestions = search_space_strings[:suggest_limit]
            suggestion_text = (
                " Did you mean one of these?"
                + "\n\t•\t"
                + "\n\t•\t".join(suggestions)
                if suggestions
                else ""
            )
            location_hint = (
                f" within {specific_remote}"
                if specific_remote is not None
                else ""
            )
            raise ValueError(
                "I can't find a remote directory matching "
                f"{exp_type_path.as_posix()}{location_hint}."
                + suggestion_text
            )

        logger.debug(
            f"search_for exp_type={exp_type_path} inside"
            f" specific_remote={specific_remote}"
        )
        if not self.has_run:
            self.grab_dirlist(specific_remote=specific_remote)
        if specific_remote is not None:
            remote_path = PureWindowsPath(specific_remote)
            logger.debug(
                "drilling down, I'm looking for for"
                f" {(remote_path / exp_type_path).as_posix()} inside"
                f" {self.dirlist}"
            )
            target = (remote_path / exp_type_path).as_posix().casefold()
            remote_casefold = remote_path.as_posix().casefold()
            potential_hits = []
            remote_space = []
            for candidate in self.dirlist:
                candidate_posix = candidate.as_posix()
                candidate_key = candidate_posix.casefold()
                if candidate_key.startswith(remote_casefold):
                    remote_space.append(candidate)
                if candidate_key == target:
                    logger.debug(
                        "found exactly the right exp_type inside the"
                        " specific_remote"
                    )
                    potential_hits.append(candidate_posix)
            if not potential_hits:
                return _raise_no_match(remote_space)
            return potential_hits
        potential_hits = []
        for candidate in self.dirlist:
            candidate_posix = candidate.as_posix()
            if exp_type_casefold in candidate_posix.casefold():
                potential_hits.append(candidate_posix)
        logger.debug(
            strm(
                "found potential hits",
                potential_hits,
                "when looking for",
                exp_type_casefold,
            ),
            "inside",
            self.dirlist,
        )
        if not potential_hits:
            return _raise_no_match(self.dirlist)
        return potential_hits


cached_searcher_instance = cached_searcher()


def rclone_search(fname, exp_type, dirname):
    logger.debug(
        strm("rclone search called with fname:", fname, "exp_type:", exp_type)
    )
    exp_type_path = PureWindowsPath(exp_type)
    exp_type_casefold = exp_type_path.as_posix().casefold()
    remotelocation = pyspec_config.get_setting(
        exp_type_casefold, section="RcloneRemotes"
    )
    if (
        remotelocation is None
        and pyspec_config._config_parser is not None
        and pyspec_config._config_parser.has_section("RcloneRemotes")
    ):
        for stored_key, stored_value in pyspec_config._config_parser.items(
            "RcloneRemotes"
        ):
            normalized_key = PureWindowsPath(stored_key).as_posix().casefold()
            if normalized_key == exp_type_casefold:
                remotelocation = stored_value
                pyspec_config.set_setting(
                    "RcloneRemotes", exp_type_casefold, stored_value
                )
                break
    result = []
    if remotelocation is None:
        # {{{ first see the exp_type is contained inside something else
        partial_remote, partial_prefix_len = find_registered_prefix_match(
            exp_type_path,
            section="RcloneRemotes",
            allow_full_match=False,
            prefer_longest=False,
        )
        if partial_remote is not None:
            remotelocation = partial_remote
            logger.debug(f"did find one level up {remotelocation}")
            unmatched_parts = exp_type_path.parts[partial_prefix_len:]
            if unmatched_parts:
                result = cached_searcher_instance.search_for(
                    PureWindowsPath(*unmatched_parts),
                    specific_remote=remotelocation,
                    require_match=True,
                )
            else:
                result = [remotelocation]
        # }}}
        # {{{ only do a general search if the previous failed
        if remotelocation is None:
            logger.debug(
                f"remote location {exp_type.lower()} not previously stored,"
                " so search for it!"
            )
            result = cached_searcher_instance.search_for(
                exp_type_path, require_match=True
            )
        # }}}
        if len(result) > 1:
            raise ValueError(
                f"the exp_type that you've selected, {exp_type}, is ambiguous"
                ", and could refer to the following remote locations:\n"
                + str(result)
            )
        else:
            remotelocation = result[0]
            logging.debug("about to write to RcloneRemotes")
            pyspec_config.set_setting(
                "RcloneRemotes", exp_type_casefold, str(remotelocation)
            )
    logger.debug(
        f"remote location previously stored {remotelocation} for"
        f"{exp_type} going to be put in {dirname}"
    )
    if not fname.startswith("*"):
        fname = "*" + fname
    if fname.endswith(r"\b"):
        fname = fname[:-2]
    elif not fname.endswith("*"):
        fname = fname + "*"
    cmd = strm(
        "rclone copy -v --include '%s' %s %s"
        % (
            fname,
            str(remotelocation),
            # dirname below needs to be replaced with path relative to current
            # directory
            os.path.normpath(os.path.join(dirname)).replace("\\", "\\\\"),
        )
    )
    cmd = cmd.replace("'", '"')
    logger.info(f"I'm about to run\n{cmd}")
    os.system(cmd)
    logger.info("... done")
    return cmd


def log_fname(logname, fname, dirname, exp_type):
    r"""logs the file name either used or missing into a csv file.

    Also, by setting the `err` flag to True, you can generate an error message
    that will guide you on how to selectively copy down this data from a
    remote source (google drive, etc.), *e.g.*:

    ``Traceback (most recent call last):
      File "proc_square_refl.py", line 21, in <module>
        directory=getDATADIR(exp_type='test_equip'))
      File "c:\users\johnf\notebook\pyspecdata\pyspecdata\core.py", line 6630,
        in __init__ check_only=True, directory=directory)
      File "c:\users\johnf\notebook\pyspecdata\pyspecdata\core.py", line 1041,
        in h5nodebypath +errmsg)
    AttributeError: You're checking for a node in a file (200110_pulse_2.h5)
      that does not exist
    I can't find 200110_pulse_2.h5 in C:\Users\johnf\exp_data\test_equip\, so
      I'm going to search for t in your rclone remotes
    checking remote g_syr:
    You should be able to retrieve this file with:
    rclone copy -v --include '200110_pulse_2.h5' g_syr:exp_data/test_equip
    C:\\Users\\johnf\\exp_data\\test_equip``
    """
    thefields = ["Filename", "Path", "exp_type"]
    therow = [fname, dirname, exp_type]
    # {{{ read all the info into a set
    if os.path.exists(logname + ".csv"):
        with open(logname + ".csv", "r", encoding="utf-8") as fp:
            reader = csv.DictReader(fp)
            all_data = []
            for row in reader:
                all_data.append(tuple(row[j] for j in thefields))
            all_data = set(all_data)
    else:
        all_data = set()
    # }}}
    # {{{ add our info to the set
    therow = dict(zip(thefields, therow))
    all_data |= set([tuple(therow[j] for j in thefields)])
    # }}}
    # {{{ write the updated csv
    with open(logname + ".csv", "w", encoding="utf-8") as fp:
        thiscsv = csv.DictWriter(fp, fieldnames=thefields)
        thiscsv.writeheader()
        for thisitem in all_data:
            thiscsv.writerow(dict(zip(thefields, thisitem)))
    # }}}
    return


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
    If the `exp_type` that you are trying to access has a slash in it, you
    should register the top-most directory.  (For example, if you want
    `UV_Vis/Liam`, then register the directory that provides `UV_Vis`).

    .. note::
        this feature was installed on 9/24/20: you need to re-run
        `setup.py` in order to get this command to work for the first time if
        you installed pyspecdata before that date.
    """
    assert len(sys.argv) == 2, "Only give one argument -- the directory!"
    exp_directory = Path(sys.argv[1]).expanduser().resolve(strict=False)
    exp_type = exp_directory.name
    exp_type_key = Path(PureWindowsPath(exp_type)).as_posix()
    logger.debug(
        strm("trying to register directory", exp_directory, "as", exp_type)
    )
    pyspec_config.set_setting("ExpTypes", exp_type_key, str(exp_directory))
    return
