These instructions are for installing Anaconda with pySpecData,
both a Python 2 and Python 3 versions.

Please note that if you downloaded the “Python 2.X” version of anaconda,
that means that your package manager, *etc.*, use Python 2.X.
As of this writing (1/13/20), those versions of the package manager are poorly
maintained, and lead to a host of headaches.
So, *yes*, the solution is to *completely uninstall* Anaconda, and then reinstall.

Note that you are *strongly encouraged* to use the Python 3.X version of
pySpecData, and to report any errors that you find.
You can convert your old scripts to Python 3 using ``2to3 -w name_of_your_python_script.py``

The basic strategy here is to install the most recent version of
anaconda (3.7) so that we are using an up-to-date “conda” package
manager, but then install a virtual environment where we can run old
Python 2 code if we so choose.

These instructions should work either on a Windows computer or a Mac (using
Anaconda -- you can also do homebrew on a Mac, which is not covered here).

**If at any point during these instructions, you don't get the expected result, stop and seek help!**

## initial

remove current anaconda

reboot to complete uninstall

go [here](https://www.anaconda.com/distribution/) to download the Python
3.7 version of Anaconda, and install

## checking for applocker

This section only applies to Windows computers.

Open the anaconda prompt -- if your you computer has applocker (*i.e.*,
if it's a university computer), and you don't see (base) as part of the
command prompt (command prompt is just `C:\` , you need to do the
following steps:

run ``secpol`` from admin command line (right click “Command Prompt”, Run as ...), navigate to “application control”
and “executable rules”


create a new rule based on path -- manually enter
`C:\ProgramData\Anaconda3\`

execute all following commands inside anaconda prompt unless noted
otherwise

## changing permissions to allow package installation

Make the Anaconda3 folder writeable by all users:
*   On Windows, you achieve this by
    opening the
    `C:\ProgramData` folder in File Explorer right click on
    `Anaconda3`→ security tab→ edit→
    users→ click “full control” checkbox on bottom →
    apply (this takes a few minutes to run)
*   On Mac, locate where Anaconda was installed (default was /opt directory) and edit permissions via Finder window.
    You can accomplish this via the following commands.  cd /opt open . (open
    new Finder window at this location) right-click on the anaconda3 directory,
    select 'Get Info' from the list. In the pop-up window, find Sharing &
    Permissions section, and in the Name list, find everyone -- change
    Privilege from 'Read only' to 'Read & Write'.

## install pySpecData inside py3 environment

If using windows, do all the following **in the anaconda prompt**.
The anaconda
prompt should read (base) indicating that you are in the base (Python 3)
environment:

Make sure that in your git repo, you have checked out a python 3 branch (as of this
writing, master is python 3, and there is a py2 branch for python 2.7, which is planned for obsolescence)

Install various python running environments
`conda install -y -c anaconda jupyter ipython spyder`

Install pySpecData prerequisites from the documentation:
for Windows: `conda install -y -c anaconda numpy scipy sympy pyqt pytables matplotlib h5py libpython mingw`
For installation on Mac, do not include mingw
(libpython may or may not be necessary, depending on details.)
If you are installing as a developer, where you will want to rebuild the documentation, also run:
`conda install -y -c conda-forge make sphinx sphinx_rtd_theme`

On windows, make sure that `where gcc` (`which gcc` for bash) returns a file **with an .exe extension**; if it does not, you need to add the directory with gcc.exe (usually Anaconda\MinGW) to your windows path.
(Or to your bash path, if you've set up for bash below, and are trying this within bash).
*This is a problem with how anaconda sets up mingw.*

install pySpecData in the python 3 (base) environment
`python setup.py develop`

If, on Windows, you have an issue where it complains about a “NoneType” `ld_version`,
like [here](https://stackoverflow.com/questions/48764602/cygwincompiler-typeerror-not-supported-between-instances-of-nonetype-and),
then you need to find the location of gcc.exe and add it to your path.
You can do this easily and temporarily by running
either with something similar to
``set PATH=%PATH%;C:\ProgramData\Anaconda3\MinGW\bin``
in DOS/Power Shell, or with
``export PATH=$PATH;C:\ProgramData\Anaconda3\MinGW\bin``
If you are having trouble, it's useful to separately run `python setup.py
build` (which calls the compilers) and `python setup.py develop` (which will
then complete the process).

Alos, if you have previously built an old version of pySpecData and building
gives syntax errors, you may need to remove the pySpecData/build directory

## create a python 2 environment, and install basic tools 

If using windows, do the following **in the anaconda prompt**:

`conda create –-name py2 python=2.7` (Windows users: do *not* do this from within git
bash! Do it from the anaconda prompt), then `conda activate py2` and install various python running
environments `conda install -y -c anaconda jupyter ipython spyder`

``conda activate base``
to switch to the base (python 3)
environment,
type ``ipython -pylab`` 
to make sure a python 3 version number is listed when ipython opens (and type
``exit`` to quit ipython),
and that numpy and matplotlib load properly.
(If they do not, try getting the `-c conda-forge` packages instead, and you can use
``python -c "import matplotlib.pyplot;print 'successful test'"`` as a more
rapid means of testing)

``conda activate py2``
to switch to the python 2
environment,
and
type ``ipython -pylab`` 
to make sure a python 3 version number is listed when ipython opens (and type
``exit`` to quit ipython),
and that numpy and matplotlib load properly.
(If they do not, try getting the `-c conda-forge` packages instead, and you can use
``python -c "import matplotlib.pyplot;print 'successful test'"`` as a more
rapid means of testing)


You may get
[this error](https://github.com/conda/conda/issues/5448) when trying to open ipython,
and needed to close the terminal window and open it again, then switch to py2
environment, in order for ipython to load.

## set up bash so it can switch environments

`conda activate base` to switch back to the base distribution
`conda update -y python-libarchive-c` then `conda init bash` -- it might
ask for admin permissions

On a Windows computer with git installed:
*   edit `~/.bashrc` to make sure `/c/ProgramData/Anaconda3` is in the path
    and that `/c/ProgramData/Anaconda3/MinGW/bin` is in the path
    *and* that you do not have an alias called `conda` defined.
*   Add the following lines to the end of `~/.bash_profile`:

    ~~~sh
    eval “$('/C/ProgramData/Anaconda3/Scripts/conda.exe' 'shell.bash' 'hook')”
    export ps1tmp=$(echo $PS1 | sed 's/\\n\>//g' )
    export PS1=“$ps1tmp ”
    ~~~

verify that you can `conda activate py2` inside bash and see a 2.7
number displayed under ipython

verify that you can `conda activate base` inside bash and that it
switches the version displayed under ipython

at this stage, you can switch to running commands inside bash for everyday use,
though we strongly recommend using the anaconda prompt for installation commands
(`python setup.py`, etc.)

## install pySpecData into py2 environment

Windows users: perform all of the following steps in an anaconda prompt (not bash).

``conda activate py2``

(in the py2 environment) install pySpecData prerequisites from the documentation into the py2
environment:
`conda install -y -c anaconda numpy scipy sympy pyqt pytables matplotlib h5py libpython mingw`

Mac users: In a manner similar to before, from terminal ``cd /opt/anaconda3/env``,
``open .`` to open a Finder window at this location, and locate the py2 folder.
Right-click on this, and allow 'Read & Write' privileges to everyone if it is
not allowed already.

in the pySpecData git distro, check out a python 2 branch (as of
this writing, that will be py2)

Check that “which gcc” (bash) or “where gcc” (dos) points to a command
inside the anaconda `envs\py2` folder (if you experience a 127 error
during linking, it's due to this issue)

install pySpecData in the python 2 environment
``python setup.py develop`` (if you have installed before on this
computer, could be good to add “`build_ext –force`” to the end of this
command line, after deleting the “build” subdirectory, just to be sure)

## Using pySpecData

Note that whenever you run pySpecData, the git branch that's checked out must match the version of python that you are currently running
(if you have activated the base environment, you must have a python 3
pySpecData checked out; if you have activated the py2 branch, you must have a
python 2 branch checked out).
This is because we are using the “develop” setup, where python is actually reading the modules out of our git repo (which is good, because then you just pull the most recent version of pySpecData to update the code).

