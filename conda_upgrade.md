These instructions are for installing Anaconda with PySpecData,
both a Python 2 and Python 3 versions.

Please note that if you downloaded the “Python 2.X” version of anaconda,
that means that your package manager, *etc.*, use Python 2.X.
As of this writing (1/13/20), those versions of the package manager are poorly
maintained, and lead to a host of headaches.
So, *yes*, the solution is to *completely uninstall* Anaconda, and then reinstall.

Note that you are *strongly encouraged* to use the Python 3.X version of
PySpecData, and to report any errors that you find.
You can convert your old scripts to Python 3 using ``2to3 -w name_of_your_python_script.py``

### set up anaconda

The basic strategy here is to install the most recent version of
anaconda (3.7) so that we are using an up-to-date “conda” package
manager, but then install a virtual environment where we can run old
Python 2 code if we so choose.

**If at any point during these instructions, you don't get the expected result, stop and seek help!**

#### initial

remove current anaconda

reboot to complete uninstall

go [here](https://www.anaconda.com/distribution/) to download the Python
3.7 version of Anaconda, and install

#### checking for applocker

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

#### changing permissions to allow package installation

make the Anaconda3 folder writeable by all users.
On Windows, you achieve this by
opening the
`C:\ProgramData` folder in File Explorer right click on
`Anaconda3`→ security tab→ edit→
users→ click “full control” checkbox on bottom →
apply (this takes a few minutes to run)

#### install pyspecdata inside py3 environment

make sure the git distro is set to a python 3 branch (as of this
writing, master is python 2, and there is a py3 branch, but that will
change)

install various python running environments
`conda install -y -c anaconda jupyter ipython spyder`

install pyspecdata prerequisites from the documentation:
`conda install -y -c anaconda numpy scipy sympy pyqt pytables matplotlib h5py libpython mingw`

install pyspecdata in the python 3 (base) environment
`python setup_paramset.py install`

if this gives syntax errors, remove the pyspecdata/build directory

and `python setup.py develop`

#### create a python 2 environment, and install basic tools 

`conda create –name py2 python=2.7` (do *not* do this from within git
bash! Do it from the anaconda prompt), then `conda activate py2` and install various python running
environments `conda install -y -c anaconda jupyter ipython spyder`

``conda activate base``
to switch to the base (python 3)
environment,
type ``ipython`` to make sure a python 3 version number is listed when ipython opens (and type ``exit`` to quit ipython)

``conda activate py2``
to switch to the python 2
environment,
type ``ipython`` to make sure a python 2 version number is listed when ipython opens (and type ``exit`` to quit ipython)

#### set up bash so it can switch environments

`conda activate base` to switch back to the base distribution
`conda update -y python-libarchive-c` then `conda init bash` -- it might
ask for admin permissions

On a Windows computer with git installed:
*   edit `~/.bashrc` to make sure `/c/ProgramData/Anaconda3` is in the path
    and that `/c/ProgramData/Anaconda3/MinGW/bin` is in the path
    *and* that you do not have an alias called `conda` defined.
*   Add the following lines to the end of ` /.bash_profile`:

    ~~~sh
    eval “$('/C/ProgramData/Anaconda3/Scripts/conda.exe' 'shell.bash' 'hook')”
    export ps1tmp=$(echo $PS1 | sed 's/\\n\>//g' )
    export PS1=“$ps1tmp ”
    ~~~

verify that you can `conda activate py2` inside git bash and see a 2.7
number displayed under ipython

verify that you can `conda activate base` inside git bash and that it
switches the version displayed under ipython

at this stage, you can switch to running commands inside bash if you
like

#### install pyspecdata into py2 environment

``conda activate py2``

(in the py2 environment) install pyspecdata prerequisites from the documentation into the py2
environment:
`conda install -y -c anaconda numpy scipy sympy pyqt pytables matplotlib h5py libpython mingw`

make sure the pyspecdata git distro is set to a python 2 branch (as of
this writing, master is python 2, but that will change)

check that “which gcc” (bash) or “where gcc” (dos) points to a command
inside the anaconda `envs\py2` folder (if you experience a 127 error
during linking, it's due to this issue)

install pyspecdata in the python 2 environment
``python setup_paramset.py install``

and ``python setup.py develop`` (if you have installed before on this
computer, could be good to add “`build_ext –force`” to the end of this
command line, after deleting the “build” subdirectory, just to be sure)
