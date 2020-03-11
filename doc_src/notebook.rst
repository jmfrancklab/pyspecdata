LaTeX Notebook Functionality
============================

This package contains tools for running a LaTeX notebook with
embedded python code.  One can then, for instance, keep an
electronic lab notebook where the plots and data are generated
and processed in place, immediately from the *raw data*.  (The author of the project has kept an
electronic notebook for several years in this way with great
success.)

Please note that this is a very different "notebook" than the one
supplied by jupyter notebooks (which pyspecdata also supports).
In this section, we aim to provide a laboratory notebook with
publication-ready figures and a complete path from raw data to those
figures.

* With minimal effort, it is possible to design a notebook that
  works in a similar fashion with HTML, Markdown, *etc.*  (The main
  author just has no interest in doing this, since the PDF output
  looks very nice.)
* It's highly recommended to keep a notebook as a series of
  different files that have no preamble or ending (*i.e.* only
  the part that goes inside the ``document`` environment), which
  can then be collected into a gigantic master document with the
  ``\input{...}`` command, or compiled individually (*e.g.* while
  you are actually in the lab working on a particular section).
* It's also highly recommended to store notes organized by
  project, which can then be cross-referenced in a separate
  chronological document (or *vice versa*) with the
  ``\\ref{...}`` command.

Setting up the notebook
-----------------------

Requirements
------------

1.  latex packages

    You need to put the texmf tree in a location where it can be found by your latex installation.

    * Under Windows, add the texmf tree to "miktex settings" under the "roots" tab.
    * (To do) We should use `this guide
      <http://ctan.math.washington.edu/tex-archive/info/dtxtut/dtxtut.pdf>`_ or
      or `this package <https://ctan.org/pkg/makedtx>`_ to package the code and include it here.

    Once you've done this, the shell command ``kpsewhich mypython.sty``
    should return a result
    (if you have miktex installed on windows, this should work from either the git
    bash prompt or the dos or powershell prompt).
2.  The pyspecdata package.
    
    Proves the commands `pdflatex_notebook_wrapper` and
    `update_notebook_pythonscripts`, described below under "Running the
    notebook."

    Also provies the command `pdflatex_notebook_view_wrapper`, which is used to
    determine the output PDF and call an appropriate viewer.
3.  A standard latex compilation system:

    You can use latexmk (shipped with miktex) with `Sumatrapdf <https://www.sumatrapdfreader.org/free-pdf-reader.html>`_
    (Sumatrapdf allows you to edit the PDF while it's open in Sumatrapdf, while Adobe Acrobat *does not*).
    Here is a ``~/.latexmkrc`` file that works on windows:

    .. code-block:: perl

        $pdflatex=q/pdflatex_notebook_wrapper %O -synctex=1 --xelatex %S/;
        $pdf_previewer=q/pdflatex_notebook_view_wrapper/;#calls the wrapviewer function

    **It should also be possible to use TeXworks** by adding pdflatex_notebook_wrapper to
    preferences → typesetting → processing tools.

    .. todo::

        Alec, can you check this out and update documentation??
4.  The `'paramset_pyspecdata'` module.
    Just run ``python setup_paramset.py install`` from repository root directory (the same directory where you run
    ``python setup.py install``
    or
    ``python setup.py develop``)

    (This external module is simply used
    to store the context in which the code is called -- *i.e.*, from within
    python *vs.* from the command line.)

5.  It's assumed that your latex files are stored in a "notebook directory."
    In some cases, during the first run, an explanatory error will appear -- just follow the instructions.

Running the notebook
--------------------

If you create figures with the `figlist_var` class,
you should simply be able to write a latex file with embedded
``python`` environments (``\\begin{python}`` ... ``\\end{python}``)
replace the ``pdflatex`` command with
``pdflatex_notebook_wrapper`` when compiling your latex notebook,
to drop the code and plots in place.
For clarity, the code output is a slightly different color (a
dark brown) than the standard text.

A synctex "jump to source" on the resulting portion of the PDF
will send you to the tex output, which is stored in
``scripts/*.tex``, where ``*`` is a sequential number
corresponding to the script, and the python source used to
generate it is stored in ``scripts/*.py``.

Each snippet of unique code is run **only once**, **ever** making the
notebook fast and efficient.
For now, an important drawback to this is that if you change modules or
libraries called by the script, the resulting PDF output will not
change.
To get around this, a command is provided that forces scripts to
be re-run.  You use it like this: ``update_notebook_pythonscripts
flush 10 21``  -- which will flush script numbers 10 to 21.
Manually deleting the ``.py`` files inside the scripts directory
will **not** have the same effect.

.. todo:: 
    To limit downtime for the PDF, pdflatex_notebook_wrapper currently copies
    the final pdf to a truncated filename (assuming that the filename consists
    of words separated by underscores, it drops the last word).

    It would be much better to copy the source PDF into a subdirectory, build it there, and then copy the pdf back into the main directory.
    This would entail changing the paths of the various files

    ``\RequirePackage[log]{snapshot}`` might be helpful to log files here.

    probably we will just want to add commands to renewcommand for input as well as the graphicx root.


How it works
------------

* Note that the code works a bit differently than in previous
  versions -- while it previously required LaTeX to be run with shell-escape enabled,
  ``pdflatex_notebook_wrapper`` doesn't require this.
* ``pdflatex_notebook_wrapper`` just calls ``pdflatex`` followed
  by ``update_notebook_pythonscripts`` 

The LaTeX end
`````````````

The file `mypython.sty` looks for the `python` environment, it pulls the relevant
code, outputs it to ``scripts/*.py`` and then writes a command to
the ``.aux`` file that tells LaTeX where to find the
``scripts/*.tex`` output.
The ``scripts/*.tex`` output is only updated once
``update_notebook_pythonscripts`` (without arguments) is run.

The python end
``````````````

``update_notebook_pythonscripts`` runs through the various
``scripts/*.py`` files, checks and checks whether or not they
have been previously cached (under the same or a different script
number).  If the python code is in the cache, it just pulls the
cached output.  If not, it runs the file, and stores the result
in the cache.

.. toctree::
    :maxdepth: 2

    pyspecdata.rst
    latexscripts.rst
