LaTeX Notebook Functionality
============================

This package contains tools for running a LaTeX notebook with
embedded python code.  One can then, for instance, keep an
electronic lab notebook where the plots and data are generated
and processed in place, immediately from the *raw data*.  (The author of the project has kept an
electronic notebook for several years in this way with great
success.)

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

Each snippet of unique code is run **only once**, making the
notebook fast and efficient.
For now, an important drawback to this is that if you change modules or
libraries called by the script, the resulting PDF output will not
change.
To get around this, a command is provided that forces scripts to
be re-run.  You use it like this: ``update_notebook_pythonscripts
flush 10 21``  -- which will flush script numbers 10 to 21.
Manually deleting the ``.py`` files inside the scripts directory
will **not** have the same effect.

How it works
------------

* Note that the code works a bit differently than in previous
  versions -- while it previously required LaTeX to be run with 
* ``pdflatex_notebook_wrapper`` just calls ``pdflatex`` followed
  by ``update_notebook_pythonscripts`` 

### The LaTeX end

The file `mypython.sty` looks for the `python` environment, it pulls the relevant
code, outputs it to ``scripts/*.py`` and then writes a command to
the ``.aux`` file that tells LaTeX where to find the
``scripts/*.tex`` output.
The ``scripts/*.tex`` output is only updated once
``update_notebook_pythonscripts`` (without arguments) is run.

### The python end

``update_notebook_pythonscripts`` runs through the various
``scripts/*.py`` files, checks and checks whether or not they
have been previously cached (under the same or a different script
number).  If the python code is in the cache, it just pulls the
cached output.  If not, it runs the file, and stores the result
in the cache.
