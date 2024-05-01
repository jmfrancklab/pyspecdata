History/Roadmap
---------------

(Current version in bold -- future plans below) 

**0.9.5.6**
    - Updated sphinx gallery with DCCT plots
      accompanying domain colored coherence transfer
      preprint.
0.9.5.5
    - Implement a reader for Cary (Varian) UV-Vis files with examples.
0.9.5.1
    - 0.9.5.1.1
      Some important debugging, and also added `pyspecdata.ipy` → executing the following at the top of a jupyter notebook:

        .. code-block:: python

            %pylab inline
            %load_ext pyspecdata.ipy

      will cause nddata to "display" as labeled plots.

    - 0.9.5.1.2
      added ability to load power saturation 2D data from Bruker

    - 0.9.5.1.3
      XEpr data loaded with dBm units rather than W units

      added ``to_ppm`` function for Bruker files

    - 0.9.5.1.4
      Improved internal logging, and started to remove gratuitous dependencies,
      ``%load_ext pyspecdata.ipy`` includes
      ``%pylab inline``, so that only

        .. code-block:: python

            %load_ext pyspecdata.ipy

        is required for jupyter.
  
    - 0.9.5.1.6

        - Removed several legacy modules, and added docstrings for the remaining modules.

        - Begin phasing out earlier `CustomError` class.

        - Make `numpy` pretty printing available from the `general_functions` module.

        - Add xelatex support to the notebook wrapper.

        - Start to move file search routines away from demanding a single "data directory."

        - Improved support for 2D Bruker XEPR

        - Made it possible to call standard trig functions with `nddata` as an argument.
    - 0.9.5.1.7
        - ILT (Tikhonov regularization) with SVD Kernel compression
          (1 and 2 dimensions)
        - ``smoosh`` and ``chunk`` deal with axes properly

0.9.5.3 **Current Version**
    upgrade to Python 3 and begin to flesh out documentation

0.9.5.4
    - 0.9.5.4.1
      - ``to_ppm`` should only be a method of inherited class
      - 1.5 and 2.5 D ILT
0.9.5
    First version distributed on pypi.python.org.

Future Plans
^^^^^^^^^^^^

0.9.5.7
    - Implement an ``nddata_placeholder`` class for quickly loading and
      searching through datasets in *e.g.* UV-Vis files or Bruker directories
      without actually loading all the data from each dataset.
1.0
    We are working on four major upgrades relative to the 0.9 sequence:

    - Axes as objects rather than a set of separate attributes of nddata.
    - Remove dependence on pytables in favor of h5py.
    - Replace figure lists with “plotting contexts,” which will still
      enable PDF vs. GUI plotting, but will better integrated with Qt and
      more object-oriented in nature
    - Comma-separated indexing to work correctly with all indexing types.
      (0.9.5 requires sequential brackets rather than comma-separated
      indexing for some combined range selections.)

1.0.2
    GUI for setting configuration directories.

    Means for dealing with non-linearly spaced data in image plots
    (0.9.5 auto-detects log spacing in 1D plots,
    but pretends that image plots are linear -- we will implement linear spline
    interpolation algorithm)

1.0.3
    Bruker DSP phase correction for raw data from newer versions of Topspin that is in sync with the code from nmrglue.

1.0.4
    Package a make-less copy of lapack to allow a cross-platform build of density matrix propagation routines.

1.1.0
    Integrate with ACERT NLSL Python package for simulation and fitting of ESR spectra.

1.2.0
    Implement a version of figure list that can be interfaced with Qt.

