ACERT spectrometer data
=======================

ACERT functions for processing spectrometer data

.. (comment) here, I use a search and replace to generate the list below, and then running sphinx-autogen -o generated acert_hdf.rst generates stub files for all my functions 
    I have to run make twice for the result to build
    I also sed replaced all the titles to remove the module name
.. (comment) 
    initially had a 
    :toctree: generated
    line under autosummary, which I used to generate the stubs, but then realize I don't want all these in the toctree

.. currentmodule:: pyspecdata.acert_hdf5
.. autosummary::
    :toctree: generated

    automagical_phasecycle
    cw
    diagnostic_plot
    echo_display
    find_attenuation
    find_file
    gen_composite
    load_and_format
    load_nutation_curve
    open_cw_file
    oscilloscope_data
    plot_coherence_diagram
    plot_comparison
    plot_oned_v_field
    plot_oned_v_offset
    postproc_b1_fid_file
    postproc_blank
    postproc_echo_T2
    postproc_eldor_3d
    postproc_eldor_file
    postproc_eldor_old
    postproc_generic
    search_freed_file
    show_coherence_pathway
    show_pathways

.. (comment) the :members: properties is needed so it pulls all the members of the class
