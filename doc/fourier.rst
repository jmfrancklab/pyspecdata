#. Set the `FT` property appropriately, and changing the units of the
   axis appropriately.

#. Identify whether the :math:`u` or :math:`v` domain is the original
   “source” of the signal, and which is derived from the source. By
   default assume that the source is not aliased. In this way, I can
   automatically marked whether an axis is assumed to be “safe” (i.e.
   “not aliased”) or not. This is relevant when performing time-
   (frequency-)shifts that are not integral multiples of
   :math:`\Delta u` (:math:`\Delta v`).

#. Pull the :math:`u`-axis and determine :math:`\Delta u`.

#. Calculate (only) the initial :math:`v`-axis (starting at 0), the
   post-shift needed to get the :math:`v`-axis I want (based on
   `FT_start_v`), any full-SW aliasing that’s needed to get the
   :math:`v`-axis that I want, and the post-shift discrepancy (not
   accounted for by the combination of the integral post-shift and the
   full-SW aliasing) needed to get the :math:`v`-axis I want (this
   becomes a :math:`u`-dependent linear phase shift)

#. If I do a traditional shift (*i.e.*, like `fftshift`), mark as
   `FT_v_not_aliased`.

#. If there is a post-shift discrepancy, check that the :math:`u`-axis
   is balanced, then apply the post-shift discrepancy as a linear phase
   shift along :math:`u`.

#. Zero-fill before any pre-shifting, since zeros should be placed at
   large positive frequencies.

#. Pre-shift to place the origin at the beginning of the axis, aliasing
   the negative frequencies to the right of the largest positive
   frequencies, and store any discrepancy that’s not accounted for by a
   the shift (which is an integral multiple of :math:`\Delta u`).

#. Perform the FFT and replace the axis with the initial :math:`v`.

#. Since the data has already been :math:`v`-shifted using the linear
   phase shift along :math:`u` above, change the :math:`v`-axis to
   reflect this.

#. Adjust the normalization of the data (this depends on whether we are
   doing `.ft()` or `.ift()`).

   -  As of now, the ft\ :math:`\Rightarrow`\ ift is not invertible,
      because I have defined them as an integral over the exponential
      only; I might later consider dividing the ft by the record length
      (:math:`T`) to return the original units.

#. If there was any pre-shift discrepancy, apply it as a phase-shift
   along the :math:`v`-axis.
