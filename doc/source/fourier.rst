Fourier Functions
=================

#. Use the `FT` property to check that I'm not trying to ft frequency-domain
   data or to ift time domain data.  To start with the data is marked as
   “neither” by having `FT` set to `None`, and the first operation marks it as
   time- or frequency-domain data by setting `FT` to `False` or `True`,
   respectively.

#. Identify whether the :math:`u` or :math:`v` domain is the original
   “source” of the signal, and which is derived from the source. By
   default assume that the source is not aliased. In this way, I can
   automatically marked whether an axis is assumed to be “safe” (i.e.
   “not aliased”) or not. This is relevant when performing time-
   (frequency-)shifts that are not integral multiples of
   :math:`\Delta u` (:math:`\Delta v`).

#. Get the previous :math:`u`-axis and change the units of the axis
   appropriately, *i.e.* s→(cyc/s), (cyc/s)→s.

#. Determine the padded length of the new axis.  If keyword argument `pad` is
   simply set to `True`, round up to the nearest power of 2, otherwise set it
   to `pad`.

#. Use :math:`\Delta u`.

#. Use :math:`\Delta u` and the padded length to calculate (only) the initial
   :math:`v`-axis, which starts at 0.  Then calculate:

   #. Any full-SW aliasing that’s needed to get the :math:`v`-axis that I want.

   #. The (sub-SW) post-transform-shift needed to get the :math:`v`-axis I want (based on `FT_start_v`). This is broken down into:

         #. An integral “post-transform-shift” that will be applied *after* the transform.
         #. A “post-transform-shift discrepancy” (between 0 and 1), which will be applied
                     as a :math:`u`-dependent phase shift *before* the
                     transform.

   - If I apply a traditional shift (*i.e.*, like `fftshift`), mark as
     `FT_[v]_not_aliased` (where *[v]* is time or frequency), *i.e.* “safe,” since the
     :math:`v`-domain is balanced about 0.

#. If there is a post-transform-shift discrepancy, deal with it before I start to mess with the :math:`u`-axis:

   #. check that the :math:`u`-axis is “safe”

   #. apply the post-transform-shift discrepancy as a linear phase shift along :math:`u`.

#. Zero-fill before any pre-transform-shifting, since zeros should be placed at
   large positive frequencies;
   :math:`u` needs to be re-calculated here
   based on original starting :math:`u` and :math:`\Delta u`.

#. Since :math:`u` doesn't necessarily start at 0, calculate the pre-transform-shift
   that's needed to make it start at 0.  Then, apply the integral part of the
   pre-transform-shift, and store the pre-transform-shift discrepancy, which will be applied as a
   phase shift along :math:`v`.

   - Note that the negative frequencies to the right of the largest positive
     frequencies.

#. Perform the FFT and replace the axis with the initial :math:`v`.

#. Apply the post-transform-shift:

   #. Apply the (previously stored) integral part of the post-transform-shift.

   #. Apply the (previously stored) full-SW aliasing.

   #. If the data has experienced a non-integral :math:`v`-shift (*i.e.*
         non-zero post-transform-shift discrepancy) using the linear
         phase shift along :math:`u` above, change the :math:`v`-axis to
         reflect this.

#. Adjust the normalization of the data (this depends on whether we are
   doing `.ft()` or `.ift()`).

   -  As of now, the ft\ :math:`\Rightarrow`\ ift is not invertible,
      because I have defined them as an integral over the exponential
      only; I might later consider dividing the ft by the record length
      (:math:`T`) to return the original units.

#. If there was any pre-transform-shift discrepancy, apply it as a phase-shift
   along the :math:`v`-axis.
