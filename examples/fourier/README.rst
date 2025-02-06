Fourier Transform Functionality
===============================

In this discussion, we refer to the two domains as *time-like* and *frequency-like*, regardless of whether their units match those of time and frequency. The *Fourier transformation (FT)* takes us from the time-like domain to the frequency-like domain, while the *inverse Fourier transformation (IFT)* reverses this process.

To deal cleanly but flexibly with Fourier transformations and the associated axis coordinates, we remember the following facts:

- The relationship between the number of points, the acquisition time, the spectral width, and the spacing between points in both domains is fixed by the Nyquist theorem.
- The signal in both the time-like and frequency-like domains is infinitely periodic, repeating:
  
  - Once every acquisition length (plus one sampling interval) in the time-like domain.
  - Once every spectral width in the frequency-like domain.

- On the computer, we choose a window over which to view these infinitely periodic signals. Because the width of this viewport is fixed by the spacing between our axis coordinates, we simply need to choose a start point for our viewport. The two most common examples are:
  
  - Starting at an axis coordinate of zero.
  - Choosing a start point such that zero frequency is exactly in the middle of our axis coordinates, especially in the frequency-like domain.

- There are other cases where one might want to choose a start point that neither corresponds to zero nor to the perfectly symmetric shifted case. For example, for echo-like data in magnetic resonance, one typically wants to set the origin of the time coordinates to the center of the echo, which in general is neither at the exact beginning of the signal nor the exact center of the signal.

- The choice of different start points can be achieved, quite trivially, by an axis-coordinate-dependent phase shift that is applied subsequent to the Fourier transformation or inverse Fourier transformation. For example, in NMR, shifting the origin of the time axis coordinates results in a frequency-dependent phase shift.

  - If the choice of start point is an integral multiple of the spacing between the axis coordinate points, this will always work, because the axis-coordinate-dependent phase shift that we need to apply after we Fourier transform completes an integral number of cycles from the beginning to the end of our viewport on the data -- in other words, from the start point through to the start point plus the acquisition time or the spectral width.

  - **An important caveat** is that if the start point is a non-integral multiple of the spacing between the axis coordinate points, then we must know that the domain that we are transforming into does not contain data that has been aliased.
    
    - For example, if the domain that we are transforming into is the domain in which we originally acquired signal, then we are almost definitely OK, because that data seldom involves the superposition of multiple Fourier aliases.

    - However, say that we have a free induction decay that includes two signals, one of which exceeds the sampling rate. Then, when we try to Fourier transform, the signal whose frequency exceeds the sampling rate will be aliased and will show up as a peak along the frequency axis coordinate that is the difference of a multiple of the Nyquist frequency and the true frequency of the signal. When we apply our frequency-dependent phase shift, we will get in trouble because the wrong phase will be applied to this aliased signal.

  - Therefore, if we have explicit knowledge that a domain does not include aliased signal, we must label it as such. If not, pyspecdata will throw an error when attempting to apply the relevant corrective phase.

