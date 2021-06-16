"""
Frequency and Time Shifting
===========================

Here we show that relabeling an axis before or after a Fourier Transform
generates the expected result.
"""
# from JF notebook sec:task4126
from pylab import *
from pyspecdata import *
fl = figlist_var()
t = r_[2.5:12:1024j]
data = nddata(empty_like(t),[-1],['t']).setaxis('t',t)
data.set_units('t','s') # set the units to s, which are automatically converted to Hz upon FT
data = data.fromaxis('t',lambda x: where(logical_and(x > 3,x < 6),1,0))

default_plot_kwargs = dict(alpha = 0.5, linewidth = 2)

fl.next('time domain -- positive starting point')
fl.plot(data,**default_plot_kwargs)
fl.plot(data.runcopy(imag),**default_plot_kwargs)
expand_y()

fl.next('and shift $t_{start}\\rightarrow 0$')
data.ft('t',shift = True)
#control = data.copy()
data.ft_clear_startpoints('t',t = 0,f = 'current')
data.ift('t')
fl.plot(data,**default_plot_kwargs)
fl.plot(data.runcopy(imag),**default_plot_kwargs)
expand_y()

#fl.next('diagnose frequencies')
#fl.plot(abs(control),label = 'control',**default_plot_kwargs)
#data.ft('t')
#fl.plot(abs(data),label = 'after ift',**default_plot_kwargs)
t = r_[0:10:32j]
dt = t[1] - t[0]
t -= 3.5*dt # to ensure that no point passes through zero
data = nddata(empty_like(t),[-1],['t']).setaxis('t',t)
data.set_units('t','s') # set the units to s, which are automatically converted to Hz upon FT
data = data.fromaxis('t',lambda x: where(logical_and(x > 3,x < 6),1,0))


default_plot_kwargs.update(dict(marker = 'o'))

fl.next('strange time axis with nothing passing through zero')
fl.plot(data,label = 'R: original',**default_plot_kwargs)
fl.plot(data.runcopy(imag),label = "I: original",**default_plot_kwargs)
expand_y()

data.ft('t',shift = True)
data.ft_clear_startpoints('t',t = 0,f = 'current')
data.ift('t',pad = 1024)
default_plot_kwargs.update(dict(marker = None))
fl.plot(data,label = "R: new",**default_plot_kwargs)
fl.plot(data.runcopy(imag),label = "I: new",**default_plot_kwargs)
expand_y()

fl.show('ft_demo_weird_startpoints_151030.pdf')
