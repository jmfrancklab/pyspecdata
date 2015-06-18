from pyspecdata.fornotebook import *
name = '100729_400nm_140CAT1_VitC'
power_file = '400nm_vitC'
path = '/home/ahammack/Data/'
dbm = auto_steps(path+power_file,
   t_minlength = 2.8*60, # minimum length of each step is 2.5 minutes
   t_maxlen = 3.5*60, # set a max length of 3.5 minutes, after which it will go back to the largest jump in power, and cut off there
   tolerance = 2, # tolerance before it declares a new step is tolerance x the standard deviation
   t_start = 4.1*60) # in this instance, I'm not starting until 4.1 minutes in, since the first step is -35
dbm = dbm[:-6] # eliminate the powers for the last 6 scans, since they are crap
dbm = r_[-999,-35,dbm] # add in the zero power, and manually decided 35 dBm steps
figure(1)
lplot('powerlog_raw'+name+'.pdf')
figure(2)
lplot('powerlog'+name+'.pdf')
dnp_for_rho(path,name,dbm,
   expno = r_[5,9:25], # the experiments to process $\Rightarrow$ in this instance, would be r_[5:34], except, I don't have a power for the first 3 with power, and also the three low-power scans at the end, then would be r_[5,9:31], but I eliminated the last 6
   h5file = None, # don't save in an HDF file
   t1expnos = r_[4,35,36], # the T1 experiments to process
   integration_width = 300,
   peak_within = 1e3)
