from pyspecdata.h5nmr import *
import textwrap
grab_data_directory()
###########################{{{
# change the following
name = 'dna_cs19_unbound_120418' # replace with the name of the experiment
chemical_name = 'dna_cs19_unbound' # this is the name of the chemical i.e. 'hydroxytempo' or 'DOPC', etc.
run_number = 120418 # this is a continuous run of data
SL_concentration = 200e-6 # the concentration of spin label
dontfit = False # only set to true where you don't expect enhancement
extra_t1_problem = True # always start with this false, and if it complains, about the number of T1 powers not matching with the number of T1 experiments, then turn it to True
bad_t1_autovals = [] # leave this empty if you're not excluding any. If you want to toss exp x and y out, [x,y]
t1_starttime = 4.9 # start time for the T1 power log (in minutes)
t1_steprange = [5.0,40] # the range (in minutes) of possible times for the T1 experiment
starttime = 4.6 # start time for the E(p) power log
extra_time = 6.0 # used to slice up the E(p) experiment -- the extra time per scan
phcycsteps_separate = True # are the phase cycled steps stored separately?
path = 'reference_data/nmr/' # the path where the data is stored, relative to the data directory you selected during startup
refresh_database = True # this option clears the data in the database generally, could be a good idea to leave this as False, unless you have reloaded the data
###########################}}}
# leave the rest of the code relatively consistent
if refresh_database:
    search_delete_datanode('dnp.h5',name)
if phcycsteps_separate:
    extra_kwargs = {'phnum':[4],'phchannel':[-1]}
else:
    extra_kwargs = {}
#{{{ generate the powers for the T1 series
print('\n\nNext, check the $T_1$ powers:\n\n')
fl = figlistl()
power_threshold = -35
t1_dbm,fl.figurelist = auto_steps(getDATADIR()+path+name+'/t1_powers.mat',
    threshold = power_threshold,t_minlength = t1_steprange[0]*60,
    t_maxlen = t1_steprange[1]*60, t_start = t1_starttime*60.,
    t_stop = inf,first_figure = fl.figurelist)# note that auto_steps is old-style, so I pull the list out of the figlist class
t1mask = bool8(ones(len(t1_dbm)))
# the next line will turn off select (noisy T1
# outputs) enter the number of the scan to remove --
# don't include power off
if extra_t1_problem == True:
    t1mask[-1] = 0
t1_autovals = r_[2:2+len(t1_dbm)] # so that I can exclude experiments
for j in bad_t1_autovals:
    t1mask = logical_and(t1mask,t1_autovals != j)
fl = check_autosteps(power_threshold,t1_dbm,figure_list = fl,mask = t1mask)
fl.show('T1powers_'+name+'.pdf')
print(r't1\_dbm is:',lsafen(t1_dbm),'\n\n')
#}}}
dnp_for_rho(getDATADIR()+path,name,integration_width = 160,
        peak_within = 500, show_t1_raw = True,
        t1_autovals = r_[2:2+len(t1_dbm)][t1mask],
        t1_powers = r_[t1_dbm[t1mask],-999.],
        power_file = name+'/power.mat',t_start = starttime,
        chemical = chemical_name,
        concentration = SL_concentration,
        extra_time = extra_time,
        dontfit = dontfit,
        run_number = run_number,
        threshold = -50.,
        **extra_kwargs)
print(r'\subparagraph{Noise test}')
standard_noise_comparison(name,path = path)
