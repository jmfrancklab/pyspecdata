from pyspecdata.h5nmr import *
grab_data_directory()
fl = figlistl()
noe_roe = {'ubq_F4C':-0.5,
    'ubq_T14C':-0.41,
    'ubq_T66C':-0.35,
    'ubq_T66C_N14':-0.35,
    'ubq_S57C':-0.27,
    'ubq_G35C':0,
    'ubq_R72C':-0.50}
# rerun
data,data_nice = retrieve_DNP_set(['ubq_T66C','ubq_F4C','ubq_G35C','ubq_T66C_N14','ubq'],divide_klow_by = 6.0,verbose = True,
    t10subst = {'ubq_T66C':'ubq','ubq_F4C':'ubq','ubq_G35C':'ubq','ubq_T66C_N14':'ubq'})
#{{{ first, just plot the nice data
data_nice = lambda_rec(data_nice,'run',lambda x: '%0.1f'%x,'run_number') # make a text field
data_nice = drop_fields(data_nice,['run_number'])
print('after formatting nicely, and taking the average of available $T_{1,0}$ measurements:\n\n')
lrecordarray(data_nice,resizebox = 0.5)
print('\n\n')
fl.next('bargraph')
clf()
textlabel_bargraph(data_nice,spacing = 0.1,tickfontsize = 5)
axis('tight')
expand_x()
ylim((-0.1,1.0))
autolegend() # redo this, where I apply the labels
print('\n\n')
#}}}
#{{{ now, do some manual analysis
print('data is:\n\n')
data = sorted(reorder_rec(data,['chemical','run_number']))
lrecordarray(data,resizebox = 0.8)
print('\n\n')
print('First, I average over the ones with undecided \\tonen\n\n')
print('and add the coupling factor\n\n')
data = applyto_rec(mean,data,['run_number','chemical','fit_type'])
data = lambda_rec(data,
    r'\xi',
    (lambda x,y: x/y),
    ['ksmax','krho'])
lrecordarray(data,resizebox = 0.8)
print('\n\nNow plot the \\ksigma vs. NOE/ROE\n\n')
noe_axis = r'$\sigma_{NOE}/\sigma_{ROE}$'
mask = data['run_number'] == 130410.2
mask = logical_and(mask,data['chemical'] == 'ubq_T66C_N14')
for thisdata in [
        ('ksigma',r'$k_\sigma$ / s$^{-1}$M$^{-1}$',data['ksmax'],data['chemical']),
        ('krho',r'$k_\rho$ / s$^{-1}$M$^{-1}$',data[~mask]['krho'],data[~mask]['chemical']),
        ('xi',r'$\xi$',data[~mask][r'\xi'],data[~mask]['chemical']),
        ('tau',r'$\tau_c$ / ps',interptau(data[~mask][r'\xi'],14.85,simple = True)/1e-12,data[~mask]['chemical'])]:
    ksmax_data = nddata(thisdata[2],[-1],[noe_axis]).labels(noe_axis,
        array([noe_roe[j] for j in thisdata[3]]))
    fl.next(thisdata[0]+'_correlation')
    plot(ksmax_data,'o')
    ylabel(thisdata[1])
    expand_x()
    expand_y()
    ksmax_data.plot_labels([j[4:].replace('_N14',' ($^{14}$N)') for j in thisdata[3]])
#{{{ Manual analysis of data
print('\n\nNow take the mean and drop a bunch of fields\n\n')
data = drop_fields(data,['run_number','fit_type'])
data_mean_new = meanstd_rec(data,['chemical'])
sfo1 = 14.85
data_mean_new = lambda_rec(data_mean_new,
    r'\tau',
    lambda x: interptau(x,sfo1,simple = True),
    [r'\xi'])
#{{{ calculate tau error by finite differences
#{{{ calculate one thousandth of xi then use it to calculate derivative by finite differences
data_mean_new = lambda_rec(data_mean_new,
    r'\tau_ERROR',
    lambda x: x*1e-3,
    [r'\xi'])
onethousandthofxi = data_mean_new[r'\tau_ERROR'].copy()
tau_at_plusonethousandthofxi = interptau(data_mean_new[r'\xi']+onethousandthofxi,sfo1,simple = True)
data_mean_new[r'\tau_ERROR'] = abs(tau_at_plusonethousandthofxi - data_mean_new[r'\tau'])/onethousandthofxi
#}}}
data_mean_new[r'\tau_ERROR'] *= data_mean_new[r'\xi_ERROR']
#}}}
#{{{ scale things to reasonable units
data_mean_new = lambda_rec(data_mean_new,
    r'\xi',
    (lambda x: x/0.01))
data_mean_new = lambda_rec(data_mean_new,
    r'\xi_ERROR',
    (lambda x: x/0.01))
data_mean_new = lambda_rec(data_mean_new,
    'concentration',
    (lambda x: x/1e-6))
data_mean_new = lambda_rec(data_mean_new,
    r'\tau',
    (lambda x: x/33.3e-12))
data_mean_new = lambda_rec(data_mean_new,
    r'\tau_ERROR',
    (lambda x: x/33.3e-12))
data_mean_new = rename_fields(data_mean_new,{r'\xi':r'\xi/0.01',
    r'\xi_ERROR':r'\xi/0.01_ERROR',
    r'\tau':r'\retardationfactor',
    r'\tau_ERROR':r'\retardationfactor_ERROR',
    'ksmax':r'\ksm',
    'ksmax_ERROR':r'\ksm_ERROR',
    'klow':r'\klow',
    'klow_ERROR':r'\klow_ERROR'})
    #}}}
data_mean_new = drop_fields(data_mean_new,['T_{1,0}','T_1','concentration_ERROR','T_1_ERROR','T_{1,0}_ERROR','krho','krho_ERROR','concentration'])
lrecordarray(data_mean_new,resizebox = 0.8,std_sf = 2)
#}}}
#}}}
fl.show('ubiquitin_summary130425.pdf')
