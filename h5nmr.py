from matlablike import *
from nmr import *
from nmrfit import *
from interptau import interptau
#from fornotebook import lrecordarray
import time
from numpy.lib.recfunctions import rename_fields,drop_fields
bulk_klow_name = r'k_{low}/k_{low,bulk}'
bulk_ksmax_name = r'k_{\sigma}/k_{\sigma,bulk}'
bulk_xi_name = r'\xi/\xi_{bulk}'
bulk_krho_value = 353.4
bulk_xi_value = 0.27
bulk_correltime = 54e-12
bulk_ksigma_value = bulk_xi_value * bulk_krho_value
bulk_klow_value = 5./3.*bulk_krho_value - 7./3.*bulk_ksigma_value
#{{{ Classes to load and store stuff to the HDF5 file
##{{{ the parent
class store_integrals:
    def __init__(self,h5filename,
            chemical_id,
            fid_args,
            first_figure = None,
            run_number = double(0),
            auxiliary_args = {},integration_args = {},integrate_function = integrate,
            lplotfigurekwargs = {},
            type_str = 'integral'):
        self.figurelist = figlistini_old(first_figure)
        self.run_number = run_number
        self.type_str = type_str
        self.name = fid_args["name"]
        self.compilationroot_name = '/compilations'
        self.expno = fid_args['expno']
        self.integralnode_name = genexpname(self.name,fid_args['expno']) # this will be "auto" for Emax, but never T1
        self.chemical_id = chemical_id
        self.h5filename = h5filename
        ###{{{ save the Emax integrals only as necessary
        try:
            #print lsafen('DEBUG: Looking for the integral data')
            h5file,datanode = h5nodebypath(self.h5filename + '/integrals' + '/' + self.integralnode_name,check_only = True)
            obs('Found the integral data',lsafe(self.integralnode_name),' in the HDF5 file, not re-integrating\n\n')
            h5file.close()
            self.newrecord = False
            self.catalog()
            # next, load the data into nddata
            #print lsafen('DEBUG: Successfully found the integral data and closed the file')
            return
        except CustomError:
            #if 'test.h5' in listdir('.'):
            #    print 'DEBUG outer: file exists'
            #else:
            #    print 'DEBUG outer: file does not exist'
            #print '\n\n'
            ####{{{ load the self and accompanying data
            self.auxiliary(fid_args,auxiliary_args) # tack on powers, delay times, etc.
            self.integral,self.figurelist = integrate_function(dirformat(fid_args["path"])+fid_args["name"]+'/',self.expno,first_figure = self.figurelist,pdfstring = self.integralnode_name,**integration_args)
            if type(self.expno) in [list,ndarray]:
                begin_and_end_files = [self.expno[0],self.expno[-1]]
            else:
                begin_and_end_files = [self.expno]
            for thisexp in begin_and_end_files:
                thisfile = dirformat(fid_args["path"]) + dirformat(fid_args["name"])+'%d'%thisexp
                thisfiledate = load_acqu(thisfile)['DATE']
                thisfiledate = time.localtime(thisfiledate)
                thisfiledate = time.strftime('%d %b %Y %I:%M:%S %p',thisfiledate)
                print '\n\n',r'{\bf Start time of %s exp %d}'%(type_str,thisexp),thisfiledate,'\n\n'
            ####}}}
            self.integral.set_prop('name',self.integralnode_name) # sets the name for subsequent saving
            try:
                self.attach_aux()
            except:
                lplotfigures(self.figurelist,'error_figures.pdf')
                raise CustomError('Look at the power meter graph to see if you can reset the extra_time parameter, or if there is something wrong with the power meter.')
            self.integral.hdf5_write(self.h5filename + '/integrals')
            self.newrecord = True
            self.catalog()
            return
        ###}}}
    def auxiliary(self,fid_args,auxiliary_args):
        pass
        return
    def attach_aux(self):
        pass
        return
    def catalog(self):
        pass
        return
##}}}
##{{{ store Emax
class store_Emax (store_integrals):
    def catalog(self):
        h5file,compilationroot_node = h5nodebypath(self.h5filename + '/compilations')
        search_string = '(integrals == \'%s\')'%self.integralnode_name 
        if self.newrecord:
            Emax_table,self.index = h5addrow(compilationroot_node,'Emax_curves',[long(self.chemical_id),self.integralnode_name,expand_expno(self.expno),self.run_number],['chemical_id','integrals','experiments','run_number'],match_row = search_string)
            print 'DEBUG loaded row:\n\n',lrecordarray(Emax_table.readWhere('index == %d'%self.index),smoosh=True),'\n\n'

        else: # these things are probably not defined yet.
            Emax_table,self.index = h5addrow(compilationroot_node,'Emax_curves',[],[],match_row = search_string)
            self.expno = expand_expno(Emax_table.readWhere('index == %d'%self.index)['experiments'].flatten())
            print 'DEBUG read row\n\n',lrecordarray(Emax_table.readWhere('index == %d'%self.index)),'\n\n'
        h5file.close()
        return
    def auxiliary(self,fid_args,auxiliary_args):
        ###{{{ auxiliary part to load the power
        #print lsafen('DEBUG: I didn\'t find the self data, but it\'s important that the file is closed now.')
        obs('Didn\'t find the self data in the HDF5 file, pulling integrals from the FID\n\n')
        path = fid_args['path']
        if 'dbm' in fid_args.keys(): # if I have manually specified the dbm values
            self.expno = fid_args['expno']
            dbm = fid_args['dbm']
        else:
            power_file = auxiliary_args.pop('power_file')
            if fid_args['expno'] == 'auto' or fid_args['expno'] == 'Auto':
                self.expno,dbm,self.figurelist = parse_powers(path,fid_args['name'],power_file,first_figure = self.figurelist,**auxiliary_args) # the 1 is because I'm not interested in the expnos it uses
                fl = figlist_var()
                fl.figurelist = self.figurelist
                mythreshold = -36
                if 'threshold' in auxiliary_args.keys():
                    mythreshold = auxiliary_args['threshold']
                fl = check_autosteps(mythreshold,dbm,figure_list = fl)
                self.figurelist = list(fl.figurelist)
                del fl
            else:
                self.expno = fid_args['expno']
                junk,dbm,self.figurelist = parse_powers(path,fid_args['name'],power_file,expno = self.expno,first_figure = self.figurelist,**auxiliary_args) # the 1 is because I'm not interested in the expnos it uses
                fl = figlist_var()
                fl.figurelist = self.figurelist
                mythreshold = -36
                if 'threshold' in auxiliary_args.keys():
                    mythreshold = auxiliary_args['threshold']
                fl = check_autosteps(mythreshold,dbm,figure_list = fl)
                self.figurelist = list(fl.figurelist)
                del fl
            print ''%(),r'\begin{minipage}{3in}','\npower log:\n\n',', '.join(map((lambda x: r'$%0.2f\;dBm$'%x),dbm)),r'\end{minipage}','\n\n'
            dbm = r_[-999,dbm] # 6/20/11, I get rid of the -10 here to make things less confusing, though I now for sure don't know what's right
        print 'len(dbm)',len(dbm)
        print '\n\nstart integrate\\Emax\n\n'
        self.powers= dbm_to_power(dbm)
        return
    def attach_aux(self):
        self.integral.labels(['power'],[self.powers])
        return
        ###}}}
##}}}
##{{{ store T1
class store_T1 (store_integrals):
    def __init__(self,h5filename,
            chemical_id,
            fid_args,
            first_figure = None,
            run_number = double(0),
            auxiliary_args = {},integration_args = {},integrate_function = process_t1,
            lplotfigurekwargs = {'numwide':4,'width':2}):
        store_integrals.__init__(self,h5filename,
            chemical_id,
            fid_args,
            first_figure = first_figure,
            run_number = run_number,
            auxiliary_args = auxiliary_args,integration_args = integration_args,integrate_function = integrate_function,lplotfigurekwargs = lplotfigurekwargs,type_str = r'$T_1$')
        return
    def auxiliary(self,
            fid_args,
            auxiliary_args):
        if 't1power' not in auxiliary_args:
            self.expno = fid_args['expno']
            firstfile = format_listofexps([dirformat(fid_args["path"])+fid_args["name"],self.expno])[0]
            ###{{{ search for the self.power in the title string
            titlestr = load_title(firstfile)
            a = re.compile('(?:^| |\$)([0-9\.]+).*dB')
            m = a.search(titlestr)
            fordebugstring = 0.0
            if m:
                g = m.groups()
                self.power = -10.0-double(g[0]) # 6/20/11 this is some type of guess based on the input -- I should calibrate it later
                fordebugstring = double(g[0])
            else:
                print "\n\n\\o{found no self.power in the title titlestr}"
                self.power = -999.
            #print lsafen('DEBUG: saved with power',self.power,'because the power string is',fordebugstring)
            ###}}}
        else:
            self.power = auxiliary_args['t1power']
        if 'temperature' in auxiliary_args:
            self.temperature = auxiliary_args['temperature']
        else:
            self.temperature = -1.
    def catalog(self):
        if self.newrecord:
            mydata = [long(self.chemical_id),self.integralnode_name]
            mynames = ['chemical_id','integrals']
            mycoeff = self.integral.output()
            mynames.append('run_number')
            mydata.append(self.run_number)
            mynames += ['power']
            mydata += [self.power]
            mynames += ['temperature']
            mydata += [self.temperature]
            mynames += ['num_coeff']
            mydata += [len(mycoeff.dtype.names)]
            mynames += mycoeff.dtype.names
            mydata += mycoeff.tolist()[0]
            covarmatrix = self.integral.covarmat().flatten()
            mynames += ['covarelement%d'%x for x in range(0,len(covarmatrix))]
            mydata += list(covarmatrix)
        h5file,compilationroot_node = h5nodebypath(self.h5filename + '/compilations')
        search_string = '(integrals == \'%s\')'%self.integralnode_name 
        if self.newrecord:
            try:
                concentration_table,self.index = h5addrow(compilationroot_node,'T1_fits',mydata,mynames,match_row = search_string)
            except:
                raise CustomError("Problem saving row",lsafen(zip(mynames,mydata)))
            #print 'DEBUG loaded row\n\n',lrecordarray(concentration_table.readWhere('index == %d'%self.index),smoosh = False),'\n\n'
        else:
            concentration_table,self.index = h5addrow(compilationroot_node,'T1_fits',[],[],match_row = search_string)
            #print 'DEBUG read row\n\n',lrecordarray(concentration_table.readWhere('index == %d'%self.index),smoosh = False),'\n\n'
        h5file.close()
        return
##}}}
#}}}
#{{{ standard power parsing
def parse_powers(path,name,power_file,expno = None,
        scanlength_estimation_error = 0.2,
        t_minlength = 10.0, # in seconds
        t_start = 4.9, # in minutes
        t_stop_extra = 0.0, # in minutes
        extra_time = 7.86, # amount extra per scan
        tolerance = 2.0,
        first_figure = None,
        threshold = -36):
    figurelist = figlistini_old(first_figure)
    fullpath = dirformat(dirformat(path)+name)+'%d/'%1
    fileinfo = load_acqu(fullpath)
    if expno == None:
        expno = r_[5:5+3+1+int(fileinfo['CNST'][6])] # get the number of experiments stored by the au program
        obs(r'automatically determined that you\'re using expno %d $\Rightarrow$ %d'%(expno[0],expno[-1]))
    try:
        fullpath = dirformat(dirformat(path)+name)+'%d/'%expno[0]
    except:
        raise CustomError('type of expno is',expno.dtype)
    fileinfo = load_acqu(fullpath)
    if det_type(fullpath)[1] == True:
        td1 = load_acqu(fullpath,whichdim='2')['TD']
    else:
        td1 = 1
    scan_length = td1*fileinfo['NS']*(fileinfo['D'][1]+fileinfo['TD']/fileinfo['SW_h']/2.0+fileinfo['D'][11]) + 11.0 + extra_time
    t_stop = t_start*60.0+t_stop_extra*60.0+(len(expno)-1)*scan_length
    obs(r'scan\_length is',scan_length,'and there are %d scans for t$_{stop}$=%0.3f\n\n'%(len(expno)-1,t_stop/60.0))
    scanlength_estimation_error *= scan_length # convert from percentage to time
    templen = scan_length-scanlength_estimation_error
    if templen > t_minlength:
       t_minlength = templen
    t_maxlen = scan_length+scanlength_estimation_error
    dbm,lastspike,figurelist = auto_steps(path+power_file,t_start = t_start*60,t_stop = t_stop,t_minlength = t_minlength,t_maxlen = t_maxlen,tolerance = tolerance, threshold = threshold,return_lastspike = True,first_figure = figurelist)
    downby = r_[-5:5]
    timetomakeup = (len(expno)-len(dbm))*60.0*scan_length+lastspike
    obs(lsafen(r"Assuming you have less scans than you need or just one too many, you are down by expno-1-dbm=%d and the last spike is %0.2f s from the end, which means your scan length is off by %0.2f, set extra_time from %0.2f to %0.2f"%(len(expno)-1-len(dbm),lastspike,timetomakeup,extra_time,extra_time + timetomakeup/(len(expno)-1))))
    #obs('If you have one scan too many, set extra\_time to',extra_time+lastspike/(len(expno)-1),'$s$\n\n')
    return figlistret(first_figure,figurelist,expno,dbm)
#}}}
#{{{ functions to pull info out of file names
##{{{ chemical
def pull_chemical_from_filename(name):
    a = re.compile('(.*?)_?([0-9]+)([mu]M)')
    b = re.compile('(.*?)(_|[0-9])')
    m = a.search(name)
    if m:
        g = m.groups()
        chemical_name = g[0]
        print "\n\n\\o{found chemical name} \\verb|",g[0],'|\n\n'
    elif b.search(name):
        chemical_name = b.search(name).groups()[0]
        print "\n\n\\o{found chemical name} \\verb|",chemical_name,'|\n\n'
    else:
        chemical_name = 'not_found'
        print "\n\n\\o{found no chemical name}\n\n"
    return chemical_name
##}}}
##{{{ concentration
def pull_concentration_from_filename(name):
    a = re.compile('([0-9]+)_?([mu]M)')
    m = a.search(name)
    if m:
        g = m.groups()
        if g[1]=='uM':
            concentration = 1e-6
        elif g[1]=='mM':
            concentration = 1e-3
        concentration *= double(int(g[0]))
    else:
        concentration = 0
    return concentration
##}}}
##{{{ date
def pull_date_from_filename(name):
    a = re.compile('([0-9]{2})([0-9]{2})([0-9]{2})')
    m = a.search(name)
    g = m.groups()
    date_of_file = mktime(datetime(int(g[0])+2000,int(g[1]),int(g[2])).timetuple())
    return date_of_file
##}}}
#}}}
#{{{ Convenience functions for different nodes
def search_delete_datanode(filename,name = None,chemical = None,concentration = None):
    'calls delete_datanode on all nodes matching the search criteria'
    if (chemical is None) and (concentration is None) and (name is not None):
        try:
            h5file,childnode =  h5nodebypath(filename+'/integrals/',check_only = True)
        except:
            print lsafen("The node /integrals/ doesn't exist in",filename,": continuing")
            return
        names_of_children = childnode._v_children.keys()
        for thisname in names_of_children:
            if thisname.find(name) > -1:
                delete_datanode(filename,thisname)
        h5file.close()
    else:
        chemical_ids = get_chemical_index(filename+'/compilations/chemicals',chemical = chemical)
        for chemical_id in chemical_ids:
            # find the nodes that contain this chemical
            searchstring = 'chemical_id == %d'%chemical_id
            for thistablename in ['Emax_curves','T1_fits']:
                h5file,compilationroot_node =  h5nodebypath(filename+'/compilations',check_only = True)
                thistable = h5table(compilationroot_node,thistablename,None)
                thistable_result = thistable.readWhere(searchstring)
                integral_names = thistable_result['integrals']
                print '\n\n'+r'search\_delete\_datanode found datasets: '+lsafen(', '.join(list(integral_names)))+'for chemical %d'%chemical_id+lsafen('(%s)'%chemical)
                h5file.close()
                for j in integral_names:
                    delete_datanode('dnp.h5',j)
def delete_datanode(filename,integralnodename):
    "Deletes a data node (by name) and all references to it"
    try:
        h5file,childnode =  h5nodebypath(filename+'/integrals/'+integralnodename,check_only = True)
    except:
        print lsafen("The node /integrals/",integralnodename,"doesn't exist in",filename,": continuing")
        return
    h5file.removeNode(childnode,recursive = True)
    print lsafe('removed node',integralnodename)
    compilationroot_node = h5file.root.compilations
    for thistablename in ['Emax_curves','T1_fits']:
        numrem,removed_table_data = h5remrows(compilationroot_node,thistablename,"integrals == '%s'"%integralnodename)
        if numrem:
            print lsafe('and %d associated lines in %s'%(numrem,thistablename))
            if thistablename is 'Emax_curves':
                for j in removed_table_data['index']:
                    for fittablename in ['xismax_fits','Emax_fits']:
                        numrem,junk = h5remrows(compilationroot_node,fittablename,"Emax_curve == %d"%j)
                        print lsafe('and %d associated lines in %s'%(numrem,'xismax_fits'))
    h5file.close()
    print lsafen()
    return
def print_chemical_by_index(compilationroot_name,indexnumber):
    h5file,chemicalnode = h5nodebypath(compilationroot_name)
    if type(indexnumber) in [list,ndarray]:
        lrecordarray(chemicalnode.chemicals.readWhere(h5inlist('index',indexnumber)),format = '%g')
    else:
        lrecordarray(chemicalnode.chemicals.readWhere('index == %d'%indexnumber),format = '%g')
    h5file.close()
    return
def get_chemical_index(compilationroot_name,*args,**kwargs):
    r'OLD format: get_chemical_index(location of HDF5 node, chemical name,concentration)' + '\n' + 'OR get_chemical_index(location of HDF5 name, string with chemical + conc)' + '\n' + 'NEW FORMAT: pass no arguments, but pass kwargs chemical, concentration, and name'
    ##{{{ manual kwargs
    verbose = False
    if 'verbose' in kwargs.keys():
        verbose = kwargs.pop('verbose')
    ##}}}
    #{{{ parse both old and new format arguments
    chemical_name = None
    concentration = None
    name = None
    only_last = True
    if len(kwargs) > 0:
        if 'only_last' in kwargs.keys():
            only_last = kwargs.pop('only_last')
        if 'chemical' in kwargs.keys():
            chemical_name = kwargs.pop('chemical')
        if 'concentration' in kwargs.keys():
            concentration = kwargs.pop('concentration')
        if 'name' in kwargs.keys():
            name = kwargs.pop(name)
        if len(kwargs) > 0:
            raise CustomError('kwargs not understood!:',kwargs)
    else:
        if len(args) == 1: # a filename was given
            raise CustomError('You should NO LONGER be letting this function determine the concentrations!!  Determine them externally, then pass to this function!!!')
            if verbose: print 'DEBUG: only one arg given',lsafen(args)
            chemical_name = pull_chemical_from_filename(args[0])
            concentration = pull_concentration_from_filename(args[0])
        elif len(args) == 2:
            chemical_name = args[0]
            concentration = args[1]
    #}}}
    if chemical_name == None:
        raise CustomError('Chemical name is "None"!')
    ##{{{ identify the id in the solutions table corresponding to the chemical and concentration
    if concentration is None:
        search_string = '(chemical == \'%s\')'%chemical_name 
    else:
        search_string = '(chemical == \'%s\') &'%chemical_name + gensearch('concentration','%0.5g',concentration,1e-6)
    ###{{{ grab (or create) the table
    h5file,compilationroot_node = h5nodebypath(compilationroot_name)
    concentration_table,concentration_id = h5addrow(compilationroot_node,'chemicals',[chemical_name,double(concentration)],['chemical','concentration'],match_row = search_string,verbose = True,only_last = only_last)
    h5attachattributes(compilationroot_node.chemicals,['TITLE'],['(indexed list of solutions)'])
    if verbose: print 'concentration table (returned row %d):\n\n'%concentration_id,lrecordarray(concentration_table.read()),'\n\n'
    h5file.close()
    ###}}}
    ##}}}
    return concentration_id
def expand_expno(expno):
    "Converts the expno list in the HDF5 file into a numpy array"
    expno = list(expno)
    if -1 in expno:
        return array([j for j in expno if j != -1])
    else:
        retval = -1*ones(200,dtype = 'int')
        retval[0:len(expno)] = expno
        return retval
def genexpname(name,expnos):
    "Generates the node name associated with a series of experiments in a standard way"
    if type(expnos) is str:
        retval = name + '_%s'%expnos
    elif isscalar(expnos):
        retval = name + '_exp%d'%(expnos)
    else:
        retval = name + '_exp%dto%d'%(expnos.min(),expnos.max())
    return retval
#}}}
#{{{ base functions for dealing with dnp data
##{{{ linear dnp plot
def emax_linearandplot(integral,first_figure = None,pdfstring = '',power_axis = 'power',max_invpower = inf,color_pair = ['r','k'],show_Evip_plot = True,**kwargs):
    "generates the 1/(1-E) plots for DNP"
    figurelist = figlistini_old(first_figure)
    lineardata = integral.linear()
    lineardata = lineardata['1 / power',lambda x: x < max_invpower]
    x = lineardata.getaxis('1 / power').copy()
    if show_Evip_plot:
        ###{{{ show the linear plot
        nextfigure(figurelist,'Evip'+pdfstring)
        plot_updown(lineardata,'1 / power',color_pair[0],color_pair[1],nosemilog = True)
        ax = gca()
        ax.set_xlim([0,array(ax.get_xlim()).max()])
        Evipplot = integral.linear(r_[1/array(ax.get_xlim()).max(),1/x.min()],**kwargs)
        plot(Evipplot,color_pair[0]+'-',alpha = 0.1)
        ###}}}
    return figlistret(first_figure,figurelist)
##}}}
##{{{ DNP for rho
def dnp_for_rho(path,
        name,
        dbm = None,
        expno = [],
        h5file = 'dnp.h5',
        t10name = None,
        t10chemical = None,
        t10conc = 0.0,
        t1powers = [], # for this and the next, remember that the no power scan is always added to the end
        t1_powers = None, # new format
        t1expnos = [],
        t1_autovals = r_[2,3], 
        power_file = None,
        guessonly = False,
        chemical = None,
        concentration = None,# the spin label concentration
        show_t1_raw = True,
        t1_offset_corr = None,
        dontfit = False,
        t1_phnum = None,
        t1_phchannel = None,
        t1_bad_powers = None,
        expnos = None,
        clear_nodes = False,
        first_figure = None,
        #simul_types = ['ernst angle','ernst angle clipped','5 x T1'], # simulate an ernst angle experiment, etc
        simul_types = [], # simulate an ernst angle experiment, etc
        pdfstring = '',ignore_emax_error = False,
        color_pair = ['k','r'],
        run_number = 0,
        run_number_for_conc = None,
        plot_vs_A = False,
        ksp_symbol = '',
        verbose = False,
        T1w = 2.4691,
        DeltaT1w = 0.10267,
        **kwargs):
    #{{{ pull T1max
    temp = load_acqu(dirformat(dirformat(path)+name)+'5',return_s = False)
    d1val = temp['D'][1]
    temp = load_acqu(dirformat(dirformat(path)+name)+'1',return_s = False)
    t1max = temp['CNST'][9]
    sfo1 = temp['SFO1']
    print '\n\n'
    #}}}
    #{{{ inital kwarg processing
    run_number = double(run_number)
    print '{\\bf\\color{blue}This data is part of series %f:}\n\n'%run_number
    if t1_powers is not None:
        t1powers = t1_powers
    if expnos is not None and expno == []:
        expno = expnos
    figurelist = figlistini_old(first_figure) # when I convert this to a new function, definitely convert to a new format
    if (type(t1powers) is list and len(t1powers) > 0) or (type(t1powers) is ndarray):
        t1powers = array(t1powers,dtype = 'double')
    if type(clear_nodes) is list:
        clear_nodes = array(clear_nodes)
    if type(clear_nodes) is ndarray:
        for j in clear_nodes:
            delete_datanode(h5file,genexpname(name,j))
        clear_nodes = False
    if clear_nodes:
        if len(expno) > 0:
            delete_datanode(h5file,genexpname(name,expno))
        if len(t1expnos) > 0:
            for j in t1expnos:
                delete_datanode(h5file,genexpname(name,j))
    if expno == [] and t1expnos == []:
        expno = 'auto'
        t1expnos = 'auto'
        if clear_nodes:
            delete_datanode(h5file,genexpname(name,expno))
            #}}}
    ###{{{ process the kwargs to divvy them up onto the functions they refer to
    fid_args = {'name':name,'path':path,'expno':expno}
    if dbm is not None:
        fid_args.update({'dbm':dbm})
    arglist = ['expno','t1expnos']
    fid_args.update(dict([(x,kwargs.pop(x)) for x in kwargs.keys() if x in arglist]))
    arglist = ['temperature']
    t1_auxiliary_args = dict([(x,kwargs.pop(x)) for x in kwargs.keys() if x in arglist])
    integration_args = {}
    arglist = []
    for func in [integrate]:
        arglist += getargspec(func)[0][-len(getargspec(func)[-1]):]
    integration_args.update(dict([(x,kwargs.pop(x)) for x in kwargs.keys() if x in arglist]))
    if verbose: print lsafen('DEBUG integration args:',integration_args)
    if power_file is None:
        auxiliary_args = {'power_file':'mat_data/'+name+'.mat'}
    else:
        auxiliary_args = {'power_file':power_file}
    arglist = []
    for func in [parse_powers]:
        arglist += getargspec(func)[0][-len(getargspec(func)[-1]):]
    auxiliary_args.update(dict([(x,kwargs.pop(x)) for x in kwargs.keys() if x in arglist]))
    if len(kwargs) > 0:
        raise CustomError('I didn\'t understand the arguments',kwargs)
    ###}}}
    ###{{{ find the concentration id
    if chemical is None:
        if verbose: print 'DEBUG: found no chemical name, so pulling from name\n\n'
        if concentration is not None:
            raise ValueError('You passed no chemical name, but you did pass a concentration.  If you passed neither, I would try to pull from the file, but now I don\'t know what to do!!')
        chemical = pull_chemical_from_filename(name)
        concentration = pull_concentration_from_filename(name)
    else:
        if concentration is None:
            if verbose: print lsafen('DEBUG: pulling from chemical %s and getting the concentration from previously stored ESR data'%(chemical))
            if run_number_for_conc is None:
                run_number_for_conc = run_number
            if type(run_number_for_conc) is array:
                run_number_for_conc = run_number_for_conc.tolist()
            if type(run_number_for_conc) is not list and type(run_number_for_conc) is not slice:
                run_number_for_conc = [run_number_for_conc]
            h5fp,chemical_node = h5nodebypath(h5file+'/compilations/chemicals')
            try:
                thisnode = chemical_node.default_chemical_id
            except:
                raise CustomError('You passed a me a chemical name but no concentration, which would seem to imply you want me to pull a concentration from the table generated by ESR (default_chemical_id), but I can\'t find this table!')
            temp = "(chemical == '%s') & %s"%(chemical,h5inlist('run_number',run_number_for_conc))
            #print "search string is:",temp
            thistable = thisnode.readWhere(temp)
            concentration_id = thistable['chemical_id']
            if len(concentration_id) == 0:
                thistable = thisnode.readWhere("(chemical == '%s')"%(chemical))
                if len(thistable) > 0:
                    raise ValueError('There is no ESR data for '+chemical+' for runs '+repr(run_number_for_conc)+' -- you might try one of the following: '+repr(thistable['run_number']))
                else:
                    raise ValueError('There is no ESR data for %s'%(chemical))
            else:
                thisnode = chemical_node.chemicals
                concentration = thisnode.readWhere(h5inlist('index',concentration_id))['concentration'].mean()
            h5fp.close()
            print "\n\nI automatically identified a concentration of $%0.3f\uM$ based on the ESR of %s for run number %s\n\n"%(concentration/1e-6,chemical,repr(run_number_for_conc))
        else:
            if verbose: print lsafen('DEBUG: pulling from chemical %s and conc %f'%(chemical,concentration))
    concentration_id = get_chemical_index(h5file+'/compilations/chemicals',chemical,concentration)
    ###}}}
    ###{{{ store the data as necessary
    #{{{ retrieve the store_Emax object emax_in_file
    emax_figlist_start = len(figurelist)
    figurelist.append({'print_string':r'\subparagraph{retrieve Emax data}'+'\n\n'})
    if len(expno) > 0:
        emax_in_file = store_Emax(h5file,concentration_id,
                fid_args,
                first_figure = figurelist,# because I append this later!
                auxiliary_args = auxiliary_args,
                integration_args = integration_args)
        figurelist = emax_in_file.figurelist
        if verbose: print lsafen('DEBUG: figurelist after Emax load is',figurelist)
        if verbose: print lsafen('DEBUG: index is',emax_in_file.index,'for',emax_in_file.integralnode_name,'with experiments',emax_in_file.expno)
    emax_figlist_stop = len(figurelist)
    #}}}
    ####{{{ grab any possible t1 args
    t1kwargs = {}
    for argname in ['phchannel','phnum','offset_corr']:
        if 't1_' + argname in integration_args.keys() and (isscalar(integration_args['t1_' + argname]) or len(integration_args['t1_'+argname])>0):
            t1kwargs[argname] = integration_args.pop('t1_'+argname) # does the part before "or" even work???
    ####}}}
    ####{{{ loop through the T_1 experiments, processing and saving as necessary
    figurelist.append({'print_string':r'\subparagraph{retrieve $T_1$ data}'+'\n\n'})
    integration_args.update(t1kwargs) # copy over the t1 args
    if t1expnos == 'auto':
        t1expnos = r_[emax_in_file.expno[-1]+t1_autovals,304] #right now, set to three, but can easily add others later
        t1names_forerror = map(lambda x: 'autoval %d'%x,t1_autovals)+['exp 304']
    else:
        t1names_forerror = map(lambda x: 'manually entered exp %d'%x,t1expnos)
    if len(t1expnos) > 0:
        if clear_nodes:
            for j in t1expnos:
                delete_datanode(h5file,genexpname(name,j))
    #{{{ manually enter certain kwargs
    integration_args_save = dict(integration_args)
    integration_args.update({'show_image':show_t1_raw})
    if t1_offset_corr is not None:
        integration_args.update({'offset_corr':t1_offset_corr})
    if t1_phnum is not None:
        integration_args.update({'phnum':t1_phnum})
    if t1_phchannel is not None:
        integration_args.update({'phchannel':t1_phchannel})
    #}}}
    if len(t1powers) > 0 and len(t1powers) != len(t1expnos):
        raise CustomError("You didn't pass the same number t1powers (%d) as there are T1's (%d)!"%(len(t1powers),len(t1expnos)))
    auxiliary_args.update(t1_auxiliary_args)
    for j,t1expno in enumerate(t1expnos):
        #{{{ should be able to get rid of this line by including with t1 auxiliary args above
        if len(t1powers) > 0:
            auxiliary_args.update({"t1power":t1powers[j]})
        elif "t1power" in auxiliary_args.keys():
            auxiliary_args.pop("t1power")
        #}}}
        fid_args.update({'expno':t1expno})
        if t1expno == 304:
            figurelist.append({'print_string':'power off\n\n'})# t1expnos and t1_autovals should have the same indeces
        else:# t1expnos and t1_autovals should have the same indeces
            figurelist.append({'print_string':'$T_1$ expno %d\n\n'%(t1expno)})
            if len(t1_autovals) > j:
                figurelist.append({'print_string':'$T_1$ autoval %d\n\n'%(t1_autovals[j])})
            if len(t1powers) > 0:
                figurelist.append({'print_string':'$T_1$ power %g\n\n'%dbm_to_power(t1powers[j])})
        try:
            t1info = store_T1(h5file,concentration_id,
                    fid_args,
                    first_figure = figurelist,
                    auxiliary_args = auxiliary_args,
                    integration_args = integration_args,
                    run_number = run_number)
            figurelist = t1info.figurelist
        except:
            test = list(figurelist)
            lplotfigures(figurelist,fid_args['name']+'.pdf') # in case it craps out, so I can still see the power plots
            raise CustomError("Couldn't load T1 for",t1names_forerror[j],'you may just want to skip it; right now manual t1expnos',t1expnos,'and auto t1_autovals',t1_autovals)
        figurelist.append({'print_string':'end $T_1$ dataset\n\n'})
    if "t1power" in auxiliary_args.keys(): auxiliary_args.pop("t1power")
    integration_args = integration_args_save
    ####}}}
    ###}}}
    if len(expno) > 0: # only do this if I have acquired enhancement scans
        if t10chemical == None:
            if t10name != None:
                t10dataset = retrieve_T1series(h5file,t10name,show_result = 'for $T_{1,0}$')
            elif chemical != None: # find the same chemical at 0 conc
                t10dataset = retrieve_T1series(h5file,None,chemical,0.0,run_number = run_number,show_result = 'for $T_{1,0}$')
            else:
                print "{\\bf Warning!!} I can't figure out what your $T_{1,0}$ data is!"
                t10dataset = None
        else:
            if t10name != None and t10conc != None:
                t10dataset = retrieve_T1series(h5file,t10name,t10chemical,t10conc,show_result = 'for $T_{1,0}$')
            elif t10conc != None:
                t10dataset = retrieve_T1series(h5file,None,t10chemical,t10conc,run_number = run_number,show_result = 'for $T_{1,0}$')
            else:
                print "{\\bf Warning!!} I can't figure out what your $T_{1,0}$ data is, even though you passed a t10 chemical name!"
                t10dataset = None
        ###{{{ search for and retrieve the previously stored T_1(p) values
        t1dataset = None
        if chemical == None and concentration == None:
            t1dataset,M0dataset,Minfdataset = retrieve_T1series(h5file,name,run_number = run_number,show_result = 'for $T_{1}$',t1max = t1max,retcheckval = True)
        else:
            t1dataset,M0dataset,Minfdataset = retrieve_T1series(h5file,None,chemical,concentration,run_number = run_number,show_result = 'for $T_{1}$',t1max = t1max,d1val = d1val,retcheckval = True)
        print '\n\n'
        #{{{ store the T1 data with error, for later error propagation
        if t1dataset is None:
            t1_w_error = NaN
        else:
            t1_at_min_power = t1dataset['power',lambda x: x == x.min()]
            t1_w_error = t1_at_min_power.copy()
            if len(t1_w_error.data) > 1:
                t1_w_error.data = t1_w_error.data.mean()
                temp = t1_w_error.get_error().mean()
                t1_w_error.set_error(temp)
        print 'The $T_1$ dataset looks like this:',lsafen(t1dataset),'and the value at minimum power is',t1_w_error
        print '\n\n'
        #}}}
        #{{{ store the T10 data with error, for later error propagation
        if t10dataset is None:
            t10_w_error = nddata([NaN],[1],['power']).set_error([NaN])
        else:
            t10_at_min_power = t10dataset['power',lambda x: x == x.min()]
            t10_w_error = t10_at_min_power.copy()
            if len(t10_w_error.data) > 1:
                t10_w_error.data = r_[t10_w_error.data.mean()]
                temp = t10_w_error.get_error().mean()
                t10_w_error.set_error(r_[temp])
        print 'The $T_{1,0}$ dataset looks like this:',lsafen(t10dataset),'and the value at minimum power is',t10_w_error
        print '\n\n'
        #}}}
        ###}}}
        ###{{{ calculate the linear function for T_1(p)
        #{{{ Set all the T10 related values
        t10_at_zero_power = 2.5
        if t10dataset is not None:
            t10_at_zero_power = t10dataset['power',0.0].data[0]
        T1w_inv_f = lambda x: nddata(1.0/(T1w+x*DeltaT1w),[-1],['power']).labels('power',x)
        print '\n\nlinear interpolation using $T_{1,0}(0)$ of %0.3f for $k_\\sigma$ analysis, which is'%t10_at_zero_power
        if t10dataset is None:
            print 'just randomly estimated\n\n'
        else:
            print 'taken from previous data\n\n'
        Ck_HH = 1./t10_at_zero_power-1./T1w
        print '\n\nThis implies that $C_{macro}k_{HH}$ is %0.3f\n\n'%(Ck_HH)
        #}}}
        #{{{ T_1 interpolation
        #{{{ if there are any scans that are bad (i.e. close to the signal null, remove them by power here
        if t1_bad_powers is not None:
            mask = zeros(size(t1dataset.getaxis('power')),dtype = 'bool8')
            for j in range(0,len(t1_bad_powers)):
                mask |= abs(t1dataset.getaxis('power') - t1_bad_powers[j]) < abs(0.05*t1_bad_powers[j])
                figurelist.append({'print_string':
                    r"\textbf{I'm leaving out $T_1$'s with power %g, because you marked them as bad (there are %d of these)}"%(t1_bad_powers[j],sum(mask))})
            new_t1dataset = nddata(t1dataset.data[~mask],[sum(~mask)],['power'])
            new_t1dataset.labels('power',t1dataset.getaxis('power')[~mask])
            new_t1dataset.set_error(t1dataset.get_error()[~mask])
            t1dataset = new_t1dataset
        #}}}
        t1_at_zero_power = t1dataset['power',0.0].data[0]
        print '\n\nlinear interpolation using $T_1(0)$ of %0.3f for $k_\\sigma$ analysis\n\n'%t1_at_zero_power
        rho = (1./t1_at_zero_power) - (1./t10_at_zero_power)
        print '\n\nlinear interpolation using $\\rho$ of %0.5f for $k_\\sigma$ analysis\n\n'%rho
        #{{{ find list of independent powers
        list_of_t1_powers = t1dataset.getaxis('power').copy()
        list_of_different_t1_powers = array([])
        for j in range(0,len(list_of_t1_powers)):
            if not any(abs(list_of_different_t1_powers -
                    list_of_t1_powers[j])
                    < 0.05*list_of_t1_powers[j]):
                list_of_different_t1_powers = r_[
                        list_of_different_t1_powers,
                        list_of_t1_powers[j]]
        figurelist.append({'print_string':
            r'I found %d different $T_1$ powers:'%len(list_of_different_t1_powers)
            +lsafen(list_of_different_t1_powers)})
        #}}}
        F_nonlinear = 1./(1./t1dataset.copy().set_error(None)
                - Ck_HH
                - T1w_inv_f(t1dataset.getaxis('power')))
        #{{{ threshold, in case the inverse was too small
        threshold_number = 1./0.5e-3 # note that this is 2000s!
        def threshold_F_nonlinear(myx):
            mask = myx.data > threshold_number # if the inverse was very small
            mask |= myx.data <= 0.0 # if the inverse was negative
            mask |= isnan(myx.data) # if it was exactly zero
            myx.data[mask] = threshold_number
        threshold_F_nonlinear(F_nonlinear)
        #}}}
        if len(list_of_different_t1_powers) > 3:
            F_nonlinear_to_order = 2
        else:
            F_nonlinear_to_order = 1
        c_F_nonlinear,F_nonlinear_line = F_nonlinear.polyfit('power',order = F_nonlinear_to_order)
        print '\n\n$F\\_{nonlinear}$ is',lsafen(F_nonlinear),'and its fit is',F_nonlinear_line
        # now, linearly interpolate w.r.t. power
        # in the same step, convert back to T1(p)
        t1err = t1dataset['power',0:1]
        t1err /= t1err.data[0] # so this way, multiplying by this will just give an nddata with teh correct error.
        DeltaT10 = c_F_nonlinear[1]
        def F_nonlinear_f(x):
            x_data = nddata(x, len(x), ['power'])
            #{{{ evaluate c_F_nonlinear over the nonlinear polynomial
            retval = c_F_nonlinear[0]
            for j in range(1,F_nonlinear_to_order+1):
                retval += x_data**j*c_F_nonlinear[j]
            #}}}
            retval.labels('power',x)
            threshold_F_nonlinear(retval)
            return retval
        t1f = lambda x: t1err/(1./F_nonlinear_f(x) + Ck_HH + T1w_inv_f(x))
        powers_forplot = linspace(0,t1dataset.getaxis('power').max(),100)
        #{{{ just diagnosis
        #{{{ plot F_nonlinear directly
        #{{{ plot stuff
        nextfigure(figurelist,'Fnonlinear' + pdfstring)
        plot(F_nonlinear,'o',label = '$F_{nonlinear}$ points')
        #}}}
        powers_fordiagnose = linspace(0,t1dataset.getaxis('power').max(),100)
        #print '\n\n$T_1(p)$ is',lsafen(t1f(powers_fordiagnose))
        closely_spaced_F = F_nonlinear_f(powers_fordiagnose)
        plot(F_nonlinear,'-',alpha=0.8,label = '$F_{nonlinear}(p)$')
        plot(F_nonlinear_line,'-',alpha=0.8,label = r'$F_{nonlinear,origfit}(p)$')
        plot(closely_spaced_F,'-',alpha=0.8,label = '$F_{nonlinear,interpolated}(p)$')
        title(name)
        autolegend()
        #}}}
        #{{{ plot inverse of F_nonlinear
        nextfigure(figurelist,'invFnonlinear' + pdfstring)
        plot(T1w_inv_f(t1dataset.getaxis('power'))+rho,
                '-',label = '$T_{1,w}^{-1}(p)+\\rho(0)$')
        plot(1.0/F_nonlinear,'o',alpha=0.8,label = '$F_{nonlinear}^{-1}(p)$')
        plot(1.0/F_nonlinear_line,'-',alpha=0.8,label = '$F_{nonlinear}^{-1}(p)$')
        plot(1./t1dataset.copy().set_error(None),
                'o', label = '$T_1^{-1}$')
        plot(powers_forplot,1.0/t1f(powers_forplot).set_error(None),
                '-', alpha = 0.2, label = '$T_1^{-1}(p) = F_{nonlinear}^{-1}(p)+(T_{1,w}^-1+p \\Delta T_{1,w}^{-1})^{-1} + C_{macro} k_{HH}$')
        myones = F_nonlinear.copy()
        myones.data[:] = 1.0
        plot(rho*myones,
                '--',
                label = '$\\rho(0)$')
        plot((1./t10_at_zero_power+rho)*myones,
                '--',
                label = '$T_{1,0}^{-1} + \\rho(0)$')
        expand_y()
        expand_x()
        gridandtick(gca())
        title(name)
        autolegend()
        #}}}
        #}}}
        DeltaT10 = c_F_nonlinear[1]
        ###}}}
        #}}}
        ###{{{ grab data for f and plot it
        if t1dataset == None or t10dataset == None:
            fdata_exists = False # if both aren't present, I can't calculate the leakage
        else:
            fdata_exists = True
        #{{{ do the T1 stuff first, regardless
        figurelist.append({'print_string':r'\subparagraph{$T_1$ data}'})
        nextfigure(figurelist,'t1data' + pdfstring)
        save_color_t1 = plot_color_counter()
        plot(1.0/t1dataset,'o', label = '$T_1^{-1}$ / $s^{-1}$')
        ylabel('various (see legend)')
        title(r'%s, run %0.2f'%(chemical,run_number))
        plot_color_counter(save_color_t1)
        plot(powers_forplot,1.0/t1f(powers_forplot).set_error(None),'-',linewidth = 3, alpha = 0.5)
        nextfigure(figurelist,'M0data' + pdfstring)
        plot_updown(-1*M0dataset,'power','k','r',symbol = 'o', label = r'$-M(0)$',alpha = 0.5)
        plot_updown(Minfdataset,'power','b','g',symbol = 'o', label = r'$M(\infty)$',alpha = 0.5)
        autolegend()
        ylabel('net magnetization')
        expand_y()
        expand_x()
        nextfigure(figurelist,'t1data' + pdfstring)		
        #}}}
        if fdata_exists: # if I can calculate the leakage factor
            save_color_t10 = plot_color_counter()
            plot(1.0/t10dataset,'o', label = '$T_{1,0}^{-1}$ / $s^{-1}$')
            ct10,t10line = t10dataset.polyfit('power')
            t10err = t10dataset['power',0:1]
            t10err /= t10err.data[0] # so this way, multiplying by this will just give an nddata with teh correct error.
            if size(t10dataset.data) < 2:
                print r'{\bf Warning:} $T_{1,0}$  dataset only has a single power point'
                if size(t1dataset.data) > 1:
                    print r'estimating power dependence from $T_1(p)$'
                    t10f = lambda x: 1.0/(1.0/t1f(x) - 1.0/t1f(r_[0.0]) + 1.0/t10dataset['power',0:1])
                else:
                    print r"only get one power point for the data:",lsafe(t1dataset)
                print '\n\n'
            else:
                t10f = lambda x: t10err*(nddata(x,len(x),['power'])*ct10[1]+ct10[0])
            f = lambda x: 1.0 - t1f(x)/t10f(x)
            plot_color_counter(save_color_t10)
            plot(powers_forplot,1.0/t10f(powers_forplot).set_error(None),'-', alpha = 0.2)
            plot(powers_forplot,f(powers_forplot).set_error(None),'-',label = '$f$', alpha = 0.2)
            #{{{ show k_\rho in the plot, along with Delta rho (overall change in rho)
            print '$T_1$ power axis is',lsafen(t1dataset.getaxis('power'))
            t1_with_power_off = t1dataset['power',lambda x: argmin(x)].data
            t1_at_max_power = t1dataset['power',lambda x: argmax(x)].data
            print r'$T_1$ for $\rho$ is',t1_with_power_off
            print r'$T_1$ for $\rho_{max}$ is',t1_at_max_power
            rho_forprint = 1.0/t1_with_power_off
            rho_atmax = 1.0/t1_at_max_power
            t10_with_power_off = t10dataset['power',lambda x: argmin(x)].data.mean()
            t10_at_max_power = t10dataset['power',lambda x: argmax(x)].data.mean()
            print r'$T_{1,0}$ for $\rho$ stored as',t10_with_power_off
            print r'$T_{1,0}$ for $\rho_{max}$ stored as',t10_at_max_power
            rho_forprint -= 1.0/t10_with_power_off
            rho_atmax -= 1.0/t10_at_max_power
            print r'$T_{1,0}$ for $\rho$ is',t10dataset['power',lambda x: argmin(abs(x))].data,'\n\n'
            print r'$T_{1,0}$ for $\rho$ is NOT',t10dataset['power',lambda x: argmin(abs(x+999.))].data,'\n\n'
            print r'$C$ for $k_\rho$ is',concentration,'\n\n'
            print r'$\rho$ is',rho_forprint,'\n\n'
            print r'$\rho_{max}$ is',rho_atmax,'\n\n'
            print r'$k_\rho$ is',rho_forprint/concentration,'\n\n'
            print r'$k_{\rho,max}$ is',rho_atmax/concentration,'\n\n'
            textstring = '$k_\\rho = %0.3f$ s$^{-1}$M$^{-1}$\n$\\Delta T_{1,0}=%0.3f$ s/W\n\Delta\\rho %0.3f %%'%(rho_forprint/concentration,DeltaT10,(rho_forprint-rho_atmax)/rho_forprint)
            ax = gca()
            text(0.5,0.5,textstring,transform = ax.transAxes,size = 'x-large', horizontalalignment = 'center',color = 'k')
            #}}}
        #{{{ set the bottom ylim to 0, which is a more sensible plot
        ax = gca()
        ylims = array(ax.get_ylim())
        ylims[ylims.argmin()] = 0.0
        ylims[ylims.argmax()] *= 1.1
        ax.set_ylim(ylims)
        #}}}
        gridandtick(gca())
        expand_x()
        autolegend()
        ###}}}
        #{{{ grab and normalize the previously stored list of integrals
        emax_unmod = nddata_hdf5(h5file + '/integrals/' + emax_in_file.integralnode_name)
        normalizer = emax_unmod.copy().sort('power')['power',0:1]
        emax_unmod /= normalizer
        #}}}
        #{{{ calculate xi*s(p) by the leakage factor method, with and without heating correction
        if fdata_exists:
            ###{{{ two methods of processing the real data
            xisp_f = (1.0 - emax_unmod.copy()) * 1.51671e-3 / f(r_[0.0])
            xisp_fp = (1.0 - emax_unmod.copy()) * 1.51671e-3 / f(emax_unmod.getaxis('power'))
            ###}}}
        #}}}
        if not dontfit:
            ###{{{ calculate and plot k*s(p), both corrected and uncorrected
            if concentration is None:
                raise ValueError('Concentration cannot be None!  Give at least a bogus concentration that I can use the calculate a relaxivity!')
            ksp_corr = (1.0 - emax_unmod.copy()) / t1f(
                    emax_unmod.getaxis('power')) / concentration * 1.51671e-3 # calculate k s(p)
            ksp_uncorr = (1.0 - emax_unmod.copy()) / t1f(
                    r_[0.0]) / concentration * 1.51671e-3 # calculate k s(p)
            figurelist.append({'print_string':r'\subparagraph{ksp data}'})
            nextfigure(figurelist,'ksp_corr')
            ksp_corr.rename('power','p')
            ksp_uncorr.rename('power','p')
            plot_updown(ksp_uncorr,'p','b','g',symbol = ksp_symbol,label = 'uncorrected')
            plot_updown(ksp_corr,'p','k','r',symbol = ksp_symbol,label = 'corrected')
            ylabel(r'$k_\sigma s(p)$')
            title(r'%s, run %0.2f'%(chemical,run_number))

            ksp_corr = ksp(ksp_corr)
            ksp_uncorr = ksp(ksp_uncorr)
            ksp_corr.fit()
            ksp_uncorr.fit()
            plot(ksp_uncorr.eval(100),'b')
            ax = gca()
            text(0.25,0.75,ksp_uncorr.latex(),transform = ax.transAxes,size = 'x-large', horizontalalignment = 'center',color = 'b')
            plot(ksp_corr.eval(100),'k')
            text(0.75,0.25,ksp_corr.latex(),transform = ax.transAxes,size = 'x-large', horizontalalignment = 'center',color = 'k')
            if verbose:
                print r'\begin{verbatim}'
                print 'DEBUG: ksp_uncorr is',ksp_uncorr
                print 'DEBUG: ksp_corr is',ksp_corr
                print r'\end{verbatim}'
            ###}}}
            ###{{{ tabulate the fits
            h5file_node,compilationroot_node = h5nodebypath(h5file + '/compilations')
            search_string = '(Emax_curve == %d)'%(emax_in_file.index)# this is the index of the Emax_curve, and is unique for each dataset
            numrem,junk = h5remrows(compilationroot_node,'ksp_fits',search_string)
            if numrem:
                print lsafen('Removed',numrem,'rows from ksp_fits, with ksmax',junk[r'ksmax'],'because they conflict with these entries')
            #{{{ add the corrected data
            if not fdata_exists:# t10_with_power_off is undefined
                t10_with_power_off = 0.0
            row_to_add = {'Emax_curve':emax_in_file.index,
                    'run_number':run_number,
                    'fit_type':'corrected',
                    'T_{1,0}':t10_with_power_off,
                    'covariance':ksp_corr.covarmat()}
            mycoeff = ksp_corr.output()
            row_to_add.update(dict(zip(mycoeff.dtype.names,mycoeff.tolist()[0])))
            concentration_table,parameters_index = h5addrow(compilationroot_node,'ksp_fits',row_to_add)
            #}}}
            #{{{ add the uncorrected data
            row_to_add = {'Emax_curve':emax_in_file.index,
                    'run_number':run_number,
                    'fit_type':'uncorrected',
                    'T_{1,0}':t10_with_power_off,
                    'covariance':ksp_uncorr.covarmat()}
            mycoeff = ksp_uncorr.output()
            row_to_add.update(dict(zip(mycoeff.dtype.names,mycoeff.tolist()[0])))
            concentration_table,parameters_index = h5addrow(compilationroot_node,'ksp_fits',row_to_add)
            #}}}
            h5file_node.close()
            ###}}}
        ###{{{ null lists so that I can collect the data
        # do this here, because I might not have any simulated data
        fit_types_simul = []
        datasets_emax_simul = []
        datasets_xisp_simul = []
        description_simul = []
        color_simul = []
        color1_simul = []
        color2_simul = []
        ###}}}
        if fdata_exists: #I also can't generate ernst without t1 info!
            if len(simul_types) > 0:
                ###{{{ Now, simulate too-fast repetition delays
                ####{{{ calculate the scaling
                power_list = emax_unmod.getaxis('power').copy()
                t1p0 = t1f(r_[0.0])
                ernstdelay = 1.5
                TR = ernstdelay * t1p0.data[0]
                TR_clipped = TR.copy()
                if TR_clipped > 0.5:
                    TR_clipped = 0.5
                expvalforernst = t1f(power_list).set_error(None)
                expvalforclipped = expvalforernst.copy()
                ed = exp(-TR/expvalforernst.data.copy())
                ed_clipped = exp(-TR_clipped/expvalforernst.data.copy())
                expvalforernst.data = (1.0-ed)/sqrt(1.0-ed**2)
                expvalforclipped.data = (1.0-ed_clipped)/sqrt(1.0-ed_clipped**2)
                expvalfor5t1 = expvalforernst.copy()
                expvalfor5t1.data = 1.0-(ed**(5.0/ernstdelay))
                del(ed)
                del(ed_clipped)
                ####}}}
                if 'ernst angle' in simul_types:
                    ####{{{ generate ernst
                    fit_types_simul.append('ernst angle simulation')
                    emax_ernst = emax_unmod.copy() * expvalforernst
                    normalizer = emax_ernst.copy().sort('power')['power',0:1]
                    emax_ernst /= normalizer
                    datasets_emax_simul.append(emax_ernst)
                    datasets_xisp_simul.append((1.0 - emax_ernst.copy()) * 1.51671e-3 / f(r_[0.0]))
                    description_simul.append('Ernst Angle\ncalc. for $%0.1f \\times T_1(p=0)$'%ernstdelay)
                    color_simul.append(['#A52A2A','y'])
                    ####}}}
                if 'ernst angle clipped' in simul_types:
                    ####{{{ same code for ernst angle clipped
                    fit_types_simul.append('ernst angle clipped simulation')
                    emax_clipped = emax_unmod.copy() * expvalforclipped
                    normalizer = emax_clipped.copy().sort('power')['power',0:1]
                    emax_clipped /= normalizer
                    datasets_emax_simul.append(emax_clipped)
                    datasets_xisp_simul.append((1.0 - emax_clipped.copy()) * 1.51671e-3 / f(r_[0.0]))
                    description_simul.append('Ernst Angle\ncalc. for $%0.1f \\times T_1(p=0)$ \nunless $T_1>0.5 s$ '%ernstdelay)
                    color_simul.append(['#8800FF','#FF00FF'])
                    ####}}}
                if '5 x T1' in simul_types:
                    ####{{{ same code to generate 5t1
                    fit_types_simul.append('5 x T1 simulation')
                    emax_5t1 = emax_unmod.copy() * expvalfor5t1
                    normalizer = emax_5t1.copy().sort('power')['power',0:1]
                    emax_5t1 /= normalizer
                    datasets_emax_simul.append(emax_5t1)
                    datasets_xisp_simul.append((1.0 - emax_5t1.copy()) * 1.51671e-3 / f(r_[0.0]))
                    description_simul.append('90$^o$ pulse\nwith $5 \\times T_1(p=0)$ delay')
                    color_simul.append(['#00aaaa','#00FFFF'])
                    ####}}}
                if len(color_simul) > 0:
                    color1_simul,color2_simul = map(list,zip(*color_simul))
                ###}}}
        ###{{{ collect into lists
        datasets = [emax_unmod] + datasets_emax_simul
        color1 = [color_pair[0]] + color1_simul
        color2 = [color_pair[1]] + color2_simul
        description = [r'Standard'] + description_simul
        ###}}}
        ###{{{ plot the Emax data
        if verbose: print lsafen('DEBUG: before Emax, figurelist',figurelist)
        ####{{{ set the order of the plots, before plotting to them
        figurelist.append({'print_string':r'\subparagraph{Emax curves}'+'\n\n'})
        nextfigure(figurelist,'Emax')
        show_Evip_plot = False
        if show_Evip_plot:
            nextfigure(figurelist,'Evip')
        figurelist.append({'print_string':r'\Emax\ integration finished'+'\n\n'})
        ####}}}
        if verbose: print lsafen('DEBUG: after Emax, figurelist',figurelist)
        ax = gca()
        for j,v in enumerate(datasets):
            ####{{{ plot the data
            nextfigure(figurelist,'Emax') # switch back to the Emax figure (after linear)
            if ignore_emax_error == True:
                v.set_error(None)
            v = emax(v)
            if guessonly:
                v.settoguess()
            else:
                v.fit()
            w = v.copy()
            if plot_vs_A:
                x = w.getaxis('power')
                x[:] *= w.output('A')
            plot_updown(w,'power',color1[j],color2[j],label = description[j] + lsafe(pdfstring))
            ax = gca()
            myylims = ax.get_ylim()
            textpos = (j+1.0)/(len(datasets)-1.0+2.0)
            w = v.eval(100)
            if plot_vs_A:
                x = w.getaxis('power')
                x[:] *= v.output('A')
                w.rename('power','$p A$')
            plot(w,color1[j])
            ax.set_ylim(myylims)
            if not guessonly:
                #{{{ show the result on the plot
                Emaxerr = v.covar(r'E_{max}')
                if Emaxerr is not None:
                    Emaxerr = sqrt(Emaxerr)
                else:
                    Emaxerr = NaN
                text(0.7,textpos,r'$E_{max}= %0.03f \pm %0.03f$'%(v.output(r'E_{max}'),Emaxerr),
                    transform = ax.transAxes,
                    size = 'x-large',
                    horizontalalignment = 'center',
                    color = color1[j])
                emax_w_error = nddata(v.output(r'E_{max}'),[1],['val']).set_error(r_[Emaxerr])
                ksigma_w_error = (1.0-emax_w_error)/t1_w_error/concentration*1.51671e-3
                krho_w_error = (1.0/t1_w_error - 1.0/t10_w_error)/concentration
                xi_w_error = ksigma_w_error / krho_w_error
                if len(xi_w_error.data) > 1:
                    raise CustomError('length of xi is greater than 1-- what\'s going on!!??')
                tau = interptau(xi_w_error.data[0],sfo1,simple = True)
                onethousandthofxi = xi_w_error.data.copy().flatten()[0]*1e-3
                tau_plusthousandth = interptau(xi_w_error.data[0]+onethousandthofxi,sfo1,simple = True)
                tau_error = abs((tau_plusthousandth-tau)/onethousandthofxi)*(xi_w_error.get_error().copy().flatten()[0])
                report_array = make_rec({'$E_{max}$' : emax_w_error.data.copy().flatten()[0],
                        '$E_{max}$_error' : emax_w_error.get_error().copy().flatten()[0],
                        '$T_1$' : t1_w_error.data.copy().flatten()[0],
                        '$T_1$_error' : t1_w_error.get_error().copy().flatten()[0],
                        '$T_{1,0}$' : t10_w_error.data.copy().flatten()[0],
                        '$T_{1,0}$_error' : t10_w_error.get_error().copy().flatten()[0],
                        r'$k_\rho$' : krho_w_error.data.copy().flatten()[0],
                        r'$k_\rho$_error' : krho_w_error.get_error().copy().flatten()[0],
                        r'$k_\sigma$' : ksigma_w_error.data.copy().flatten()[0],
                        r'$k_\sigma$_error' : ksigma_w_error.get_error().copy().flatten()[0],
                        r'$\xi$' : xi_w_error.data.copy().flatten()[0],
                        r'$\xi$_error' : xi_w_error.get_error().copy().flatten()[0],
                        r'$\tau_c$ / ps' : tau/1e-12,
                        r'$\tau_c$_error / ps' : tau_error/1e-12,
                        })
                report_array = reorder_rec(report_array,['$E_{max}$' , '$E_{max}$_error' , '$T_1$' , '$T_1$_error' , '$T_{1,0}$' , '$T_{1,0}$_error' , r'$k_\rho$' , r'$k_\rho$_error' , r'$k_\sigma$' , r'$k_\sigma$_error' , r'$\xi$' , r'$\xi$_error'])
                figurelist.append({'print_string':r'\subparagraph{Results of uncorrected analysis}'+'\n\n'+lrecordarray(report_array,resizebox = 1.0,return_only=True)})
                figurelist.append({'print_string':r'\begin{tiny}(calculated as follows:) From above, the $T_1$ is '+lsafen(t1_w_error)+r'I find the \Emax to be '+lsafen(emax_w_error)+r'giving \ksm = '+lsafen(ksigma_w_error)+'and sfo1 is %0.4f'%sfo1+r'\end{tiny}'})
                #}}}
                ####{{{ tabulate it
                h5file_node,compilationroot_node = h5nodebypath(h5file + '/compilations')
                mycoeff = v.output()
                search_string = '(Emax_curve == %d)'%(emax_in_file.index)
                mydata = [emax_in_file.index]
                mynames = ['Emax_curve']
                mydata.append(run_number)
                mynames.append('run_number')
                mynames += mycoeff.dtype.names
                mydata += mycoeff.tolist()[0]
                mynames += ['covariance']
                mydata += [v.covarmat()]
                numrem,junk = h5remrows(compilationroot_node,'Emax_fits',search_string)
                if numrem:
                    print lsafen('Removed',numrem,'rows from Emax_fits, with values',junk[r'E_{max}'],'because they conflict with these entries')
                concentration_table,parameters_index = h5addrow(compilationroot_node,'Emax_fits',mydata,mynames,match_row = search_string)
                h5file_node.close()
                ####}}}
            ####}}}
            plot_newlinear = False
            if plot_newlinear:
                ####{{{ new linear plot
                figurelist.append({'legend':True})
                nextfigure(figurelist,'newlinear')
                thisemax = v.output(r'E_{max}')
                print ndshape(v)
                newlinear = (1.-thisemax)/(v.copy()-thisemax)
                current_color = plot_color_counter()
                plot(newlinear,'o',label = description[j])
                plot_color_counter(current_color)
                maxp = lambda x: x == x.max()
                print 'shape of newlinear:',ndshape(newlinear)
                plot(r_[0.0,newlinear.getaxis('power').max()],r_[1.0,newlinear['power',maxp].data[0]],'-',alpha = 0.5,label = description[j])
                xlabel('power / $W$')
                ylabel(r'$\frac{1-E_{max}}{E-E_{max}}$')
                ####}}}
            figurelist = emax_linearandplot(v,first_figure = figurelist,max_invpower = 100, color_pair = color_pair[::-1],show_Evip_plot = show_Evip_plot)
        #gridandtick(gca())
        #autolegend()
        ###}}}
        if fdata_exists:
            process_xisp = False
            if process_xisp:
                ###{{{ process the xi s(p) data
                ####{{{ collect into lists
                datasets = [xisp_f,xisp_fp] + datasets_xisp_simul
                color1 = ['b','k'] + color1_simul
                color2 = ['g','r'] + color2_simul
                description = [r'$(1-E(p))/f(0)$',r'$(1-E(p))/f(p)$'] + description_simul
                fit_type = ['constant f','power dependent f'] + fit_types_simul
                ####}}}
                nextfigure(figurelist,'multifitxismax' + pdfstring)
                ax = gca()
                h5file_node,compilationroot_node = h5nodebypath(h5file + '/compilations')
                if not dontfit:
                    for j,v in enumerate(datasets):
                        ####{{{ plot the data
                        plot_updown(v,'power',color1[j],color2[j],label = description[j])
                        v = xismax(v)
                        try:
                            v.fit()
                        except:
                            raise CustomError('Try setting dontfit to True; it is currently',dontfit)
                        textpos = (j+1.0)/(len(datasets)-1.0+2.0)
                        plot(v.eval(100),color1[j],label = description[j])
                        text(0.7,textpos,r'$\xi s_{max}= %0.03f \pm %0.03f$'%(v.output(r'\xi s_{max}'),sqrt(v.covar(r'\xi s_{max}'))),
                                transform = ax.transAxes,
                                size = 'x-large',
                                horizontalalignment = 'center',
                                color = color1[j])
                        ####}}}
                        ####{{{ tabulate it
                        mycoeff = v.output()
                        search_string = '(Emax_curve == %d) & (fit_type == \'%s\')'%(emax_in_file.index,fit_type[j])
                        mydata = [emax_in_file.index,fit_type[j]]
                        mynames = ['Emax_curve','fit_type']
                        mydata.append(run_number)
                        mynames.append('run_number')
                        mynames += mycoeff.dtype.names
                        mydata += mycoeff.tolist()[0]
                        mynames += ['covariance']
                        mydata += [v.covarmat()]
                        numrem,junk = h5remrows(compilationroot_node,'xismax_fits',search_string)
                        if numrem:
                            print lsafen('Removed',numrem,'rows from xismax_fits, with values',junk[r'\xi s_{max}'],'because they conflict with these entries')
                        concentration_table,parameters_index = h5addrow(compilationroot_node,'xismax_fits',mydata,mynames,match_row = search_string)
                        ####}}}
                h5file_node.close()
                autolegend()
                title(r'$\xi s(p)$ by various methods')
                ###}}}
    return figlistret(first_figure,figurelist,basename = fid_args['name'])
def T1_interpolation(t1dataset,
        t10_at_zero_power = 2.5,
        T1w = 2.4691,
        DeltaT1w = 0.10267,
        t1_bad_powers = None,
        fl = None):
        if fl is None:
            fl = figlist_var()
        #{{{ Set all the T10 related values
        T1w_inv_f = lambda x: nddata(1.0/(T1w+x*DeltaT1w),[-1],['power']).labels('power',x)
        print '\n\nlinear interpolation using $T_{1,0}(0)$ of %0.3f for $k_\\sigma$ analysis, which is'%t10_at_zero_power
        Ck_HH = 1./t10_at_zero_power-1./T1w
        print '\n\nThis implies that $C_{macro}k_{HH}$ is %0.3f\n\n'%(Ck_HH)
        #}}}
        #{{{ T_1 interpolation
        #{{{ if there are any scans that are bad (i.e. close to the signal null, remove them by power here
        if t1_bad_powers is not None:
            mask = zeros(size(t1dataset.getaxis('power')),dtype = 'bool8')
            for j in range(0,len(t1_bad_powers)):
                mask |= abs(t1dataset.getaxis('power') - t1_bad_powers[j]) < abs(0.05*t1_bad_powers[j])
                figurelist.append({'print_string':
                    r"\textbf{I'm leaving out $T_1$'s with power %g, because you marked them as bad (there are %d of these)}"%(t1_bad_powers[j],sum(mask))})
            new_t1dataset = nddata(t1dataset.data[~mask],[sum(~mask)],['power'])
            new_t1dataset.labels('power',t1dataset.getaxis('power')[~mask])
            new_t1dataset.set_error(t1dataset.get_error()[~mask])
            t1dataset = new_t1dataset
        #}}}
        t1_at_zero_power = t1dataset['power',0.0].data[0]
        print '\n\nlinear interpolation using $T_1(0)$ of %0.3f for $k_\\sigma$ analysis\n\n'%t1_at_zero_power
        rho = (1./t1_at_zero_power) - (1./t10_at_zero_power)
        print '\n\nlinear interpolation using $\\rho$ of %0.5f for $k_\\sigma$ analysis\n\n'%rho
        #{{{ find list of independent powers
        list_of_t1_powers = t1dataset.getaxis('power').copy()
        list_of_different_t1_powers = array([])
        for j in range(0,len(list_of_t1_powers)):
            if not any(abs(list_of_different_t1_powers -
                    list_of_t1_powers[j])
                    < 0.05*list_of_t1_powers[j]):
                list_of_different_t1_powers = r_[
                        list_of_different_t1_powers,
                        list_of_t1_powers[j]]
        fl.text(r'I found %d different $T_1$ powers:'%len(list_of_different_t1_powers)
            +lsafen(list_of_different_t1_powers))
        #}}}
        F_nonlinear = 1./(1./t1dataset.copy().set_error(None)
                - Ck_HH
                - T1w_inv_f(t1dataset.getaxis('power')))
        #{{{ threshold, in case the inverse was too small
        threshold_number = 1./0.5e-3 # note that this is 2000s!
        def threshold_F_nonlinear(myx):
            mask = myx.data > threshold_number # if the inverse was very small
            mask |= myx.data <= 0.0 # if the inverse was negative
            mask |= isnan(myx.data) # if it was exactly zero
            myx.data[mask] = threshold_number
        threshold_F_nonlinear(F_nonlinear)
        #}}}
        if len(list_of_different_t1_powers) > 3:
            F_nonlinear_to_order = 2
        else:
            F_nonlinear_to_order = 1
        c_F_nonlinear,F_nonlinear_line = F_nonlinear.polyfit('power',order = F_nonlinear_to_order)
        print '\n\n$F\\_{nonlinear}$ is',lsafen(F_nonlinear),'and its fit is',F_nonlinear_line
        # now, linearly interpolate w.r.t. power
        # in the same step, convert back to T1(p)
        t1err = t1dataset['power',0:1]
        t1err /= t1err.data[0] # so this way, multiplying by this will just give an nddata with teh correct error.
        DeltaT10 = c_F_nonlinear[1]
        def F_nonlinear_f(x):
            x_data = nddata(x, len(x), ['power'])
            #{{{ evaluate c_F_nonlinear over the nonlinear polynomial
            retval = c_F_nonlinear[0]
            for j in range(1,F_nonlinear_to_order+1):
                retval += x_data**j*c_F_nonlinear[j]
            #}}}
            retval.labels('power',x)
            threshold_F_nonlinear(retval)
            return retval
        t1f = lambda x: t1err/(1./F_nonlinear_f(x) + Ck_HH + T1w_inv_f(x))
        powers_forplot = linspace(0,t1dataset.getaxis('power').max(),100)
        #{{{ just diagnosis
        #{{{ plot F_nonlinear directly
        #{{{ plot stuff
        fl.next('Fnonlinear')
        fl.plot(F_nonlinear,'o',label = '$F_{nonlinear}$ points')
        #}}}
        powers_fordiagnose = linspace(0,t1dataset.getaxis('power').max(),100)
        #print '\n\n$T_1(p)$ is',lsafen(t1f(powers_fordiagnose))
        closely_spaced_F = F_nonlinear_f(powers_fordiagnose)
        fl.plot(F_nonlinear,'-',alpha=0.8,label = '$F_{nonlinear}(p)$')
        fl.plot(F_nonlinear_line,'-',alpha=0.8,label = r'$F_{nonlinear,origfit}(p)$')
        fl.plot(closely_spaced_F,'-',alpha=0.8,label = '$F_{nonlinear,interpolated}(p)$')
        #title(name)
        autolegend()
        #}}}
        #{{{ plot inverse of F_nonlinear
        fl.next('invFnonlinear')
        fl.plot(T1w_inv_f(t1dataset.getaxis('power'))+rho,
                '-',label = '$T_{1,w}^{-1}(p)+\\rho(0)$')
        fl.plot(1.0/F_nonlinear,'o',alpha=0.8,label = '$F_{nonlinear}^{-1}(p)$')
        fl.plot(1.0/F_nonlinear_line,'-',alpha=0.8,label = '$F_{nonlinear}^{-1}(p)$')
        fl.plot(1./t1dataset.copy().set_error(None),
                'o', label = '$T_1^{-1}$')
        fl.plot(powers_forplot,1.0/t1f(powers_forplot).set_error(None),
                '-', alpha = 0.2, label = '$T_1^{-1}(p) = F_{nonlinear}^{-1}(p)+(T_{1,w}^-1+p \\Delta T_{1,w}^{-1})^{-1} + C_{macro} k_{HH}$')
        myones = F_nonlinear.copy()
        myones.data[:] = 1.0
        fl.plot(rho*myones,
                '--',
                label = '$\\rho(0)$')
        fl.plot((1./t10_at_zero_power+rho)*myones,
                '--',
                label = '$T_{1,0}^{-1} + \\rho(0)$')
        expand_y()
        expand_x()
        gridandtick(gca())
        #title(name)
        autolegend()
        #}}}
        #}}}
        DeltaT10 = c_F_nonlinear[1]
        ###}}}
        return t1f,powers_forplot
##}}}
#}}}
#{{{ retrieve the t1 power series results from an HDF file
def plot_t1series(path,name,t1expnos,dbm,h5file = 'temperature_paper.h5',gensvg = False,**kwargs):
    dnp_for_rho(path,name,[],expno = [],t1expnos = t1expnos,t1powers = dbm,**kwargs)
    h5file = openFile(h5file,'r+')
    T1_table = h5file.root.concentration_series.T1
    print r'\begin{verbatim}'
    search_string = '(chemical_name == "%s")'%pull_chemical_from_filename(name)
    search_string += ' & '+gensearch('concentration','%0.8g',pull_concentration_from_filename(name),1e-6)
    search_string += ' & '+gensearch('date','%0.3f',pull_date_from_filename(name),0.1)
    results = T1_table.readWhere(search_string)
    print results
    print r'\end{verbatim}'
    figure(100)
    clf()
    t1data = nddata(results['value'],[len(results),1],['power',r'$T_1$'])
    #print 'powers are ',results['power']
    t1data.labels(['power'],[dbm_to_power(results['power'])])
    t1data.sort('power')
    c,straightline = t1data.polyfit('power')
    print r'found coeff \begin{verbatim}',c,r'\end{verbatim}'
    plot(straightline,'-')
    plot(t1data,'o')
    lplot('T1vp'+name+'_%d_%d.pdf'%(t1expnos[0],t1expnos[-1]),gensvg = gensvg)
    figure(101)
    clf()
    m0data = nddata(results['Minf'],[len(results),1],['power',r'$M(\infty)$'])
    m0data.labels(['power'],[dbm_to_power(results['power'])])
    m0data.sort('power')
    m0data /= m0data['power',0]
    plot(m0data,'o')
    lplot('T1vpMinf'+name+'_%d_%d.pdf'%(t1expnos[0],t1expnos[-1]))
    c_forstore = zeros(10)
    c_forstore[0:size(c)] = c.flatten()
    foundrow = False
    T1fit_table = h5file.root.concentration_series.T1fit
    for row in T1fit_table.where(search_string):
        foundrow = True
        row['chemical_name'] = pull_chemical_from_filename(name)
        row['concentration'] = pull_concentration_from_filename(name)
        row['date'] = pull_date_from_filename(name)
        row['c_len'] = len(c)
        row['c'] = c_forstore
        row.update()
    if not foundrow:
        row = T1fit_table.row
        row['chemical_name'] = pull_chemical_from_filename(name)
        row['concentration'] = pull_concentration_from_filename(name)
        row['date'] = pull_date_from_filename(name)
        row['c_len'] = len(c)
        row['c'] = c_forstore
        row.append()
    h5file.close()
    return m0data,c,straightline
#}}}
#{{{ stuff to deal with T1 vs. p data:
#t1vp: load/retrieve T1 vs. p series from HDF5 file --> new version of previous
#t10data: retrieve T10 data for rigorous Emax processing
##{{{ parse the covariance matrix in the HDF5 file
def parse_t1_covariance(data):
    if len(data) == 0:
        raise CustomError('length of data is 0!')
    num_coeff = data['num_coeff'][0]
    covariance = zeros((num_coeff**2,len(data)))
    for j in range(0,len(data)):
        for k in range(0,num_coeff**2):
            try:
                covariance[k,j] = data['covarelement%d'%k][j]
            except:
                raise CustomError("Couldn't find the field I was looking for in",data.dtype.names)
    return covariance.reshape((num_coeff,num_coeff,len(data)))
##}}}
##{{{ load a series of t1 vs power
def t1vp(h5filename,expnos,dbm,fid_args = {},integration_args = {}, auxiliary_args = {},chem_conc = None,chem_name = None,first_figure = None):
    r'''this tries to load a series of T1 vs. power curves
    it uses expno and t1power to select the T1 files to load,
    and it calls store T1'''
    figurelist = figlistini_old(first_figure)
    if len(expnos) != len(dbm):
        raise CustomError('len(expnos)=',len(expnos),'len(dbm)=',len(dbm),dbm)
    if chem_name is None:
        if chem_conc is not None:
            raise CustomError('either BOTH chemical name and concentration must not be passed, or both must!!!')
        chem_name = pull_chemical_from_filename(fid_args['name'])
        chem_conc = pull_concentration_from_filename(fid_args['name'])
    chemidx = get_chemical_index(h5filename+'/compilations/chemicals',chem_name,chem_conc)
    figurelist.append({'print_string':r'\subparagraph{$T_1$ raw data}'+'\n\n'})
    for j in range(0,len(expnos)):
        fid_args.update({'expno':expnos[j]})
        auxiliary_args.update({'t1power':dbm[j]})
        store_T1(h5filename,chemidx,
                fid_args,
                first_figure = figurelist,
                auxiliary_args = auxiliary_args,
                integration_args = integration_args)
    h5file = tables.openFile(h5filename)
    data = h5file.root.compilations.T1_fits.readWhere('chemical_id == %d'%chemidx)
    h5file.close()
    covariance = parse_t1_covariance(data)
    errorbars = zeros(len(data))
    for j in range(0,len(data)):
        #errorbars[j] = sqrt(1.0/(pinv(covariance[:,:,j])[0,0]))
        errorbars[j] = sqrt(covariance[:,:,j][2,2])
    powers = data['power'][:]
    data = data['T_1'][:]
    data = nddata(data,[len(data)],['power'],data_error = errorbars,axis_coords = [dbm_to_power(powers)])
    print data
    plot(data,'o')
    lplot('T1vp_%s.pdf'%fid_args['name'])
##}}}
#{{{ compile all the DNP data for a given chemical, etc
def retrieve_DNP_set(chemical_list,h5file = 'dnp.h5',fit_type = 'corrected',t10_additional_search = '',verbose = False,divide_klow_by = None,t10subst = None):
    if t10subst is not None:
        obsn(r"You have chosen to make the following \tonen substitutions")
        temp = empty(len(t10subst),dtype = [('dataset','|S100'),('uses $T_{1,0}$ from','|S100')])
        temp['dataset'] = array(t10subst.keys())
        temp['uses $T_{1,0}$ from'] = array(t10subst.values())
        lrecordarray(temp,std_sf = 7,scientific_notation = False)
        print '\n\n'
    h5file = tables.openFile(h5file) # open the HDF5 file
    #{{{ retrieve chemical information
    search_string = h5inlist('chemical',chemical_list)
    chemical_data = h5file.root.compilations.chemicals.chemicals.readWhere(search_string)
    if verbose:
        print "All the relevant chemicals (chemical\\_data) I see are:\n\n"
        lrecordarray(chemical_data,format = '%g')
        print '\n\n'
    #}}}
    #{{{ join ksp data onto the chemical information
    data = h5join((h5file.root.compilations.Emax_curves,'chemical_id'),
        (chemical_data,'index'),
        pop_fields = ['experiments'])# the first time, run it verbose, so we can see how this works
    data = h5join((h5file.root.compilations.ksp_fits,'Emax_curve'),
        (data,'index'),
        additional_search = 'fit_type == "%s"'%fit_type)
    if verbose:
        print "All \\ksp data (data) I see are:\n\n"
        lrecordarray(data,resizebox = True,std_sf = 7,scientific_notation = False)
        print '\n\n'
    #}}}
    #{{{ pull the T10 data
    list_of_interesting_fields = ['chemical','run_number']
    mask = chemical_data['concentration'] == 0. # to select just the T1,0 samples
    if verbose:
        print "The indeces for the \\tonen samples are:\n\n"
        lrecordarray(chemical_data[mask],format = '%g')
        print '\n\n'
    if t10subst is not None:
        temp = data['chemical'].copy()
        data['chemical'] = [t10subst[j] if j in t10subst.keys() else j for j in data['chemical']]
        t10data = decorate_rec( (data[list_of_interesting_fields],'chemical'),
            (chemical_data[mask],'chemical'))# decorate those interesting fields with the T10 chemical id
        t10data['chemical'] = temp
        data['chemical'] = temp
    else:
        t10data = decorate_rec( (data[list_of_interesting_fields],'chemical'),
            (chemical_data[mask],'chemical'))# decorate those interesting fields with the T10 chemical id
    t10data = rename_fields(t10data,{'index':'chemical_id'})# numpy function to rename the fields
    if verbose:
        print "selecting the fields",lsafe(list_of_interesting_fields),r"from the \ksp data (\texttt{data})"
        print "and decorating with the chemical ID of the \\tonen samples gives the new t10data:\n\n"
        lrecordarray(t10data,format = '%g')
        print '\n\n'
    list_of_interesting_fields.append('T_1')
    #{{{ here, I specify only by chemical id, so I can pull the data later
    if len(t10_additional_search) > 0:
        temp = ' & '.join(map(lambda x: '('+x+')',['power == -999.',t10_additional_search]))
    else:
        temp = 'power == -999.'
    t10data = h5join((h5file.root.compilations.T1_fits,['chemical_id']),
        (t10data,['chemical_id']),
        additional_search = temp,
        select_fields = list_of_interesting_fields)
    #}}}
    # this is where I could go back and find stuff that matches the chemical ID but not the run number
    if len(t10data) == 0:
        raise ValueError("I didn't find any T10 data!")
    t10data = applyto_rec(mean,t10data,['run_number','chemical'])
    t10data = rename_fields(t10data,{'T_1':'T_{1,0}'})
    if verbose:
        print r"After joining the \tonen data onto the t10data array and averaging, it becomes:"
        print '\n\n'
        lrecordarray(t10data,std_sf = 7,scientific_notation = False)
        print '\n\n'
    #}}}
    #{{{ reorder the data, and specify this is the chemical ID for the spin labeled compound
    data = data[['run_number','chemical','fit_type','ksmax','chemical_id','concentration']]
    data = rename_fields(data,{'chemical_id':'SL_chemical_id'})
    data.sort()
    #}}}
    #{{{ add the T10 data onto the data
    newdata,missing_t10s = decorate_rec((data,['chemical','run_number']),
        (t10data,['chemical','run_number']),
        drop_rows = 'return')
    if verbose:
        print r"Decorating the \ksp data with \texttt{t10data} gives \texttt{newdata} where there was a successfully identified \tonen:"
        print '\n\n'
        lrecordarray(newdata,std_sf = 7,scientific_notation = False)
        print '\n\n'
        print r"and the dropped rows (\texttt{missing_t10s}):"
        print '\n\n'
        lrecordarray(missing_t10s,std_sf = 7,scientific_notation = False)
        print '\n\n'
    if missing_t10s is not None:
        #{{{ alert to, and take care of the missing t10 data
        obs('\n\n$T_{1,0}$ experiments are missing for the following chemical / run number pairs:\n\n')
        print r'{\color{red}'
        lrecordarray(missing_t10s[['chemical','run_number']],std_sf = 7,scientific_notation = False)
        t10fields = list(t10data.dtype.names)
        print '\n\n',r'To rectify this, I make the following assignments','\n\n'
        missing_t10s = rename_fields(missing_t10s,{'run_number':'DNP_run_number'})
        missing_t10s = decorate_rec((missing_t10s,'chemical'),
            (t10data,'chemical'),
            drop_rows = True)
        print r'}'
        missing_t10s = rename_fields(missing_t10s,{'run_number':'T10_run_number'})
        lrecordarray(missing_t10s,resizebox = True,std_sf = 7,scientific_notation = False)
        #{{{ put into the write format and average them
        print r'Then average for all degenerate entries of chemical and DNP run number:','\n\n'
        missing_t10s = rename_fields(missing_t10s,{'DNP_run_number':'run_number'})
        missing_t10s = missing_t10s[t10fields] # go back to fields used in the t10data function, because otherwise I have a string and it freaks out, and also because I need to put them together, below
        missing_t10s = reorder_rec(missing_t10s,t10fields)# because apparently, the above does not accomplish a reorder
        applyto_rec(mean,missing_t10s,['chemical','run_number'])
        lrecordarray(missing_t10s,std_sf = 7,scientific_notation = False)
        print '\n\n'
        #}}}
        #{{{ Finally, add the second-choice T10's in, and retrieve the resulting data
        try:
            t10data = r_[t10data,missing_t10s]
        except:
            raise CustomError('There is a problem with either t10data=',t10data,'or missing t10s=',missing_t10s)
        # try again
        newdata,missing_t10s = decorate_rec((data,['chemical','run_number']),
            (t10data,['chemical','run_number']),
            drop_rows = 'return')
        if missing_t10s is not None:
            obsn('$T_{1,0}$ experiments are STILL missing for the following chemical / run number pairs!!')
            print r'{\color{red}'
            lrecordarray(missing_t10s[['chemical','run_number']],std_sf = 7,scientific_notation = False)
            print r'}'
        #}}}
        #}}}
        print r"{\color{red} $\Rightarrow$ done substituting $T_{1,0}$'s}",'\n\n'
    data = newdata
    #}}}
    #{{{ Tack on the T1 data.
    # Here we note that we are only interested in the fields I have right now
    list_of_interesting_fields = list(data.dtype.names)
    # and the T_1
    list_of_interesting_fields.append('T_1')
    list_of_interesting_fields.remove('SL_chemical_id')
    data = h5join((h5file.root.compilations.T1_fits,['run_number','chemical_id']),
        (data,['run_number','SL_chemical_id']),
        additional_search = 'power == -999.',
        select_fields = list_of_interesting_fields)
    #}}}
    h5file.close() # close the HDF5 file, since all the rest is just spreadsheet math
    if verbose:
        print r"After joining the data for the spin-labeled $T_1$'s onto \texttt{data}, we have","\n\n"
        lrecordarray(data,resizebox = 0.8,std_sf = 7,scientific_notation = False)
        print '\n\n'
    if verbose: print r"Now, I use lambda\_rec~\ref{codelabel:lambda_rec} to calculate $k_\rho$, then $k_{low}$:"+'\n\n'
    #{{{ take the mean of the T10 values, so I only get one krho for data_nice
    nont10fields = list(data.dtype.names)
    try:
        temp = nont10fields.index('T_{1,0}')
    except:
        raise CustomError('list is',nont10fields)
    nont10fields.pop(temp)
    data_nice = applyto_rec(mean,data,nont10fields)
    #}}}
    #{{{ calculate krho, then klow; do for both tables, since I get many krho for one and only one for the other
    data = lambda_rec(data,
        r'krho',
        (lambda x,y,C: (1.0/x - 1.0/y)/C),
        ['T_1','T_{1,0}','concentration'])
    data_nice = lambda_rec(data_nice,
        r'krho',
        (lambda x,y,C: (1.0/x - 1.0/y)/C),
        ['T_1','T_{1,0}','concentration'])
    # 5/3rho-7/3sigma
    data = lambda_rec(data,
        r'klow',
        (lambda x,y: 5./3.*x - 7./3.*y),
        ['krho','ksmax'])
    data_nice = lambda_rec(data_nice,
        r'klow',
        (lambda x,y: 5./3.*x - 7./3.*y),
        ['krho','ksmax'])
    data_nice = lambda_rec(data_nice,
        r'\xi',
        (lambda x,y: x/y),
        ['ksmax','krho'])
    #}}}
    #{{{ reference against free spin label values, and give it nicer names
    data_nice = data_nice[['chemical','run_number','ksmax','klow',r'\xi']]
    data_nice = reorder_rec(data_nice,['chemical','run_number','ksmax','klow',r'\xi'])
    data_nice = lambda_rec(data_nice,
        r'ksmax',
        (lambda x: x/bulk_ksigma_value),
        ['ksmax'])
    extra_factor = 1
    if divide_klow_by is not None:
        extra_factor = divide_klow_by
    data_nice = lambda_rec(data_nice,
        r'klow',
        (lambda x: x/bulk_klow_value/extra_factor),
        ['klow'])
    data_nice = lambda_rec(data_nice,
        r'\xi',
        (lambda x: x/bulk_xi_value),
        [r'\xi'])
    #{{{ use clumsier but nicer names for my table
    if divide_klow_by is None:
        bulk_klow_name_new = bulk_klow_name
    else:
        bulk_klow_name_new = '%s/%d'%(bulk_klow_name,extra_factor)
    xi_name = r'\xi/\xi_{bulk}'
    data_nice = rename_fields(data_nice,{r'\xi':xi_name,
                                'ksmax':bulk_ksmax_name,
                                'klow':bulk_klow_name_new})
    data_nice = reorder_rec(data_nice,['chemical','run_number',bulk_klow_name_new,bulk_ksmax_name,xi_name])
    data_nice.sort()
    #}}}
    #}}}
    return data,data_nice
#}}}
##{{{ retrieve a series of T1 measurements from HDF5
def retrieve_T1series(h5filename,name,*cheminfo,**kwargs):
    'cheminfo -- either file name, or empty to use name'
    #{{{ kwargs processing
    verbose = False
    indirect_dim = 'power'
    run_number = None
    show_result = False
    retcheckval = False
    if 'retcheckval' in kwargs.keys():
        retcheckval = kwargs.pop('retcheckval')
    if 'indirect_dim' in kwargs.keys():
        indirect_dim = kwargs.pop('indirect_dim')
    if 'verbose' in kwargs.keys():
        verbose = kwargs.pop('verbose')
        if verbose == True and show_result == False:
            show_result = True
    if 'run_number' in kwargs.keys():
        run_number = double(kwargs.pop('run_number'))
    if 'show_result' in kwargs.keys():
        show_result = kwargs.pop('show_result')
    t1max = None
    if 't1max' in kwargs.keys():
        t1max = kwargs.pop('t1max')
    d1val = None
    if 'd1val' in kwargs.keys():
        d1val = kwargs.pop('d1val')
    if verbose: print lsafen('DEBUG: retrieve T1series called with name=',name,'cheminfo',cheminfo)
    #}}}
    if len(cheminfo) == 2:
        chemical_name = cheminfo[0]
        chemical_concentration = cheminfo[1]
    elif len(cheminfo) == 0:
        chemical_name = pull_chemical_from_filename(name)
        chemical_concentration = pull_concentration_from_filename(name)
    else:
        raise CustomError('You called retrieve_T1series with neither no info about the chemical nor a chemical, concentration pair')
    chemidx = get_chemical_index(h5filename+'/compilations/chemicals',chemical_name,chemical_concentration)
    h5file = tables.openFile(h5filename)
    search_string = 'chemical_id == %d'%chemidx
    if run_number is not None:
        search_string = '('+search_string+') & (run_number == %0.3f)'%run_number
    if verbose: print "search string is",lsafen(search_string)
    data = h5file.root.compilations.T1_fits.readWhere(search_string)
    if verbose: print "search data is",lsafen(data)
    h5file.close()
    if name is not None: # only run this code if I pass an explicit name
        data = [x.reshape(1) for x in data if re.match('^'+name,x['integrals'])]
        if len(data) == 0:
            if show_result:
                raise RuntimeError("found no T1 data "+show_result+" for chemical "+name)
            else:
                raise RuntimeError("found no T1 data for"+name)
        data = concatenate(data) # this concatenates the nddata
    if len(data) == 0:
        print 'No $T_1$ data found for chemidx=\n\n'
        print_chemical_by_index(h5filename+'/compilations/chemicals',chemidx)
        print '\n\n'
        return None
    else:
        print '$T_1$ data was found for chemidx=\n\n'
        print_chemical_by_index(h5filename+'/compilations/chemicals',chemidx)
        print '\n\n'
    if show_result:
        print r'\subparagraph{$T_1$ series}','\n\n',r'{\bf retrieved $T_1$ series',show_result,':}','\\begin{tiny}\n\n',lrecordarray(data),'\\end{tiny}\n\n'
        temp = [x[0:x.find('_exp')] for x in data['integrals'].tolist() if x.find('_exp')>-1]
        temp = unique(array(temp))
        if len(temp)>1:
            print r"{\Large{\color{red} Warning!!}} You are using more than one dataset (",lsafe(temp),") for your $T_1$ experiments; unless you are doing this on purpose, it's a mistake.  You probably need to change your \\texttt{run\_number}, which should be the same for all samples of a \\textit{different} chemical composition in a run, but should \\textit{not} be the same for repeats with the same chemical composition (that would be a different run...)\n\n"
        else:
            print "$T_1$ datasets unique\n\n"
    if t1max is not None:
        actual_t1max = data['T_1'].max()
        print r'You set $T_{1,max}$ to %0.3f s'%t1max,"($\\frac{d1}{5}=%0.3f$ s) which compared to %0.2f"%(d1val/5.,actual_t1max),
        if actual_t1max>t1max:
            print r'{\Large{\color{red} is NOT OK!!!!}}'
        else:
            print r'{\large is OK}'
    if verbose: print '\n\n',lsafe('matching indeces:',data['index']),'\n\n'
    myerrors = sqrt(parse_t1_covariance(data)[2,2,:])# the T1 is the third thing, not the first!!
    if indirect_dim == 'power':
        powers = data[indirect_dim][:]
        powers = dbm_to_power(powers)
    elif indirect_dim in data.dtype.names:
            powers = data[indirect_dim][:]
    else:
            raise CustomError('indirect dim',indirect_dim,'not in',data.dtype.names)
    Minfdata = nddata(data[r'M(\infty)'],[len(data)],[indirect_dim],data_error = myerrors,axis_coords = [powers])		
    M0data = nddata(data[r'M(0)'],[len(data)],[indirect_dim],data_error = myerrors,axis_coords = [powers])		
    data = data['T_1'][:]
    retval = nddata(data,[len(data)],[indirect_dim],data_error = myerrors,axis_coords = [powers])
    if indirect_dim == 'power':
        retval.sort(indirect_dim)
    retval.name('T_1')
    if retcheckval:
	    return retval,M0data,Minfdata
    else:
        return retval
##}}}
#}}}
#{{{ standard EPR processing
def standard_epr(dir = None,
        files = None,
        linewidth = 1.8,
        normalize_field = False,
        find_maxslope = False,
        maxslope_within = [0,0.5],# a fractional amount the maxslope is within
        offset_spectra = False,
        subtract_first = False,
        normalize_peak = False,
        background = None,
        grid = False,
        h5file = 'dnp.h5',
        run_number = None,
        figure_list = None,#the next (conc_calib) gives the concentration in mM/(signal units / mm)
        conc_calib = None):
    """the files argument is either:
    a single string for the file
    a list of files
    a list of tuples with sample length (to calculate a
        "concentration") and/or a chemical_name, like this:
        (filename,sample_length,chemical_name)"""
    if type(run_number) in [int,long,int32,int64]:
        run_number = double(run_number)
    if type(files) is str or type(files) is tuple:
        files = [files]
    if background is not None:
        if type(background) is list:
            if len(background) > 1:
                raise CustomError("multiple backgrounds not yet supported")
            else:
                background = background[0]
        files = [background] + files
        subtract_first = True
    dir = dirformat(dir)
    mu_B = 9.27400915e-24
    h = 6.62606957e-34
    if subtract_first:
        firstdata = 2*load_indiv_file(dir+files.pop(0))
    for index,file in enumerate(files):
        if type(file) is tuple:
            if len(file) > 2:
                chemical_name = file[2]
            else:
                chemical_name = file[0]
            sample_len = file[1]
            file = file[0]
        else:
            sample_len = None
            chemical_name = file
        try:
            data = load_indiv_file(dir+file)
        except:
            raise CustomError('Error loading file: type of dir',type(dir),'type of file',type(file),'file=',file)
        if subtract_first:
            data -= firstdata
        field = r'$B_0$'
        neworder = list(data.dimlabels)
        data.reorder([neworder.pop(neworder.index(field))]+neworder) # in case this is a saturation experiment
        data -= data.copy().run_nopop(mean,field)
        figure_list.next('epr',legend = True,boundaries = False)
        v = winepr_load_acqu(dir+file)
        if index == 0:
             fieldbar = data[field,lambda x: logical_and(x>x.mean(),x<x.mean()+10.)]
             fieldbar.data[:] = 0.5
             fieldbar.data[0] = 0.6
             fieldbar.data[-1] = 0.6
             fxaxis = fieldbar.getaxis(field)
        xaxis = data.getaxis(field)
        centerfield = None
        if normalize_field:
            xaxis /= v['MF']
            if index == 0:
                fxaxis /= v['MF']
            newname = r'$B_0/\nu_e$'
        else:
            deriv = data.copy()
            deriv.run_nopop(diff,field)
            deriv.data[abs(data.data) > abs(data.data).max()/10.] = 0 # so it doesn't give a fast slope, non-zero area
            deriv = abs(deriv)
            thisstart = deriv.getaxis(field)[0]
            thiswidth = deriv.getaxis(field)[-1]-thisstart
            deriv = deriv[field,lambda x: logical_and(
                (x-thisstart)>thiswidth * maxslope_within[0]
                ,
                (x-thisstart)<thiswidth * maxslope_within[1]
                )]
            deriv.argmax(field,raw_index = True)
            centerfield = mean(xaxis[int32(deriv.data)]) # "centerfield" is the center of peak, where the derivative is maximized
            if find_maxslope:
                xaxis -= centerfield
                if index == 0:
                    fxaxis -= centerfield
                newname = r'$\Delta B_0$'
            else:
                newname = field
        data.rename(field,newname)
        if index == 0:
            fieldbar.rename(field,newname)
        mask = data.getaxis(newname)
        mask = mask > mask[int32(len(mask)-len(mask)/20)]
        snr = abs(data.data).max()/std(data.data[mask])
        integral = data.copy()
        integral.data -= integral.data.mean() # baseline correct it
        integral.integrate(newname)
        figure_list.next('epr_int',legend = True,boundaries = False)
        pc = plot_color_counter()
        plot_color_counter(pc)
        legendtext = chemical_name
        legendtext += '\n'
        if centerfield != None:
            legendtext += ', %0.03f $G$'%centerfield
            gfactor = h*v['MF']*1e9/mu_B/(centerfield*1e-4)
            legendtext += ', g=%0.03f'%gfactor
        doubleintegral = integral.copy()
        doubleintegral.integrate(newname)
        if sample_len is None:
            if conc_calib is not None:
                raise ValueError("If you want to set the concentration calibration, you need to give a sample length!")
            legendtext += r', SNR %0.2g $\int\int=$ %0.3g sig. units'%(snr,doubleintegral[newname,-1].data[-1])
        else:
            if conc_calib is None:
                legendtext += r', SNR %0.2g $\int\int=$ %0.3g sig. units / mm'%(snr,doubleintegral[newname,-1].data[-1]/double(sample_len))
            else:
                thisconc = doubleintegral[newname,-1].data[-1]/double(sample_len)*double(conc_calib)
                legendtext += r', SNR %0.2g $\int\int=$ %0.3g mM'%(snr,thisconc)
                if run_number is not None:
                    chem_id = get_chemical_index(h5file+'/compilations/chemicals',chemical_name,thisconc*1e-3)
                    h5fp,thisnode = h5nodebypath(h5file+'/compilations/chemicals')
                    h5remrows(thisnode,
                            'default_chemical_id',
                            "(run_number == %0.3f) & (chemical == '%s')"%(run_number,chemical_name))# clear any data for this run number and chemical name that was previously entered (since there can be only one entry)
                    h5addrow(thisnode,
                            'default_chemical_id',
                            {'run_number':run_number,'chemical':chemical_name,'chemical_id':chem_id},
                            match_row = "(chemical_id == %d) & (run_number == %0.3f) & (chemical == '%s')"%(chem_id,run_number,chemical_name))
                    h5attachattributes(thisnode.default_chemical_id,['TITLE'],['(default concentration, determined by ESR, for a particular chemical name and run number)'])
                    h5fp.close()
        plot(integral,alpha=0.5,linewidth=linewidth,label = legendtext)
        figure_list.next('epr')
        if normalize_peak:
           normalization = abs(data).run_nopop(max,newname)
           data /= normalization
        ax = gca()
        if offset_spectra and index>0:
            myoffset = array(ax.get_ylim()).min()-array(ax.get_ylim()).max() # so I set the zero at the bottom of the plot, which is negative, and also allow space for the top of the spectrum (which is positive
        else:
            myoffset = 0.
        plot_color_counter(pc)
        plot(data+myoffset,alpha=0.5,linewidth=linewidth,label=chemical_name)
        if not normalize_field:
            minval = abs(data.getaxis(newname)-centerfield).argmin()
            centerpoint = data[newname,minval]
            plot_color_counter(pc)
            if not find_maxslope:
                plot(centerfield,centerpoint.data+myoffset,'o',markersize = 5,alpha=0.3)
        axis('tight')
        if index == 0:
            fieldbar *= array(ax.get_ylim()).max()
    plot_color_counter(pc)
    show_errorzone = False
    if show_errorzone:
        plot(data.getaxis(newname)[mask],zeros(shape(data.getaxis(newname)[mask])),'k',alpha=0.2,linewidth=10)
    figure_list.next('epr')
    plot(fieldbar,'k',linewidth = 2.0)
    fieldbar[newname,4:5].plot_labels(['10 G\n'],va = 'bottom',color = 'k',alpha = 1.0,size = 'x-small')
    if grid:
        gridandtick(gca())
    axis('tight')
    figure_list.next('epr_int')
    axis('tight')
    return figure_list
#}}}
