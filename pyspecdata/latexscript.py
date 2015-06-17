import sys
import hashlib
import numpy
import textwrap
import os
from matplotlib.mlab import rec2csv, csv2rec
from difflib import ndiff
from datadir import grab_data_directory
from subprocess import Popen,PIPE
def main(show_error = True):
    script_catalog_type = numpy.dtype([('number','>u8'),('hash','|S72')])
    grabbed_datadir_from_file,datadir_error = grab_data_directory()
    def determine_file_catalog_number(filename):
        r'a function that returns the index number associated with a given file; if necessary, it will add the file to the underlying csv file'
        return
    def read_dependency_dict():    
        r'''return a dict of form {key1,list1,...} based on dep_catalog.csv, which is stored in the format:
        key1,list1[0],...,list1[N]
        key2,list2[0],...,list2[N]
        ...
        keyM,listM[0],...,listM[N]'''
        return
    def add_dependency(scriptnumber,latexfile):
        r'mark down that script "filename" uses the script number "scriptnumber"'
        latexfile_index = determine_file_catalog_number(latexcatalog)
        dependencydict = read_dependency_dict()
    def read_whole_file(filename):
        r'returns the whole file, as a list of lines'
        if type(filename) is str:
            fp = open(filename,'r')
        else:
            fp = filename
        retval = fp.readlines()
        fp.close()
        return retval
    def make_shell_escape_warning():
        fp = open('scripts/warn_shell_escape.tex','w')
        fp.write('{\\color{red}{\\bf For the code to work, you need to run pdflatex with shell escape!}}')
        fp.close()
    def sha_string(script):
        'convert the sha hash to a string'
        s = hashlib.sha256()
        s.update(script)
        hasharray = numpy.fromstring(s.digest(),'>u8')
        del s
        return ''.join(map(lambda x: '%016x'%x,list(hasharray)))
    def get_catalog(filename = 'script_catalog.csv',this_dtype = script_catalog_type):
        'get the catalog file'
        #{{{ just generate auto figures here
        if not os.path.exists('auto_figures/'):
            os.mkdir('scripts')
        #}}}
        if os.path.exists('scripts/'):
            if not os.path.exists('scripts/warn_shell_escape.tex'):
                make_shell_escape_warning()
            if os.path.exists('scripts/'+filename):
                return csv2rec('scripts/'+filename)
            else:
                return numpy.array([],this_dtype)
        else:
            os.mkdir('scripts')
            make_shell_escape_warning()
            return numpy.array([],this_dtype)
    def determine_last_index(catalog):
        'return the last index used to this point'
        try:
            indeces = catalog['number']
        except:
            raise RuntimeError('Problem finding indeces for'+repr(catalog))
        if len(indeces) > 0:
            max_index = indeces.max()
            if max_index > len(indeces): #then I have skipped one
                for j in range(1,max_index+1):
                    if j not in indeces:
                        return j-1
            else:
                return max_index
        else:
            return 0
    def see_if_run(script):
        'see if the file is cataloged, and add to the catalog, then return a flag that tells if it has been run already, the script number, and the catalog'
        catalog = get_catalog()
        hashstring = sha_string(script)
        if hashstring in catalog['hash']:
            retval = [catalog['number'][x] for x in range(0,len(catalog)) if catalog['hash'][x] == hashstring]
            print "I found",hashstring,"in catalog as number",retval
            if len(retval) > 1:
                raise RuntimeError('For some reason, there are two values with this hash number!!')
            return True,retval[0],catalog,hashstring
        else: # this file has not been run before, and needs to be added
            print "I didn't find ",hashstring,"in catalog, running it"
            thisindex = determine_last_index(catalog) + 1
            catalog = catalog.tolist()
            catalog.append(numpy.array((thisindex,hashstring),script_catalog_type))
            #{{{ write the script
            fp_script = open('scripts/%d.py'%thisindex,'w')
            fp_script.write('import pyspecdata.paramset;pyspecdata.paramset.myparams[\'figlist_type\'] = \'figlistl\';from pyspecdata.h5nmr import *;from pyspecdata.acert_hdf5 import *\n')
            fp_script.write(script)
            fp_script.close()
            #}}}
            return False,thisindex,numpy.array(catalog,script_catalog_type),hashstring
    arglist = sys.argv[1:]
    if arglist[0] == 'NOERR':
        show_error = False
        arglist = arglist[1:]
    if len(arglist) == 4:
        texfile,jobnumber,scriptfile,file_for_script_name = arglist
        showcode = False
    elif len(arglist) == 5:
        texfile,jobnumber,scriptfile,file_for_script_name,showcode = arglist
    else:
        raise ValueError('''you need 4 [5] command line arguments to run this:
            texfile,jobnumber,scriptfile,file_for_script_name[,showcode]''')
    script = read_whole_file(scriptfile)
    has_been_run,script_number,catalog,hashstring = see_if_run(''.join(script))
    script_fname = 'scripts/%d.py'%script_number
    output_fname = script_fname + '.out'
    intel_envvariables = os.environ # copy all environment variables, because otherwise some versions of windows cant find the DLL
    if has_been_run:
        #{{{ check that the script matches what's stored here
        stored_script = read_whole_file(script_fname)[1:]
        if not stored_script == script: # because the inclusions are in the first line
            fp_out = open(output_fname,'w')
            fp_out.write(r"{\color{red} Warning!} the previously saved script file for \#%d doesn't match what's in the latex file; possibly, you edited the script file (via synctex), but this is not yet supported"%script_number)# this goes to normal standard output, not into the file!
            fp_out.write('\nThe diff is:\\begin{verbatim}\n%s\\end{verbatim}\n'%(''.join(ndiff(script,stored_script))).replace('\\','\\\\'))
            fp_out.close()
        #}}}
        print "read the output"
    else:
        fp_out = open(output_fname,'w')
        #fp_out.write(r'{\color{red}script: %d}'%script_number+'\n')
        if showcode:
            fp_out.write(r'\begin{lstlisting}'+'\n')
            fp_out.write(''.join(script))
            fp_out.write(r'\end{lstlisting}'+'\n')

        if os.name == 'posix':
            temp = intel_envvariables
            temp.update({'PYTHONPATH':os.getcwd(),
                        'PYTHONDATADIR':os.environ['PYTHONDATADIR']})
            proc = Popen(['python','-W','ignore',script_fname],
                    stdout = PIPE,
                    stdin = PIPE,
                    stderr = PIPE,
                    env = temp)
        else: #windows should give os.name == 'nt'
            temp = intel_envvariables
            temp.update({'PYTHONPATH':os.getcwd(),
                        'MPLCONFIGDIR':os.getcwd()+'/.matplotlib',
                        'PYTHONDATADIR':os.environ['PYTHONDATADIR']})
            proc = Popen(['python','-W','ignore',script_fname],
                    stdout = PIPE,
                    stdin = PIPE,
                    stderr = PIPE,
                    env = temp)
        stdoutdata,stderrdata = proc.communicate()
        if os.name != 'posix':
            stdoutdata = stdoutdata.replace('\r','')
        #fp_out.write('File has not been run, output:\n\n')
        fp_out.write(stdoutdata)
        #fp_out.write('\n\n$\\Rightarrow$ end of output\n\n')
        if show_error:
            if datadir_error or (stderrdata is not None and len(stderrdata) > 0):
                fp_out.write("\\quad\\\\ {\\small {\\color{red} {\\tt ERRORS---------------------} \\fn{%s}}}\\\\\n"%script_fname)
                fp_out.write("\n\nThe current directory is \\verb|%s|\n\n"%os.getcwd())
                if datadir_error:
                    fp_out.write(datadir_error)
                else:
                    fp_out.write("\n\nThe data directory was"+"not"*(not grabbed_datadir_from_file)+"read from .datadir\n\n")
                    fp_out.write("\\begin{tiny}\n")
                    fp_out.write("\\begin{verbatim}\n")
                    #fp_out.write('...\n\t'.join(textwrap.wrap(stderrdata,80)))
                    fp_out.write(stderrdata)
                    fp_out.write("\\end{verbatim}\n")
                    fp_out.write("\\end{tiny}\n")
                fp_out.write("{\\small {\\color{red} {\\tt ---------------------------}}}\\\\\n")
        fp_out.close()
        print "ran and read the output"
    rec2csv(catalog,'scripts/script_catalog.csv')
    fp_scriptname = open(file_for_script_name,'w')
    fp_scriptname.write(script_fname)
    fp_scriptname.close()
    fp_pynum = open('pythonjobname.txt','w')
    fp_pynum.write(hashstring)
    fp_pynum.close()
    #print '%s\t%s'%(m.hexdigest(),thisfile)
    #del m

def mv(src,dest):
    if os.path.isfile(dest):
        os.remove(dest)
    if os.path.isfile(src):
        os.rename(src,dest)
    else:
        fp = open(dest,'w')
        fp.write('')
        fp.close()
def cp(src,dest): # it really seems like there should be an easier way to do this, but I can't find it in os
    fp2 = open(dest,'w')
    if os.path.isfile(src):
        fp = open(src,'r')
        fp2.write(fp.read())
        fp.close()
    else:
        fp2.write('')
    fp2.close()
def cycle_files():
    script_name, whichop = sys.argv[1:]
    if whichop == 'presentoutput':
        mv(script_name,script_name+'.old')
        cp(script_name + '.out',script_name)
    elif whichop == 'storeoutput':
        cp(script_name+'.old',script_name)
        #{{{ this way, by default, the script name is set to the shell escape warning
        fp = open('scriptname.txt','w')
        fp.write('scripts/warn_shell_escape.tex')
        fp.close()
        #}}}
    else:
        raise ValueError('This script must be called with a second argument that\'s either presentoutput or storeoutput')
