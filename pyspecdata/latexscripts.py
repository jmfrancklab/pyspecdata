r'''Provides the ``pdflatex_notebook_wrapper`` shell/dos command, which you run
instead of your normal Latex command to build a lab notebook.
The results of python environments are **cached** and **only re-run if the code changes**,
even if the python environments are moved around.
This makes the compilation of a Latex lab notebook extremely efficient.
'''
from .datadir import getDATADIR, pyspec_config
#from .datadir import get_notebook_dir
from distutils.spawn import find_executable
import os.path
import hashlib
import numpy
import sys
import time
import platform
from shutil import copy2 as sh_copy2
from subprocess import Popen,PIPE,check_output
haswatchdog = False
if haswatchdog:
    from watchdog.events import FileSystemEventHandler
    from watchdog.observers import Observer
if haswatchdog:
    class handle_dir_change (FileSystemEventHandler):
        def on_modified(self,event):
            #print "called 'on_modified'"
            if event is None:# if I called it manually
                runornot = True
                print("Scons signaled manually")
            else:
                #print "in event:",dir(event)
                #print event.event_type
                #print event.is_directory
                #print event.key
                #print event.src_path
                if event.src_path in ['.' + os.path.sep + j for j in self.dependencies]:
                    print("Scons signaled on",event.src_path,", which is a dependency")
                    runornot = True
                else:
                    print("Modified",event.src_path,"which is not a dependency")
                    runornot = False
            if runornot:
                outstring = check_output(['scons','--tree=all'] + sys.argv[1:],shell = True)
                #print "outstring:\n",outstring
                self.dependencies = []
                for line in outstring.split(os.linesep):
                    # here, I grab the tree, and then I should be able to use watchdog to implement latexmk-like functionality
                    x = line.find('+-')
                    if x > -1:
                        x += 1
                        while line[x] == '-':
                            x += 1
                        self.dependencies.append(line[x:])
                self.dependencies.pop(0) # this is the file itself
                return self.dependencies
def det_new_pdf_name(thisargv):
    'based on an original tex or pdf name, determine the original basename (i.e., no extension), as well as one with the final word after the underscore removed'
    tex_basename = list(filter(lambda x: x[0] != '-',
            thisargv))[-1]
    tex_basename = os.path.basename(tex_basename)
    if tex_basename[-4:] == '.tex':
        tex_basename = tex_basename[:-4]
    elif tex_basename[-4:] == '.pdf':
        tex_basename = tex_basename[:-4]
    orig_tex_basename = tex_basename
    tex_basename = tex_basename.split('_')
    if (len(tex_basename) > 1) and (tex_basename[-1] in ['basic','fancy','georgia']):
        new_pdf_basename = '_'.join(tex_basename[:-1])
    else:
        new_pdf_basename = '_'.join(tex_basename)
    return orig_tex_basename,new_pdf_basename
def genconfig():
    '''creates a template configuration directory'''
    pyspec_config._config_parser = None # this supresses the output
    if platform.platform().startswith('Windows'):
        hide_start = '_' # the default hidden/config starter for vim, mingw applications, etc
    else:
        hide_start = '.'
    filename = os.path.join(os.path.expanduser('~'),hide_start+'pyspecdata')
    with open(filename,'w',encoding='utf-8') as fp:
        fp.write('[General]\n')
        fp.write('# replace the following with your default data location (this is just a suggestion)\n')
        possible_data = [x for x in next(os.walk(os.path.expanduser('~')))[1] if 'data' in x.lower() and 'app' not in x.lower()]
        fp.write('data_directory = ')
        if len(possible_data) > 0:
            fp.write(os.path.join(os.path.expanduser('~'),possible_data[0])+'\n')
        else:
            fp.write(os.path.join(os.path.expanduser('~'),'???')+'\n')
        fp.write('\n')
        fp.write('[ExpTypes]\n')
        fp.write('# in this section, you can place specific subdirectories (for faster access, or if they live elsewhere\n') 
    print("now edit "+filename)

def wraplatex():
    '''runs the python scripts after running latex also creates a copy of latex without the final portion under the underscore
    This prevents the viewer from hanging while it's waiting for a refresh.
    This can be used in combination with wrapviewer() and latexmk by using a ``~/.latexmkrc`` file that looks like this:

    If you pass the ``--xelatex`` argument, xelatex is used instead of pdflatex
    (note that if you're using latexmk, you need to add this in the latexmkrc file).

    .. code-block:: perl

        $pdflatex=q/pdflatex_notebook_wrapper %O -synctex=1 %S/;# calls this function
        $pdf_previewer=q/pdflatex_notebook_view_wrapper/;# calls the wrapviewer function
    '''
    proc_args = list(sys.argv)
    if '--xelatex' in proc_args:
        proc_args.pop(proc_args.index('--xelatex'))
        use_xelatex = True
    else:
        use_xelatex = False
    print("about to update the python script outputs....")
    orig_tex_basename,new_pdf_basename = det_new_pdf_name(proc_args)
    with open(orig_tex_basename+'.tex','r',encoding='utf-8') as fp:
        thisline = fp.readline()
        while thisline.startswith('%!'):# in case we want to allow multiple directives
            if 'xelatex' in thisline:
                use_xelatex = True
            elif 'pdflatex' in thisline:
                use_xelatex = False
            thisline = fp.readline()
    if use_xelatex:
        shellcmd = ' '.join(['xelatex']+proc_args[1:])
    else:
        shellcmd = ' '.join(['pdflatex']+proc_args[1:])
    print("executing:",shellcmd)
    os.system(shellcmd)
    print("executing:",'update_notebook_pythonscripts')
    os.system('update_notebook_pythonscripts')
    if orig_tex_basename != new_pdf_basename:
        print("preparing to:",'cp '+orig_tex_basename+'.pdf '+new_pdf_basename+'.pdf')
        os.system('cp '+orig_tex_basename+'.pdf '+new_pdf_basename+'.pdf')
        os.system('cp '+orig_tex_basename+'.synctex.gz '
                +new_pdf_basename+'.synctex.gz')
    return
def wrapviewer():
    'see :func:`wraplatex <pyspecdata.latexscripts.wraplatex>`'
    print("started wrapviewer!")
    pdf_basename = list(filter(lambda x: x[0] != '-',
            sys.argv))[-1]
    orig_tex_basename,new_pdf_basename = det_new_pdf_name(sys.argv)
    if os.name == 'posix':
        # {{{ this plays the role of the function that I used to call "new_evince" with argument "b"
        full_pdf_name = new_pdf_basename+'.pdf'
        full_tex_name = orig_tex_basename+'.tex'
        cmd = ['zathura --synctex-forward']
        cmd.append(f"1:0:{full_tex_name} {full_pdf_name}")
        print(' '.join(cmd))
        os.system(' '.join(cmd))
        # }}}
    else:
        os.system('start sumatrapdf -reuse-instance '+new_pdf_basename+'.pdf')
    if new_pdf_basename == 'lists':
        os.system('cp lists.pdf "'+os.path.join(os.path.expanduser('~'),'Dropbox','lists.pdf'))
    return
def script_filename(scriptnum_as_str):
    return get_scripts_dir()+scriptnum_as_str+'.py'
def cached_filename(hashstring,returndir = False):
    r'''this sets the format for where the cached file is stored
    we use the first two characters as directory names (so there are just 16 of them'''
    return get_scripts_dir() + 'cache' + os.path.sep + hashstring[0] + os.path.sep + hashstring[1] + os.path.sep + hashstring[2:] + '.tex'
def grab_script_string(scriptnum_as_str):
    fp_script = open(script_filename(scriptnum_as_str),encoding='utf-8')
    #print "opening",script_filename(scriptnum_as_str)
    script_string = ''
    reading = False
    for scriptline in fp_script.readlines():
        if reading:
            script_string += scriptline
        if scriptline.strip() == '### start checksum ###':
            reading = True
    fp_script.close()
    return script_string
def check_image_path():
    image_path = os.path.sep.join([os.getcwd(),'auto_figures'])
    if not os.path.exists(image_path):
        os.mkdir(image_path)
    return
def get_scripts_dir():
    script_path = os.path.sep.join([os.getcwd(),'scripts',''])
    if not os.path.exists(script_path):
        os.mkdir(os.path.sep.join([os.getcwd(),'scripts']))
    return script_path
def sha_string(script):
    'convert the sha hash to a string'
    s = hashlib.sha256()
    s.update(script.encode('utf-8'))
    hasharray = numpy.fromstring(s.digest(),'>u8')
    del s
    return ''.join(['%016x'%x for x in list(hasharray)])
def cache_output_if_needed(scriptnum_as_str,hashstring,showcode = False,show_error = True):
    r'if needed, run the python script numbered by scriptnum_as_str that hashes to hashstring, and output the result to the cache ONLY'
    output_fname = cached_filename(hashstring)
    script_fname = script_filename(scriptnum_as_str)
    # {{{ interpret the "NOerr" directive
    with open(script_fname,'r',encoding='utf-8') as fp:
        firstline = fp.readline()
        if firstline.startswith('### NOerr'): show_error = False
    # }}}
    if os.path.exists(output_fname):
        print(output_fname,"already exists, so I will just use it")
    else:
        print("no cached file")
        if not os.path.exists(os.path.dirname(output_fname)):
            os.makedirs(os.path.dirname(output_fname))
        fp_out = open(output_fname,'w',encoding='utf-8')
        #fp_out.write(r'{\color{red}script: %d}'%script_number+'\n')
        if showcode:
            fp_out.write(r'\begin{lstlisting}'+'\n')
            fp_out.write(''.join(script))
            fp_out.write(r'\end{lstlisting}'+'\n')
        temp = os.environ
        print("about to run python")
        python_name = 'python'
        if os.name == 'posix':
            temp.update({'PYTHON_DATA_DIR':getDATADIR()})
            # on mac, we frequently want to use python2
            if find_executable('python2'):
                python_name = 'python2'
            proc = Popen([python_name,'-W','ignore',script_fname],
                    stdout = PIPE,
                    stdin = PIPE,
                    stderr = PIPE,
                    env = temp)
        else: #windows should give os.name == 'nt'
            temp.update({'MPLCONFIGDIR':os.getcwd()+'/.matplotlib',
                        'PYTHON_DATA_DIR':getDATADIR()})
            proc = Popen([python_name,'-W','ignore',script_fname],
                    stdout = PIPE,
                    stdin = PIPE,
                    stderr = PIPE,
                    env = temp)
        stdoutdata,stderrdata = [j.decode('utf-8') for j in proc.communicate()]
        print("ran python")
        if os.name != 'posix':
            stdoutdata = stdoutdata.replace('\r','')
            if stderrdata is not None:
                stderrdata = stderrdata.replace('\r','')
        #fp_out.write('File has not been run, output:\n\n')
        fp_out.write("%%% Generated automatically by python for script {:s} hashstring {:s}".format(scriptnum_as_str,hashstring))
        fp_out.write("\n")
        fp_out.write(stdoutdata)
        #fp_out.write('\n\n$\\Rightarrow$ end of output\n\n')
        if show_error:
            if stderrdata is not None and len(stderrdata) > 0:
                fp_out.write("\\quad\\\\ {\\small {\\color{red} {\\tt ERRORS---------------------} "+r'\makeatletter\fn{scripts/\thepy@codenum.py}\makeatother'+"}}\\\\\n")
                fp_out.write("\n\nThe current directory is \\verb|%s|\n\n"%os.getcwd())
                fp_out.write("\n\nThe data directory was set to: \\verb|"+getDATADIR()+"|\n\n")
                #fp_out.write("\n\nThe notebook directory was set to: \\verb|"+get_notebook_dir()+"|\n\n")
                fp_out.write("\\begin{tiny}\n")
                fp_out.write("\\begin{verbatim}\n")
                #fp_out.write('...\n\t'.join(textwrap.wrap(stderrdata,80)))
                fp_out.write(stderrdata)
                fp_out.write("\\end{verbatim}\n")
                fp_out.write("\\end{tiny}\n")
                fp_out.write("{\\small {\\color{red} {\\tt ---------------------------}}}\\\\\n")
        fp_out.close()
        return
def flush_script(number):
    tex_name = os.path.normpath(os.path.join(os.getcwd(),'scripts',number+'.tex'))
    print("removing:",tex_name)
    if os.path.exists(tex_name):
        os.remove(tex_name)
    file_name = cached_filename(sha_string(grab_script_string(number)))
    print("removing:",file_name)
    if os.path.exists(file_name):
        os.remove(file_name)
        print("does it still exist?",os.path.exists(file_name))
    else:
        print("there is no cache for script",number)
    return
if haswatchdog:
    def repeat_scons():
        r'This is just a function to run scons over and over again (e.g. so it can replace latexmk)'
        for j in os.listdir('.'):
            if j.endswith('.sconslog'):
                os.remove(j)# by removing these, we force the PDF reader to open
        observer = Observer()
        myhandler = handle_dir_change()
        dependencies = myhandler.on_modified(None) # run once manually
        for j in dependencies:
            fname = '.'+ os.path.sep + j
            #print "adding",os.path.dirname(fname),"for",fname
            observer.schedule(myhandler,path = os.path.dirname(fname),recursive = False)
        observer.start()
        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            observer.stop()
        observer.join()
        return
def main():
    r'''This looks for `scripts/scriptsUsed.csv` inside the notebook directory, and checks whether or not it should be run
    if a command line argument of "flush" is passed, it flushes that script number from the cache'''
    check_image_path()
    if len(sys.argv) > 2:
        if sys.argv[1] == 'flush':
            if len(sys.argv) == 3:
                flush_script(sys.argv[2])
                exit()
            elif len(sys.argv) == 4:
                for j in range(int(sys.argv[2]),int(sys.argv[3])+1):
                    flush_script(str(j))
                exit()
            else:
                raise RuntimeError("What did you pass me???")
    if not os.path.exists(get_scripts_dir()): return
    fp = open(get_scripts_dir() + 'scriptsUsed.csv')
    for line in fp.readlines():
        scriptnum_as_str = line.strip()
        outname    = scriptnum_as_str + '.tex'
        hashstring = sha_string(grab_script_string(scriptnum_as_str))
        print("reading script",script_filename(scriptnum_as_str),"which has hash",hashstring)
        cache_output_if_needed(scriptnum_as_str,hashstring)
        # always pull the  output from the cache
        sh_copy2(cached_filename(hashstring),get_scripts_dir()+outname)
    fp.close()
