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
import configparser
import importlib
import subprocess
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
    '''creates or edits the configuration directory'''
    pyspec_config._config_parser = None # this supresses the output
    if platform.platform().startswith('Windows'):
        hide_start = '_' # the default hidden/config starter for vim, mingw applications, etc
    else:
        hide_start = '.'
    filename = os.path.join(os.path.expanduser('~'),hide_start+'pyspecdata')
    qt_widgets = None
    # The GUI is optional, so the code attempts to import any available Qt binding.
    for module_name in ['PyQt5','PySide2','PyQt6','PySide6']:
        try:
            qt_widgets = importlib.import_module(module_name+'.QtWidgets')
            break
        except ImportError:
            qt_widgets = None
            continue
    if qt_widgets is None:
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
        return

    # At this point Qt is available, so load the existing configuration to build the editor.
    config_parser = configparser.ConfigParser()
    if os.path.exists(filename):
        config_parser.read(filename,encoding='utf-8')
    if not config_parser.has_section('General'):
        config_parser.add_section('General')
    if not config_parser.has_section('ExpTypes'):
        config_parser.add_section('ExpTypes')
    if not config_parser.has_section('RcloneRemotes'):
        config_parser.add_section('RcloneRemotes')

    try:
        remote_output = subprocess.check_output(['rclone','listremotes'],stderr=subprocess.STDOUT)
        remote_output = remote_output.decode('utf-8')
        remote_candidates = [j.strip() for j in remote_output.splitlines() if j.strip() != '']
        rclone_remotes = [j[:-1] if j.endswith(':') else j for j in remote_candidates]
    except Exception:
        rclone_remotes = []

    # The dialog is defined inside the function so that the module can be imported without Qt.
    class ConfigDialog(qt_widgets.QDialog):
        def __init__(self, this_config, config_path, remote_names):
            qt_widgets.QDialog.__init__(self)
            self.setWindowTitle('pyspecdata data configuration')
            self.resize(900,600)
            self.config_parser = this_config
            self.config_path = config_path
            self.remote_names = remote_names
            self.file_rows = []
            self.general_entries = []
            self.build_interface()

        def build_interface(self):
            'construct the entire dialog'
            main_layout = qt_widgets.QVBoxLayout(self)
            self.tab_widget = qt_widgets.QTabWidget()
            main_layout.addWidget(self.tab_widget)
            self.create_files_tab()
            self.create_variables_tab()
            button_layout = qt_widgets.QHBoxLayout()
            button_layout.addStretch(1)
            self.save_button = qt_widgets.QPushButton('Save')
            self.cancel_button = qt_widgets.QPushButton('Cancel')
            button_layout.addWidget(self.save_button)
            button_layout.addWidget(self.cancel_button)
            main_layout.addLayout(button_layout)
            self.save_button.clicked.connect(self.handle_save)
            self.cancel_button.clicked.connect(self.reject)

        def create_files_tab(self):
            'create the tab that holds data directories and experiment types'
            files_tab = qt_widgets.QWidget()
            files_layout = qt_widgets.QVBoxLayout(files_tab)
            path_layout = qt_widgets.QHBoxLayout()
            path_label = qt_widgets.QLabel('data_directory')
            path_layout.addWidget(path_label)
            self.data_directory_edit = qt_widgets.QLineEdit()
            if self.config_parser.has_option('General','data_directory'):
                self.data_directory_edit.setText(self.config_parser['General']['data_directory'])
            else:
                possible_data = [x for x in next(os.walk(os.path.expanduser('~')))[1] if 'data' in x.lower() and 'app' not in x.lower()]
                if len(possible_data) > 0:
                    self.data_directory_edit.setText(os.path.join(os.path.expanduser('~'),possible_data[0]))
                else:
                    self.data_directory_edit.setText(os.path.join(os.path.expanduser('~'),'???'))
            path_layout.addWidget(self.data_directory_edit)
            browse_button = qt_widgets.QToolButton()
            browse_button.setIcon(self.style().standardIcon(qt_widgets.QStyle.SP_DirOpenIcon))
            path_layout.addWidget(browse_button)
            browse_button.clicked.connect(self.select_data_directory)
            files_layout.addLayout(path_layout)

            self.files_scroll = qt_widgets.QScrollArea()
            self.files_scroll.setWidgetResizable(True)
            self.files_container = qt_widgets.QWidget()
            self.files_layout = qt_widgets.QVBoxLayout(self.files_container)
            self.files_layout.setContentsMargins(0,0,0,0)
            self.files_layout.setSpacing(6)

            existing_types = []
            for key,value in self.config_parser.items('ExpTypes'):
                existing_types.append((key,value))
            remote_info = {}
            for key,value in self.config_parser.items('RcloneRemotes'):
                remote_info[key] = value
            exp_keys = [j for j,_ in existing_types]
            for key in remote_info.keys():
                if key not in exp_keys:
                    existing_types.append((key,''))
            if len(existing_types) == 0:
                self.add_file_row('','','None','')
            else:
                for key,value in existing_types:
                    if key in remote_info:
                        if ':' in remote_info[key]:
                            remote_name,remote_path = remote_info[key].split(':',1)
                        else:
                            remote_name = remote_info[key]
                            remote_path = ''
                    else:
                        remote_name = 'None'
                        remote_path = ''
                    self.add_file_row(key,value,remote_name,remote_path)
            spacer = qt_widgets.QWidget()
            spacer.setSizePolicy(qt_widgets.QSizePolicy.Minimum,qt_widgets.QSizePolicy.Expanding)
            self.files_layout.addWidget(spacer)
            self.files_scroll.setWidget(self.files_container)
            files_layout.addWidget(self.files_scroll)
            self.add_button = qt_widgets.QPushButton('Add Entry')
            files_layout.addWidget(self.add_button)
            self.add_button.clicked.connect(lambda: self.add_file_row('','','None',''))
            self.tab_widget.addTab(files_tab,'Files')

        def create_variables_tab(self):
            'create the tab that holds the general configuration variables'
            variables_tab = qt_widgets.QWidget()
            variables_layout = qt_widgets.QVBoxLayout(variables_tab)
            self.variables_scroll = qt_widgets.QScrollArea()
            self.variables_scroll.setWidgetResizable(True)
            self.variables_container = qt_widgets.QWidget()
            self.variables_layout = qt_widgets.QVBoxLayout(self.variables_container)
            self.variables_layout.setContentsMargins(0,0,0,0)
            self.variables_layout.setSpacing(6)
            for key,value in self.config_parser.items('General'):
                if key == 'data_directory':
                    continue
                self.add_general_row(key,value)
            self.variables_scroll.setWidget(self.variables_container)
            variables_layout.addWidget(self.variables_scroll)
            self.tab_widget.addTab(variables_tab,'Variables')

        def add_general_row(self,key,value):
            'add a single general variable row'
            row_widget = qt_widgets.QWidget()
            row_layout = qt_widgets.QHBoxLayout(row_widget)
            row_layout.setContentsMargins(0,0,0,0)
            key_edit = qt_widgets.QLineEdit()
            key_edit.setText(key)
            value_edit = qt_widgets.QLineEdit()
            value_edit.setText(value)
            row_layout.addWidget(key_edit)
            row_layout.addWidget(value_edit)
            self.variables_layout.addWidget(row_widget)
            self.general_entries.append({'key_edit':key_edit,'value_edit':value_edit})

        def add_file_row(self,name,path,remote_name,remote_path):
            'add a single experiment type row to the scrolling list'
            row_widget = qt_widgets.QWidget()
            row_layout = qt_widgets.QHBoxLayout(row_widget)
            row_layout.setContentsMargins(0,0,0,0)
            name_edit = qt_widgets.QLineEdit()
            name_edit.setText(name)
            path_edit = qt_widgets.QLineEdit()
            path_edit.setText(path)
            browse_button = qt_widgets.QToolButton()
            browse_button.setIcon(self.style().standardIcon(qt_widgets.QStyle.SP_DirOpenIcon))
            remote_combo = qt_widgets.QComboBox()
            combo_entries = ['None']
            for remote_entry in self.remote_names:
                if remote_entry not in combo_entries:
                    combo_entries.append(remote_entry)
            if remote_name not in combo_entries and remote_name not in [None,'']:
                combo_entries.append(remote_name)
            for entry in combo_entries:
                remote_combo.addItem(entry)
            if remote_name in combo_entries:
                remote_combo.setCurrentIndex(combo_entries.index(remote_name))
            remote_path_edit = qt_widgets.QLineEdit()
            remote_path_edit.setText(remote_path)
            row_layout.addWidget(name_edit)
            row_layout.addWidget(path_edit)
            row_layout.addWidget(browse_button)
            row_layout.addWidget(remote_combo)
            row_layout.addWidget(remote_path_edit)
            self.files_layout.insertWidget(self.files_layout.count()-1,row_widget)
            self.file_rows.append({'name_edit':name_edit,'path_edit':path_edit,'browse_button':browse_button,'remote_combo':remote_combo,'remote_path_edit':remote_path_edit})
            browse_button.clicked.connect(lambda: self.select_path(path_edit))
            remote_combo.currentIndexChanged.connect(lambda _: self.update_remote_state(remote_combo,remote_path_edit))
            self.update_remote_state(remote_combo,remote_path_edit)

        def update_remote_state(self,combo,remote_edit):
            'enable or disable the remote path field based on the combo box'
            if combo.currentText() == 'None':
                remote_edit.setEnabled(False)
                remote_edit.setText('')
            else:
                remote_edit.setEnabled(True)

        def select_data_directory(self):
            'browse for the base data directory'
            chosen = qt_widgets.QFileDialog.getExistingDirectory(self,'Select data directory',self.data_directory_edit.text())
            if chosen:
                self.data_directory_edit.setText(chosen)

        def select_path(self,target_edit):
            'browse for an experiment directory'
            chosen = qt_widgets.QFileDialog.getExistingDirectory(self,'Select directory',target_edit.text())
            if chosen:
                target_edit.setText(chosen)

        def handle_save(self):
            'collect the contents of the dialog and write the configuration file'
            new_config = configparser.ConfigParser()
            new_config.add_section('General')
            new_config.set('General','data_directory',self.data_directory_edit.text().strip())
            for entry in self.general_entries:
                key_text = entry['key_edit'].text().strip()
                if key_text == '' or key_text == 'data_directory':
                    continue
                new_config.set('General',key_text,entry['value_edit'].text().strip())
            new_config.add_section('ExpTypes')
            remote_values = {}
            for row in self.file_rows:
                name_text = row['name_edit'].text().strip()
                path_text = row['path_edit'].text().strip()
                if name_text == '' or path_text == '':
                    continue
                new_config.set('ExpTypes',name_text,path_text)
                combo_value = row['remote_combo'].currentText()
                if combo_value != 'None':
                    remote_path_text = row['remote_path_edit'].text().strip()
                    if remote_path_text != '':
                        remote_values[name_text] = combo_value+':'+remote_path_text
            new_config.add_section('RcloneRemotes')
            for key in sorted(remote_values.keys()):
                new_config.set('RcloneRemotes',key,remote_values[key])
            for section in self.config_parser.sections():
                if section in ['General','ExpTypes','RcloneRemotes']:
                    continue
                if not new_config.has_section(section):
                    new_config.add_section(section)
                for key,value in self.config_parser.items(section):
                    new_config.set(section,key,value)
            with open(self.config_path,'w',encoding='utf-8') as fp:
                new_config.write(fp)
            self.accept()

    existing_app = qt_widgets.QApplication.instance()
    if existing_app is None:
        existing_app = qt_widgets.QApplication(sys.argv)
    dialog = ConfigDialog(config_parser,filename,rclone_remotes)
    dialog.exec_()

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
    print("executing:",shellcmd)
    os.system(shellcmd)
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
