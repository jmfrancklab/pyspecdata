import sys
import os
import re
from PyQt4 import QtGui, QtCore
from install_gui import Ui_Dialog
from py_compile import compile
from subprocess import Popen,PIPE,call
from python.datadir import getDATADIR
import os
#{{{ stuff for explicit latex generation
def dump_latex_preamble(fp):
    fp.write(r'''\documentclass[10pt]{book}
\usepackage{mynotebook}
\begin{document}
    ''')
def dump_latex_finish(fp):    
    fp.write(r'''
\end{document}
    ''')
#}}}
# found this here http://hohohuhu.blogspot.com/search/label/PyQt
def get_homedir():# taken from latexscript.py
    return os.path.expanduser('~')
app = QtGui.QApplication(sys.argv)
class my_widget_class (QtGui.QDialog):
    # here, I use the QDialog class, which has accept and reject, and I add the following custom routines, which I can call as slots
    def my_complete(self):
        print "this is the my complete function"
        if not os.path.exists('auto_figures'):
            os.mkdir('auto_figures')
        if self.datadir_changed:
            if os.path.exists(os.curdir+'/.datadir'):
                fp = open(os.curdir+'/.datadir','r')
                print '.datadir originally said:',fp.read()
                fp.close()
            fp = open(os.curdir+'/.datadir','w')
            fp.write(self.currently_displayed_datadir)
            fp = open('.datadir','r')
            print 'now, .datadir says:',fp.read()
            fp.close()
        if os.path.exists(self.currently_displayed_homedir+'/.latexmkrc'):
            print 'I see that you already have a',self.currently_displayed_homedir+'/.datadir, so I\'m not going to overwrite it'
        else:
            print "I see that no latexmkrc exists in",self.currently_displayed_homedir,"so I will create one"
            QtGui.QMessageBox.warning(self,"testing installation","Since this seems to be the first time that you are running stuff, I'm going to try to open sumatrapdf --> if it comes up, just close it.")
            if os.name == 'posix':#assume linux
                fp = open(self.currently_displayed_homedir+'/.latexmkrc','w')
                fp.write('''$pdflatex=q/pdflatex %O -synctex=1 -shell-escape %S/;
$pdf_previewer=q/okular --unique/;''')
                fp.close()
            else:#assume it's windows
                exit_fail = False
                #{{{ test for sumatrapdf
                pobj = Popen("start sumatrapdf",shell = True)
                pobj.wait()
                print "sumatrapdf return code is ",pobj.returncode
                if pobj.returncode > 0:
                    QtGui.QMessageBox.warning(self,"Can't find","Actually, I can't find sumatrapdf, so you should probably install it from <a href='http://blog.kowalczyk.info/software/sumatrapdf/download-free-pdf-viewer.html'>here</a>")
                    exit_fail = True
                #}}}
                if not exit_fail:
                    fp = open(self.currently_displayed_homedir+'/.latexmkrc','w')
                    fp.write(r'$pdflatex=q/pdflatex %O -synctex=1 -shell-escape %S/;'\
                            +'\n'\
                            +r'$pdf_previewer=q/start sumatrapdf -reuse-instance/;''')
                    fp.close()
                    print 'now, I will also compile '+os.getcwd()+'/latexscript.py...'
                    compile(os.getcwd()+'/latexscript.py')
    def my_initialize_directories(self):
        self.emit(QtCore.SIGNAL("my_change_datadir(QString)"),getDATADIR())
        self.emit(QtCore.SIGNAL("my_change_homedir(QString)"),get_homedir())
        self.currently_displayed_datadir = getDATADIR()
        self.currently_displayed_homedir = get_homedir()
        self.datadir_changed = False
        self.homedir_changed = False
        self.script_to_run = None
    def my_change_selection(self,arg):
        print "this is the my change selection function, with argument:",arg
        if arg == "Set Data Directory":
            print "you want to change the directory"
            temp = str(QtGui.QFileDialog.getExistingDirectory(self, "Select Data Directory",self.currently_displayed_datadir))
            if len(temp)>0:
                self.currently_displayed_datadir = temp
                print "and you selected",self.currently_displayed_datadir
            else:
                print "you didn't select anything"
            os.environ['PYTHONDATADIR'] = self.currently_displayed_datadir
            self.datadir_changed = True
            self.emit(QtCore.SIGNAL("my_change_datadir(QString)"),self.currently_displayed_datadir)
        elif arg == "Set Home Directory":
            print "you want to set the home directory"
            temp = str(QtGui.QFileDialog.getExistingDirectory(self, "Select Home Directory",self.currently_displayed_homedir))
            if len(temp)>0:
                self.currently_displayed_homedir = temp
                print "and you selected",self.currently_displayed_homedir
            else:
                print "you didn't select anything"
            self.homedir_changed = True
            self.emit(QtCore.SIGNAL("my_change_homedir(QString)"),self.currently_displayed_homedir)
        elif arg == "Run Script":
            print "you want to run a script"
            curdir = os.getcwd()
            print "the current directory is",curdir
            temp = str(QtGui.QFileDialog.getOpenFileName(self,
                "Select a script to run",
                directory = curdir,
                filter = "Python script (*.py);;Tex file (*.tex)"
                ))
            if len(temp)>0:
                self.script_to_run = temp
            else:
                print "you didn't select anything"
            if self.script_to_run is not None:
                filename_re = re.compile(r'(.*).py$')
                m = filename_re.match(self.script_to_run)
                if m:
                    file_to_gen_basename = m.groups()[0]
                    file_to_gen = file_to_gen_basename+'.tex'
                    print "now, I will run the script",self.script_to_run,"to generate the file",file_to_gen
                    if os.path.exists(file_to_gen):
                        QtGui.QMessageBox.warning(self,"overwrite","The file %s exists, so I'm going to overwrite it"%file_to_gen)
                    if os.name == 'posix':
                        python_arg_dict = {'PYTHONPATH':os.getcwd(),
                                    'PYTHONDATADIR':os.environ['PYTHONDATADIR']}
                    else: #windows should give os.name == 'nt'
                        python_arg_dict = {'PYTHONPATH':os.getcwd(),
                                'SystemRoot':'C:\\Windows',
                                'MPLCONFIGDIR':os.getcwd()+'/.matplotlib',
                                'PATH':'',
                                'PYTHONDATADIR':os.environ['PYTHONDATADIR']}
                    print "about to run"
                    proc = Popen(['python','-W','ignore',self.script_to_run],
                            stdout = PIPE,
                            stdin = PIPE,
                            stderr = PIPE,
                            env = python_arg_dict)
                    stdoutdata,stderrdata = proc.communicate()
                    print "done running"
                    fp_out = open(file_to_gen,'w')
                    dump_latex_preamble(fp_out)
                    fp_out.write(stdoutdata)
                    if stderrdata is not None and len(stderrdata) > 0:
                        #QtGui.QMessageBox.warning(self,"Errors!","generated the following errors:\n%s"%stderrdata)
                        mb = QtGui.QMessageBox(self)
                        mb.setText("There were errors while running your script!\nPlease check the details, which should tell you exactly what's going on, so that you know how to fix it!!!")
                        mb.setDetailedText(str(stderrdata))
                        mb.exec_()
                        fp_out.write("\\quad\\\\ {\\small {\\color{red} {\\tt ERRORS---------------------}}}\\\\\n")
                        fp_out.write("\\begin{verbatim}\n")
                        fp_out.write(stderrdata)
                        fp_out.write("\\end{verbatim}\n")
                        fp_out.write("{\\small {\\color{red} {\\tt ---------------------------}}}\\\\\n")
                    else:
                        QtGui.QMessageBox.warning(self,"Done!","Script ran without any errors!")
                    dump_latex_finish(fp_out)
                    fp_out.close()
                    proc = Popen(['pdflatex','-synctex=1','-shell-escape',file_to_gen],
                            shell = True)
                    proc.wait()
                    proc = Popen(['bibtex',file_to_gen_basename],
                            shell = True)
                    proc.wait()
                    proc = Popen(['pdflatex',
                        '-synctex=1',
                        '-shell-escape',
                        file_to_gen],shell = True)
                    proc.wait()
                    proc = Popen(['start',
                        'sumatrapdf',
                        '-reuse-instance',
                        file_to_gen_basename+'.pdf'],shell = True)
                    proc.wait()
                else:
                    raise ValueError("You selected a file that's not a python script!!")
        elif arg == "---Select Action to Perform---":
            print "detected a reset to the original combo box value"
        self.text_to_change_to = "---Select Action to Perform---"
        self.emit(QtCore.SIGNAL("my_change_list_selection(QString)"),self.text_to_change_to) # this actually doesn't work, and I have no idea why
        print "this should have changed the selection to ---Select Action to Perform---"
        return
widget = my_widget_class()
dialog = Ui_Dialog()
dialog.setupUi(widget)
def my_accept_callback():
    print 'this is my accept callback!!'
#widget.connect(my_widget_class,QtCore.SIGNAL('accept()'),my_accept_callback)
widget.my_initialize_directories()
widget.show()
sys.exit(app.exec_())
