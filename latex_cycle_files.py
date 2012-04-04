import os
import sys
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
