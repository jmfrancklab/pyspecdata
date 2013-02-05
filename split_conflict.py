#again rerun
from lxml import html,etree
import os
from matlablike import *
from unidecode import unidecode
import re
import sys
fp = open(sys.argv[1],'r')
needsspace_re = re.compile(r'(\w[):;"-\.,!}]*) +([^%}~])')
text_list = []
print 'opened',sys.argv[1]
log_in_head = True
log_in_new = True
head_title = None
new_title = None
#{{{ pull out just the part between the document text
head_textlist = []
new_textlist = []
j = 0
for thisline in fp:
    j+=1
    if thisline[-11:] == '%FIRSTSET%\n': # if the first set, treat like it's not a comment
        if log_in_head:
            head_textlist.append(thisline)
        if log_in_new:
            new_textlist.append(thisline)
    else:
        if (thisline[:7] == '<<<<<<<'):
            log_in_head = True
            log_in_new = False
            if (head_title is None): # for the first marker, store the title
                    head_title = thisline[7:]
            elif thisline[7:] == head_title:
                pass
            else:
                raise ValueError("I don't understand line %d, which seems to give an inconsistent head title.  It gave:\n%s\nvs expected:\n%s"%(j,thisline[7:],head_title))
        elif (thisline[:7] == '>>>>>>>'):
            log_in_head = True
            log_in_new = True
            if (new_title is None): # for the first marker, store the title
                    new_title = thisline[7:]
            elif thisline[7:] == new_title:
                pass
            else:
                raise ValueError("I don't understand line %d, which seems to give an inconsistent new title.  It gave:\n%s\nvs expected:\n%s"%(j,thisline[7:],new_title))
        elif (thisline[:7] == '======='):
            log_in_head = False
            log_in_new = True
        else:
            if log_in_head:
                head_textlist.append(thisline)
            if log_in_new:
                new_textlist.append(thisline)
fp.close()
#{{{ write out the result
newfile = re.sub(r"(.*)",r'\1.merge_new',sys.argv[1]) 
fp = open(newfile,'w')
new_textlist = ['#%%%%%BRANCH TITLE (This side is saved): '+new_title] + new_textlist
fp.write(''.join(new_textlist))
fp.close()
newfile = re.sub(r"(.*)",r'\1.merge_head',sys.argv[1]) 
fp = open(newfile,'w')
head_textlist = ['#%%%%%BRANCH TITLE: '+head_title] + head_textlist
fp.write(''.join(head_textlist))
fp.close()
#}}}
