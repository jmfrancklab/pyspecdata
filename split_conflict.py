#again rerun
from lxml import html,etree
import os
import re
import sys
needsspace_re = re.compile(r'(\w[):;"-\.,!}]*) +([^%}~])')
allmarkers_re = re.compile(r'(.*?)%(\[ORIG|ORIG\]\[NEW|NEW\])%(.*\n*)')# in all these, the newline at the end is just so it doesn't gobble up the newline
def parse_line(orig_log_text,new_log_text,log_in_orig,log_in_new,text_to_parse):
    match = allmarkers_re.match(text_to_parse)
    if match:
        textone,typeofmarker,texttwo = match.groups()
        if typeofmarker == 'NEW]':
            print "found a matching line, current status is (",log_in_orig,',',log_in_new,")"
            if log_in_new and not log_in_orig:
                switchto = (True,True)# log in orig, log in new
                print "in text:\n",text_to_parse,"\n--> encountered an end marker, switching to",switchto
            else:
                raise ValueError("I encountered an %NEW]% marker, but I was not leaving orig along and logging only in new (False,True), but rather "+repr(log_in_orig)+','+repr(log_in_new)+":\n"+text_to_parse)
        elif typeofmarker == 'ORIG][NEW':
            print "found a matching line, current status is (",log_in_orig,',',log_in_new,")"
            if log_in_orig and not log_in_new:
                switchto = (False,True)# log in orig, log in new
                print "in text:\n",text_to_parse,"\n--> encountered a middle marker, switching to",switchto
            else:
                raise ValueError("I encountered an %ORIG][NEW% marker, but I was not logging in orig but not in new, but rather "+repr(log_in_orig)+','+repr(log_in_new),":\n",text_to_parse)
        elif typeofmarker == '[ORIG':
            print "found a matching line, current status is (",log_in_orig,',',log_in_new,")"
            if log_in_new and log_in_orig:
                switchto = (True,False)# log in orig, log in new
                print "in text:\n",text_to_parse,"\n--> encountered an %[ORIG% marker, switching to",switchto
            else:
                raise ValueError("I encountered an %[ORIG% marker, but I was not logging in both orig and new, but rather"+repr(log_in_orig)+','+repr(log_in_new)+":\n"+text_to_parse)
    else:
        textone = text_to_parse
        texttwo = None
    #}}} check to see if I have a separator
    # regardless, dump the first group into the current bin 
    if log_in_orig:
        orig_log_text += textone
    if log_in_new:
        new_log_text += textone
    if match:
        log_in_orig,log_in_new = switchto
        print "yes, I am actually switching the binning"
        print "so that status is (",log_in_orig,',',log_in_new,")"
    # if there is a second group (if I have a separator), change which bin I'm in, and add to the end of the current line!
    if texttwo is not None:
        orig_log_text,new_log_text,log_in_orig,log_in_new = parse_line(orig_log_text,new_log_text,log_in_orig,log_in_new,texttwo)
    return orig_log_text,new_log_text,log_in_orig,log_in_new
fp = open(sys.argv[1],'r')
text_list = []
print 'opened',sys.argv[1]
log_in_orig = True
log_in_new = True
head_title = None
new_title = None
#{{{ pull out just the part between the document text
orig_textlist = []
new_textlist = []
j = 0
for thisline in fp:
    if j == 0:
        if thisline[:12] == '%ONEWORDDIFF':
            print "found %ONEWORDDIFF marker, title is:"
            head_title = 'HEAD\n'
            new_title = thisline[14:]
            print new_title
            this_is_a_onewordfile = True
        else:
            this_is_a_onewordfile = False
        if this_is_a_onewordfile:
            print "I found this to be a oneword format file"
        else:
            print "I did not find this to be a oneword format file"
    if this_is_a_onewordfile:# this is only stored if it's a onewordfile
        #new processing for oneworddiff
        #{{{ check to see if I have a separator, and set switchto, to show where I switch
        orig_log_text,new_log_text,log_in_orig,log_in_new = parse_line('','',log_in_orig,log_in_new,thisline)
        if len(orig_log_text) > 0:
            orig_textlist.append(orig_log_text)
        if len(new_log_text) > 0:
            new_textlist.append(new_log_text)
    else:
        #standard processing
        if thisline[-11:] == '%FIRSTSET%\n': # if the first set, treat like it's not a comment
            if log_in_orig:
                orig_textlist.append(thisline)
            if log_in_new:
                new_textlist.append(thisline)
        else:
            if (thisline[:7] == '<<<<<<<'):
                log_in_orig = True
                log_in_new = False
                if (head_title is None): # for the first marker, store the title
                        head_title = thisline[7:]
                elif thisline[7:] == head_title:
                    pass
                else:
                    raise ValueError("I don't understand line %d, which seems to give an inconsistent head title.  It gave:\n%s\nvs expected:\n%s"%(j,thisline[7:],head_title))
            elif (thisline[:7] == '>>>>>>>'):
                log_in_orig = True
                log_in_new = True
                if (new_title is None): # for the first marker, store the title
                        new_title = thisline[7:]
                elif thisline[7:] == new_title:
                    pass
                else:
                    raise ValueError("I don't understand line %d, which seems to give an inconsistent new title.  It gave:\n%s\nvs expected:\n%s"%(j,thisline[7:],new_title))
            elif (thisline[:7] == '======='):
                log_in_orig = False
                log_in_new = True
            else:
                if log_in_orig:
                    orig_textlist.append(thisline)
                if log_in_new:
                    new_textlist.append(thisline)
    j+=1
if this_is_a_onewordfile:
    print "I found this to be a oneword format file"
else:
    print "I did not find this to be a oneword format file"
fp.close()
#{{{ write out the result
newfile = re.sub(r"(.*)",r'\1.merge_new',sys.argv[1]) 
fp = open(newfile,'w')
new_textlist = ['#%%%%%BRANCH TITLE (This side is saved): '+new_title] + new_textlist
fp.write(''.join(new_textlist))
fp.close()
newfile = re.sub(r"(.*)",r'\1.merge_head',sys.argv[1]) 
fp = open(newfile,'w')
orig_textlist = ['#%%%%%BRANCH TITLE: '+head_title] + orig_textlist
fp.write(''.join(orig_textlist))
fp.close()
#}}}
