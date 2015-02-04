#!/bin/bash

# this is a script that uses vim to resolve conflicts, typically, put it in your home directory
# then, to resolve conflicts, you can run ~/resolve_conflicts.sh,
# which lists the conflicts without a file name, and with a file name, will allow you to use gvim diff to resolve the two conflicts into one version (this creates two temporary files, and compares them to each other)
# once the resolution is done, you run ~/resolve_markdone.sh (also in this directory), which gets rid of the temporary files and marks the conflict as resolved 

# diff is called by git with 7 parameters:
# path old-file old-hex old-mode new-file new-hex new-mode

if [ "$#" == "0" ]; then
	git status -s | sed "/^[^U]/d"
else
	rm $1.merge_head
    cat $1.merge_new | sed "/^#%%%%%BRANCH TITLE/d" > $1
    rm $1.merge_new
    git add $1
fi
