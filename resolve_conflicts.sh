#!/bin/bash
if [ "$#" == "0" ]; then
	git status -s | sed "/^[^U]/d"
else
    if ! [ -e $1.merge_head ]; then
        python ~/notebook/split_conflict.py $1
    fi
	/c/Program\ Files/vim72/vim72/gvim -d $1.merge_head $1.merge_new
fi
