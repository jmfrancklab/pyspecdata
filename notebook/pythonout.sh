echo $3 > pythonjobname.txt
# called as "bash pythonout.sh notebook.24.py pythonout.out 24"
#mv temp.py $1
# replace a "latex style comments" w/ obs
# should put a while loop around first sed, to catch all instances
cat $1 | sed "s/^\([ \t]*\)%\(.*[^\\]\)<\(.*[^\\]\)>/\1%\2\'\'\',\3,r\'\'\' /" | sed "s/^\([ \t]*\)%\(.*\)/\1obs(r\'\'\'\2 \'\'\')/" > temp.py
mv temp.py $1
#cat $1 | sed "s/^\([ \t]*\)%/\1#/" | sed '/^[#%]DIFDELCMD/d' | sed -e N -e 's/%DIFDELCMD.*\n//' > temp.py
cat $1 | sed "s/^\([ \t]*\)%/\1#/" | sed '/^[#%]DIFDELCMD/d' | sed -e :a -e N -e "s/\n/__LINEBREAK__/g" -e 's/%DIFDELCMD.*__LINEBREAK__//g' -e ta | sed "s/__LINEBREAK__/\n/g" > temp.py
mv temp.py $1
cat $1 | sed "s/\/mnt\/windows\/linuxfiles/\/home\/franck\/data/g" > temp.py
#cat $1 | sed "s/\/home\/franck\/data/\/mnt\/windows\/linuxfiles/g" > temp.py
mv temp.py $1
#rm temp.py
if [ -e "$1.old" ] && [ -e "$1" ]; then
	if ! diff $1.old $1; then
		# if the old version and the new version are different, then run python, otherwise, do nothing
		python -W "ignore" $1 > $1.out 2> $1.err
		mv $1 $1.old
	fi
else
	# also need to run if both files don't exist
	echo "$1 either old or new doesn't exist... creating"
	python -W "ignore" $1 > $1.out 2> $1.err
	mv $1 $1.old
fi
cat $1.out > $2
if [ -s "$1.err" ] # if the error file isn't empty
then
	#echo "{\\small {\\color{red} {\\tt ERRORS---------------------}$1}}\\\\" >> $2
	echo "\\quad\\\\ {\\small {\\color{red} {\\tt ERRORS---------------------} \\fn{$1.old}}}\\\\" >> $2
	echo "\\begin{verbatim}" >> $2
	cat $1.err >> $2
	echo "\\end{verbatim}" >> $2
	echo "{\\small {\\color{red} {\\tt ---------------------------}}}\\\\" >> $2
fi
