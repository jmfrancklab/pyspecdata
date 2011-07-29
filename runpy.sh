#echo "from matlablike import *" > runtest.py
#echo "def lplot(*args):" > runtest.py
#echo "   print 'not clearing the plot'" >> runtest.py
#echo "   return" >> runtest.py
#echo "def lplotfigures(*args):" >> runtest.py
#echo "   print 'not clearing the plot'" >> runtest.py
#echo "   return" >> runtest.py
echo "from matlablike import *" > runtest.py
echo "SHOWPLOTS = True" >> runtest.py
#cat notebook.$1.py.old | sed "s/\<lplot(.*//" | sed "s/##/#/" >> runtest.py
cat notebook.$1.py.old | sed "s/##/#/" | sed "s/.*\<lplot(.*/#&/" | sed "s/\(\<figlistl(\))/\1env='test')/" >> runtest.py
#cat notebook.$1.py.old | sed "s/##/#/" >> runtest.py
echo "show()" >> runtest.py
python runtest.py
