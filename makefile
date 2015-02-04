ifndef BASEFILE
  BASEFILE=diff_wrapper.tex
endif
ifndef DATE
  ifdef VER
    VERSTRING=$(VER)
  else
    VERSTRING=HEAD
  endif
else
  VERSTRING=\{$(DATE)\}
endif
#DATE=:\\\\begin{verbatim}\\n
#DATE+=$(shell svn log -r$(VERSTRING) $(FILE).tex | sed s/-*//)\\n
DATE=$(shell svn log -r$(VERSTRING) $(FILE).tex | sed 's/-\+/--/g' | sed 's/[_]/\\\\\\\\_/g' | sed 's/\//\\\//g')\\n
FNAME=$(shell echo $(FILE) | sed "s/.*\/\(.*\)/\1/")
#DATE+=\\\\end{verbatim}\\n
PLAINBIBFILES=
DIFFBIBFILES=
ifndef NOBIB
  DIFFBIBFILES+=$(FNAME)_onlyfordiff.bbl
  PLAINBIBFILES+=$(FNAME)_onlyforplain.bbl
endif

all: compilation.pdf progress.pdf
clean:
	rm notebook.aux notebook.log notebook.out notebook.*.py.out notebook.*.py.old notebook.*.py.err
compclean:
	rm compilation.aux compilation.log compilation.out compilation.*.py.out compilation.*.py.old compilation.*.py.err
notebook.aux: notebook.tex papers.tex mynotebook.sty 09*.tex notebook*.tex summary.tex inprocess/*.tex library.bib lists.tex
	pdflatex -synctex=1 --shell-escape notebook 
notebook.out: notebook.tex papers.tex mynotebook.sty 09*.tex notebook*.tex summary.tex inprocess/*.tex lists.tex
	pdflatex -synctex=1 --shell-escape notebook 
notebook.bbl: notebook.aux library.bib
	bibtex notebook
notebook_wc.pdf: notebook.tex notebook.out
#papers.tex notebook*.tex mynotebook.sty load_cpmg_emax.py notebook.out summary.tex inprocess.tex summary.tex inprocess/*.tex
	pdflatex -synctex=1 --shell-escape notebook 
	-rm notebook.*.py
	mv notebook.pdf notebook_wc.pdf
	mv notebook.synctex.gz notebook_wc.synctex.gz
	bash remove_empty.sh
	-rm *-oldtmp-*.tex
	#aplay beep-14_soft.wav
compilation_wc.pdf: compilation.aux compilation.tex papers.tex mynotebook.sty 09*.tex notebook*.tex notebook.*.py.old
	pdflatex -synctex=1 --shell-escape compilation.tex 
	-rm compilation.*.py
	-rm compilation_wc.pdf
	mv compilation.pdf compilation_wc.pdf
	bash remove_empty.sh
dnp_protocol.pdf: dnp_protocol.aux dnp_protocol.tex mynotebook.sty notebook.*.py.old 
	pdflatex -synctex=1 --shell-escape dnp_protocol.tex 
	-rm dnp_protocol.*.py
compilation.aux: compilation.tex papers.tex mynotebook.sty 09*.tex notebook*.tex
	pdflatex -synctex=1 --shell-escape compilation.tex 
progress.pdf: papers.tex progress.tex mynotebook.sty
	pdflatex -synctex=1 --shell-escape progress.tex 

$(FNAME)_onlyforplain.tex: $(FILE).tex $(BASEFILE)
	echo "BUILDING plain tex"
	cat $(BASEFILE) | sed "s|INPUTHERE|$(FILE)|" | sed "s/DATEHERE/$(DATE)/" | sed "s/PUTTITLEHERE/$(TITLE)/" > temp.tex
 ifndef NOBIB
	cat temp.tex | sed "s/BIBLIOGRAPHYHERE/\\\\bibliography{library.bib}/" > $(FNAME)_onlyforplain.tex
else
	cat temp.tex | sed "s/BIBLIOGRAPHYHERE//" > $(FNAME)_onlyforplain.tex
endif
$(FNAME)_onlyforplain.aux: $(FNAME)_onlyforplain.tex
	echo "BUILDING aux"
	pdflatex -synctex=1 --shell-escape $(FNAME)_onlyforplain 
$(FNAME)_onlyforplain.out: $(FNAME)_onlyforplain.tex
	echo "BUILDING out"
	pdflatex -synctex=1 --shell-escape $(FNAME)_onlyforplain 
$(FNAME)_onlyforplain.bbl: library.bib $(FNAME)_onlyforplain.tex
	pdflatex -synctex=1 --shell-escape $(FNAME)_onlyforplain 
 ifndef NOBIB
	echo "BUILDING bibtex for plain"
	bibtex $(FNAME)_onlyforplain
	makeindex $(FNAME)_onlyforplain
 endif
$(FNAME)_plain.pdf: $(PLAINBIBFILES) $(FNAME)_onlyforplain.out $(FNAME)_onlyforplain.aux $(FNAME)_onlyforplain.tex mynotebook.sty
	echo "building $(FILE)_plain.pdf"
	pdflatex -synctex=1 --shell-escape $(FNAME)_onlyforplain 
	pdflatex -synctex=1 --shell-escape $(FNAME)_onlyforplain 
	echo "MOVING"
	echo "verstring is $(VERSTRING)"
	-rm $(FNAME)_onlyforplain.*.py
	mv $(FNAME)_onlyforplain.pdf $(FNAME)_plain.pdf
	-mv $(FNAME)_onlyforplain.synctex.gz $(FNAME)_plain.synctex.gz
	#aplay beep-14_soft.wav
$(FNAME)_onlyfordiff.tex: $(FILE).tex $(BASEFILE)
	echo "BUILDING diff tex"
	-rm $(FILE)-diff*.tex
	echo "verstring" $(VERSTRING)
	if latexdiff-svn -r$(VERSTRING) --exclude-safecmd="o" --append-textcmd="item" $(FILE).tex;then echo "exit status good";else echo "exit status bad";cp "$(FILE).tex" "$(FILE)-diff$(VERSTRING).tex";fi
	cat $(BASEFILE) | sed "s|INPUTHERE|$(FILE)-diff$(VERSTRING)|" | sed "s/DATEHERE/$(DATE)/" > temp.tex
 ifndef NOBIB
	cat temp.tex | sed "s/BIBLIOGRAPHYHERE/\\\\bibliography{library.bib}/" > $(FNAME)_onlyfordiff.tex
else
	cat temp.tex | sed "s/BIBLIOGRAPHYHERE//" > $(FNAME)_onlyfordiff.tex
endif
$(FNAME)_onlyfordiff.bbl: library.bib $(FNAME)_onlyfordiff.tex
	echo "BUILDING bibtex for diff"
	pdflatex -synctex=1 --shell-escape $(FNAME)_onlyfordiff 
	bibtex $(FNAME)_onlyfordiff
	makeindex $(FNAME)_onlyfordiff
$(FNAME)_onlyfordiff.aux: $(FNAME)_onlyfordiff.tex
	echo "BUILDING aux"
	pdflatex -synctex=1 --shell-escape $(FNAME)_onlyfordiff 
$(FNAME)_onlyfordiff.out: $(FNAME)_onlyfordiff.tex
	echo "BUILDING out"
	pdflatex -synctex=1 --shell-escape $(FNAME)_onlyfordiff 
$(FNAME)_diff.pdf: $(DIFFBIBFILES) $(FNAME)_onlyfordiff.out $(FNAME)_onlyfordiff.aux $(FNAME)_onlyfordiff.tex mynotebook.sty
	echo "BUILDING pdf"
	pdflatex -synctex=1 --shell-escape $(FNAME)_onlyfordiff 
	pdflatex -synctex=1 --shell-escape $(FNAME)_onlyfordiff 
	echo "MOVING"
	echo "verstring is $(VERSTRING)"
	-rm $(FNAME)_onlyfordiff.*.py
	mv $(FNAME)_onlyfordiff.pdf $(FNAME)_diff.pdf
	-mv $(FNAME)_onlyfordiff.synctex.gz $(FNAME)_diff.synctex.gz
	#aplay beep-14_soft.wav
diff: $(FNAME)_diff.pdf
 ifndef NOBIB
   override DIFFBIBFILES+=$(FNAME)_onlyfordiff.bbl
 endif
plain: $(FNAME)_plain.pdf
	echo "building plain"
