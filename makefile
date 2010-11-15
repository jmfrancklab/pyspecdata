ifndef BASEFILE
  BASEFILE=diff_wrapper.tex
endif
ifndef DATE
  ifdef VER
    VERSTRING=$(VER)
    DATE=:\\\\begin{verbatim}\\n
    DATE+=$(shell svn log -r$(VERSTRING) $(FILE).tex | sed s/-*//)\\n
    DATE+=\\\\end{verbatim}\\n
  else
    VERSTRING=HEAD
    DATE=HEAD
  endif
else
  VERSTRING=\{$(DATE)\}
endif
BIBFILES=
ifndef NOBIB
  BIBFILES+=$(FILE)_onlyforplain.bbl
endif

all: compilation.pdf progress.pdf
clean:
	rm notebook.aux notebook.log notebook.out notebook.*.py.out notebook.*.py.old notebook.*.py.err
compclean:
	rm compilation.aux compilation.log compilation.out compilation.*.py.out compilation.*.py.old compilation.*.py.err
notebook.aux: notebook.tex papers.tex mynotebook.sty 09*.tex notebook*.tex summary.tex inprocess/*.tex library.bib
	pdflatex --shell-escape notebook 
notebook.out: notebook.tex papers.tex mynotebook.sty 09*.tex notebook*.tex summary.tex inprocess/*.tex
	pdflatex --shell-escape notebook 
notebook.bbl: notebook.aux library.bib
	bibtex notebook
notebook_wc.pdf: notebook.tex notebook.out
#papers.tex notebook*.tex mynotebook.sty load_cpmg_emax.py notebook.out summary.tex inprocess.tex summary.tex inprocess/*.tex
	pdflatex --shell-escape notebook 
	-rm notebook.*.py
	mv notebook.pdf notebook_wc.pdf
	bash remove_empty.sh
	aplay beep-14_soft.wav
compilation_wc.pdf: compilation.aux compilation.tex papers.tex mynotebook.sty 09*.tex notebook*.tex notebook.*.py.old
	pdflatex --shell-escape compilation.tex 
	-rm compilation.*.py
	-rm compilation_wc.pdf
	mv compilation.pdf compilation_wc.pdf
	bash remove_empty.sh
dnp_protocol.pdf: dnp_protocol.aux dnp_protocol.tex mynotebook.sty notebook.*.py.old 
	pdflatex --shell-escape dnp_protocol.tex 
	-rm dnp_protocol.*.py
compilation.aux: compilation.tex papers.tex mynotebook.sty 09*.tex notebook*.tex
	pdflatex --shell-escape compilation.tex 
progress.pdf: papers.tex progress.tex mynotebook.sty
	pdflatex --shell-escape progress.tex 

$(FILE)_onlyforplain.tex: $(FILE).tex $(BASEFILE)
	echo "BUILDING plain tex"
	cat $(BASEFILE) | sed "s/INPUTHERE/$(FILE)/" | sed "s/DATEHERE/$(DATE)/" | sed "s/PUTTITLEHERE/$(TITLE)/" > temp.tex
 ifndef NOBIB
	cat temp.tex | sed "s/BIBLIOGRAPHYHERE/\\\\bibliography{library.bib}/" > $(FILE)_onlyforplain.tex
else
	cat temp.tex | sed "s/BIBLIOGRAPHYHERE//" > $(FILE)_onlyforplain.tex
endif
$(FILE)_onlyforplain.aux: $(FILE)_onlyforplain.tex
	echo "BUILDING aux"
	pdflatex --shell-escape $(FILE)_onlyforplain 
$(FILE)_onlyforplain.out: $(FILE)_onlyforplain.tex
	echo "BUILDING out"
	pdflatex --shell-escape $(FILE)_onlyforplain 
$(FILE)_onlyforplain.bbl: library.bib $(FILE)_onlyforplain.tex
	pdflatex --shell-escape $(FILE)_onlyforplain 
 ifndef NOBIB
	echo "BUILDING bibtex"
	bibtex $(FILE)_onlyforplain
	makeindex $(FILE)_onlyforplain
 endif
$(FILE)_plain.pdf: $(BIBFILES) $(FILE)_onlyforplain.out $(FILE)_onlyforplain.aux
	echo "BUILDING pdf"
	pdflatex --shell-escape $(FILE)_onlyforplain 
	pdflatex --shell-escape $(FILE)_onlyforplain 
	echo "MOVING"
	echo "verstring is $(VERSTRING)"
	-rm $(FILE)_onlyforplain.*.py
	mv $(FILE)_onlyforplain.pdf $(FILE)_plain.pdf
	aplay beep-14_soft.wav
$(FILE)_onlyfordiff.tex: $(FILE).tex $(BASEFILE)
	echo "BUILDING diff tex"
	-rm $(FILE)-diff*.tex
	echo "verstring" $(VERSTRING)
	latexdiff-svn -r$(VERSTRING) --exclude-safecmd=o $(FILE).tex
	cat $(BASEFILE) | sed s/INPUTHERE/$(FILE)-diff$(VERSTRING)/ | sed "s/DATEHERE/$(DATE)/" > temp.tex
 ifndef NOBIB
	cat temp.tex | sed "s/BIBLIOGRAPHYHERE/\\\\bibliography{library.bib}/" > $(FILE)_onlyfordiff.tex
else
	cat temp.tex | sed "s/BIBLIOGRAPHYHERE//" > $(FILE)_onlyfordiff.tex
endif
$(FILE)_onlyfordiff.bbl: library.bib $(FILE)_onlyfordiff.tex
	echo "BUILDING bibtex"
	pdflatex --shell-escape $(FILE)_onlyfordiff 
	bibtex $(FILE)_onlyfordiff
	makeindex $(FILE)_onlyfordiff
$(FILE)_onlyfordiff.aux: $(FILE)_onlyfordiff.tex
	echo "BUILDING aux"
	pdflatex --shell-escape $(FILE)_onlyfordiff 
$(FILE)_onlyfordiff.out: $(FILE)_onlyfordiff.tex
	echo "BUILDING out"
	pdflatex --shell-escape $(FILE)_onlyfordiff 
$(FILE)_diff.pdf: $(BIBFILES) $(FILE)_onlyfordiff.out $(FILE)_onlyfordiff.aux
	echo "BUILDING pdf"
	pdflatex --shell-escape $(FILE)_onlyfordiff 
	pdflatex --shell-escape $(FILE)_onlyfordiff 
	echo "MOVING"
	echo "verstring is $(VERSTRING)"
	-rm $(FILE)_onlyfordiff.*.py
	mv $(FILE)_onlyfordiff.pdf $(FILE)_diff.pdf
	aplay beep-14_soft.wav
diff: $(FILE)_diff.pdf
 ifndef NOBIB
   override BIBFILES+=$(FILE)_onlyfordiff.bbl
 endif
plain: $(FILE)_plain.pdf
