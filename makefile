all: notebook_wc.pdf
notebook_wc.pdf: notebook.tex mynotebook.sty notebook.out
	touch notebook.0.py # to prevent it from breaking on remove --> has to be a better way of doing this
	pdflatex --shell-escape notebook.tex
	rm notebook.*.py
	mv notebook.pdf notebook_wc.pdf
notebook.aux: notebook.tex mynotebook.sty 09*.tex notebook*.tex
	bibtex notebook
	pdflatex --shell-escape notebook.tex
notebook.out: notebook.tex mynotebook.sty *notebook*.tex
	pdflatex --shell-escape notebook.tex
