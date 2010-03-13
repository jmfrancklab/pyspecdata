all: compilation.pdf progress.pdf
notebook_wc.pdf: notebook.tex papers.tex current.tex mynotebook.sty load_cpmg_emax.py notebook.out
	touch notebook.0.py # to prevent it from breaking on remove --> has to be a better way of doing this
	pdflatex --shell-escape notebook.tex
	rm notebook.*.py
	mv notebook.pdf notebook_wc.pdf
notebook.aux: notebook.tex papers.tex mynotebook.sty 09*.tex notebook*.tex
	bibtex notebook
	pdflatex --shell-escape notebook.tex
notebook.out: notebook.tex papers.tex mynotebook.sty 09*.tex notebook*.tex
	pdflatex --shell-escape notebook.tex
compilation_wc.pdf: compilation.aux compilation.tex papers.tex mynotebook.sty 09*.tex notebook*.tex
	pdflatex --shell-escape compilation.tex
	rm compilation.*.py
	rm compilation_wc.pdf
	mv compilation.pdf compilation_wc.pdf
compilation.aux: compilation.tex papers.tex mynotebook.sty 09*.tex notebook*.tex
	pdflatex --shell-escape compilation.tex
progress.pdf: papers.tex progress.tex mynotebook.sty
	pdflatex --shell-escape progress.tex
