#/bin/bash
rm $1.*.py.old $1.*.py.err
pdflatex --shell-escape $1
pdflatex --shell-escape $1
