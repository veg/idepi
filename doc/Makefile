LATEX = latex -interaction=batchmode
BIBTEX = bibtex -terse

all: plos.pdf

.PHONY: clean

plos.pdf: plos.tex plos.bib
	${LATEX} plos.tex
	-${BIBTEX} plos.aux
	${LATEX} plos.tex
	${LATEX} plos.tex
	dvipdf plos.dvi

clean:
	-rm plos.aux plos.bbl plos.blg plos.dvi plos.log plos.pdf
