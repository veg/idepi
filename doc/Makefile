all: plos.pdf

.PHONY: clean

plos.pdf: plos.tex plos.bib
	latex plos.tex
	-bibtex plos.aux
	latex plos.tex
	latex plos.tex
	dvipdf plos.dvi

clean:
	-rm plos.aux plos.bbl plos.blg plos.dvi plos.log plos.pdf
