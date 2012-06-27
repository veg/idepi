BIBTEX = bibtex -terse
DVIPDF = dvipdf
LATEX  = latex -interaction=batchmode

.PHONY : all clean open

all : plos.pdf

plos.pdf : plos.dvi
	$(DVIPDF) $<

plos.dvi : plos.tex plos.aux plos.bbl
	$(LATEX) $(firstword $<)
	$(LATEX) $(firstword $<)

plos.bbl : plos.aux
	-$(BIBTEX) $<

plos.aux : plos.tex
	$(LATEX) $<

clean :
	-rm plos.aux plos.bbl plos.blg plos.dvi plos.log plos.pdf

open : plos.pdf
	-open $<
