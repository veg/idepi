BIBTEX = bibtex -terse
DVIPDF = dvipdf
LATEX = latex -interaction=batchmode

BIBFILE = plos
# FIGUREFILES = 
LATEXFILE = plos

.PHONY : all clean open

all : $(LATEXFILE).pdf

$(LATEXFILE).pdf : $(LATEXFILE).dvi
	$(DVIPDF) $<

$(LATEXFILE).dvi : $(LATEXFILE).tex $(LATEXFILE).aux $(LATEXFILE).bbl
	$(LATEX) $(firstword $<)
	$(LATEX) $(firstword $<)

$(LATEXFILE).bbl : $(LATEXFILE).aux $(BIBFILE).bib
	-$(BIBTEX) $<

$(LATEXFILE).aux : $(LATEXFILE).tex $(FIGUREFILES)
	$(LATEX) $<

clean :
	-rm $(LATEXFILE).aux $(LATEXFILE).bbl $(LATEXFILE).blg $(LATEXFILE).dvi $(LATEXFILE).log $(LATEXFILE).pdf

open : $(LATEXFILE).pdf
	-open $<
