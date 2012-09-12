BIBTEX = bibtex -terse
TEX = pdflatex -file-line-error -interaction=batchmode -shell-escape

BIBFILE = plos
# FIGUREFILES =
LATEXFILE = plos

.INTERMEDIATE : $(LATEXFILE).aux $(LATEXFILE).bbl $(LATEXFILE).blg $(LATEXFILE).log
.PHONY : all clean view

all : $(LATEXFILE).pdf

%.pdf %.log : %.tex %.aux %.bbl $(FIGUREFILES)
	$(TEX) $<
	$(TEX) $<

%.bbl %.blg : %.aux $(BIBFILE).bib
	-$(BIBTEX) $<

%.aux %.log : %.tex $(FIGUREFILES)
	$(TEX) $<

clean :
	@rm -f $(LATEXFILE).aux $(LATEXFILE).bbl $(LATEXFILE).blg $(LATEXFILE).log $(LATEXFILE).pdf

view : $(LATEXFILE).pdf
	-open $<
