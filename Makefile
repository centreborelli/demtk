HTML = $(shell ls *.tex | sed 's/\.tex/.html/')
PDF  = $(shell ls *.tex | sed 's/\.tex/.pdf/')

all: $(HTML) $(PDF)

%.html : %.tex .latex_to_html5
	grep ^%SCRIPT $< | cut -c9- | sh
	@./.process_verbatims $<
	./.latex_to_html5 < $< > $@

%.pdf : %.tex %.html .header.tex .footer.tex
	cat .header.tex $< .footer.tex | ./.pdflatexfilter > $@

clean :
	$(RM) $(HTML) $(PDF) canvas *.png *.tiff *.g *.txt o/*.png
