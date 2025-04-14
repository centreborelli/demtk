PDF  = $(shell ls *.tex | sed 's/\.tex/.pdf/')

pdf: $(PDF)

%.pdf : %.tex
	grep ^%SCRIPT $< | cut -c9- | bash
	cat $< | ./.pdflatexfilter > $@

p:;cat p.tex | ./.pdflatexfilter > p.pdf

clean :
	$(RM) $(PDF) f/*.png f/fuji_p.npy f/fuji_q.npy f/fuji_b.npy
