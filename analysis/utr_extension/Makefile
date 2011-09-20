all: table.pdf
# all: table.eps

table.eps : table.tex
	latex table.tex 
	dvips table.dvi

table.pdf : table.tex
	pdflatex $< 

table.tex : 3utr_extension.tsf template.tex csv2latex.py
	python csv2latex.py
