NAME=second_version
TEXSRCS=second_version.tex
# BIBTEXSRCS= /home/jorgsk/phdproject/bibtex/jorgsk.bib
CLEAN_FILES+= *.log *.aux *.glo *.ilg *.dvi *.toc *.lof *.bbl *.blg *.gls \
              ${NAME}.pdf ${NAME}.ps
# FIGDIR= fig
# FIGS= ${FIGDIR}/*.mp
RM= rm -f

default: ${NAME}.dvi

# ${NAME}.pdf : ${TEXSRCS}
	# pdflatex ${NAME}

# ${NAME}.dvi : ${TEXSRCS}  ${BIBTEXSRCS} ${FIGS}
${NAME}.dvi : ${TEXSRCS}
	# $(MAKE) -C ${FIGDIR}
	latex ${NAME}
	-grep -q "LaTeX Warning: There were undefined" ${NAME}.log && bibtex ${NAME}
	(grep -q "find a database entry for" ${NAME}.blg && exit 1 || true)
	-test -f ${NAME}.glo && \
	  makeindex ${NAME}.glo -s nomencl.ist -o ${NAME}.gls; latex ${NAME}
	grep -q "Label(s) may have changed." ${NAME}.log && latex ${NAME} || true
	grep -q "Label(s) may have changed." ${NAME}.log && latex ${NAME} || true
	grep -q "Label(s) may have changed." ${NAME}.log && latex ${NAME} || true

${NAME}.ps : ${NAME}.dvi
	dvips ${NAME}

${NAME}.pdf : ${NAME}.ps
	ps2pdf ${NAME}.ps

.PHONY : clean
clean :
	$(RM) ${CLEAN_FILES}

