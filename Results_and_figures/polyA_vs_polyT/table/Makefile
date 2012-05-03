NAME=table
TEXSRCS=table.tex
# BIBTEXSRCS= /home/jorgsk/phdproject/bibtex/jorgsk.bib
CLEAN_FILES+= *.log *.aux *.glo *.ilg *.dvi *.toc *.lof *.bbl *.blg *.gls \
              ${NAME}.pdf ${NAME}.ps
# FIGDIR= fig
# FIGS= ${FIGDIR}/*.mp
RM= rm -f

default: ${NAME}.pdf

${NAME}.pdf : ${TEXSRCS}
	pdflatex ${NAME}

# The original makefile from Olaf. Stuff to do with images and bibtex has been
# commented out. Uncomment them to get the function you desire.
#
# ${NAME}.dvi : ${TEXSRCS}  ${BIBTEXSRCS} ${FIGS}
	# # $(MAKE) -C ${FIGDIR}
	# latex ${NAME}
	# -grep -q "LaTeX Warning: There were undefined" ${NAME}.log && bibtex ${NAME}
	# # (grep -q "find a database entry for" ${NAME}.blg && exit 1 || true)
	# # -test -f ${NAME}.glo && \
	  # # makeindex ${NAME}.glo -s nomencl.ist -o ${NAME}.gls; latex ${NAME}
	# grep -q "Label(s) may have changed." ${NAME}.log && latex ${NAME} || true
	# grep -q "Label(s) may have changed." ${NAME}.log && latex ${NAME} || true
	# grep -q "Label(s) may have changed." ${NAME}.log && latex ${NAME} || true

# ${NAME}.ps : ${NAME}.dvi
	# dvips ${NAME}

# ${NAME}.pdf : ${NAME}.ps
	# ps2pdf ${NAME}.ps

.PHONY : clean
clean :
	$(RM) ${CLEAN_FILES}

