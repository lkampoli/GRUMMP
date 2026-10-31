.SUFFIXES: .lyx .txt .ps .tex .html

.SILENT:

.lyx.txt: $*.lyx
	lyx -e text $*.lyx

REMOTE_HTML_DIR=cfog@tetra.mech.ubc.ca:/var/www/GRUMMP

LYX_FILES=UserGuide.lyx Changes.lyx GettingStarted.lyx License.lyx Algorithms.lyx ReadMe.lyx

TEX_FILES=$(LYX_FILES:.lyx=.tex)

TXT_FILES=$(LYX_FILES:.lyx=.txt)

distrib: $(TEX_FILES) UserGuide.ps UserGuide.pdf $(TXT_FILES) UserGuide
	cp Changes.txt ../CHANGES
	cp GettingStarted.txt ../GETTING_STARTED
	cp License.txt ../LICENSE
	cp ReadMe.txt ../README
	cp GRUMMP.lsm UserGuide

$(TEX_FILES):
	lyx -e latex UserGuide.lyx

UserGuide.dvi: $(TEX_FILES)
	# Be sure to run LaTeX sufficiently to get all the biblio
	# information right.  
	latex UserGuide.tex
	bibtex UserGuide
	latex UserGuide.tex
	latex UserGuide.tex

UserGuide.ps: UserGuide.dvi
	dvips -t letter -o UserGuide.ps UserGuide.dvi

UserGuide.pdf: UserGuide.ps
	ps2pdf UserGuide.ps

UserGuide: UserGuide.dvi
	rm -rf UserGuide
	latex2html -local_icons -long_titles 3 -split 4 UserGuide.tex
	rm -f UserGuide/index.html
	rm -rf UserGuide/l2h*

publish: UserGuide.ps UserGuide.pdf UserGuide
	gzip -9 -v UserGuide.ps
	scp UserGuide.pdf UserGuide.ps.gz UserGuide/* ../KNOWN_BUGS $(REMOTE_HTML_DIR)
	scp home_page.html $(REMOTE_HTML_DIR)/index.html

latex_junk = *.aux *.bbl *.blg *.dvi *.lof *.log *.toc 

clean:
	-rm -f *~ $(latex_junk)

distclean: clean
	-rm -f \#* Changes.txt GettingStarted.txt License.txt ReadMe.txt \
		Changes.tex GettingStarted.tex License.tex ReadMe.tex \
		Algorithms.tex Algorithms.txt UserGuide.tex UserGuide.txt \
		UserGuide.pdf UserGuide.ps*
	-rm -rf UserGuide

init:
