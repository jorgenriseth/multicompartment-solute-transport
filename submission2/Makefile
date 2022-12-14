default: 
	pdflatex main.tex

final:
	pdflatex main.tex
	bibtex paper
	pdflatex main.tex
	pdflatex main.tex

clean:
	rm *~ *.aux *.out *.toc *.blg *.bbl *.log *.tdo *.pag

diff:
	cp main.tex paper_new.tex
	#git checkout 3366c4f17f66add725e2f54501a2b2cf9c5bc8d2 -- main.tex
	#cp main.tex paper_old.tex
	#cp paper_new.tex main.tex
	latexdiff paper_old.tex main.tex > diff.tex
	pdflatex diff.tex
	bibtex biblio.aux
	pdflatex diff.tex
	pdflatex diff.tex
