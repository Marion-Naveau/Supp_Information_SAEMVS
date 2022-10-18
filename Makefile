README.md: README.Rmd
	Rscript -e 'rmarkdown::render("README.Rmd")'; \
	./gh-md-toc --insert README.md

readme: README.md
