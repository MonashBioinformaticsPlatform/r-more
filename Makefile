
RMDS=$(wildcard *.Rmd topics/*.Rmd slides/introduction.Rmd)
HTMLS=$(patsubst %.Rmd,%.html,$(RMDS))

all : $(HTMLS) r-more-files.zip

%.html : %.Rmd
	Rscript -e 'rmarkdown::render("$<", "all")'

r-more-files.zip : r-more-files/*
	zip -FSr r-more-files.zip r-more-files

clean :
	rm $(HTMLS) r-more-files.zip
