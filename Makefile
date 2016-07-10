
RMDS=$(wildcard *.Rmd */*.Rmd)
HTMLS=$(patsubst %.Rmd,%.html,$(RMDS))

all : $(HTMLS)

%.html : %.Rmd
	Rscript -e 'rmarkdown::render("$<", "all")'

clean :
	rm $(HTMLS)
