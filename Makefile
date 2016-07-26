
RMDS=index.Rmd \
     slides/introduction.Rmd \
     topics/Rmarkdown.Rmd \
     topics/programming.Rmd \
     topics/tidyverse.Rmd \
     topics/sequences_and_features.Rmd

HTMLS=$(patsubst %.Rmd,%.html,$(RMDS))

# Create stripped down versions of .Rmd files
RS=topics/programming.R \
   topics/tidyverse.R \
   topics/sequences_and_features.R

all : $(RS) $(HTMLS) r-more-files.zip

%.html : %.Rmd
	Rscript -e 'rmarkdown::render("$<", "all")'

%.R : %.Rmd
	python purify.py <$< >$@

r-more-files.zip : r-more-files/* r-more-files/fastqc-output/*
	zip -FSr r-more-files.zip r-more-files

clean :
	rm $(HTMLS) r-more-files.zip
