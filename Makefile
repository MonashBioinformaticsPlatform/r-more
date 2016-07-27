
RMDS=index.Rmd \
     slides/introduction.Rmd \
     topics/Rmarkdown.Rmd \
     topics/programming.Rmd \
     topics/tidyverse.Rmd \
     topics/sequences_and_features.Rmd

HTMLS=$(patsubst %.Rmd,%.html,$(RMDS))

# Create stripped down versions of .Rmd files
RS=r-more-files/programming.R \
   r-more-files/tidyverse.R \
   r-more-files/sequences_and_features.R

all : $(RS) $(HTMLS) r-more-files.zip

%.html : %.Rmd
	Rscript -e 'rmarkdown::render("$<", "all")'

r-more-files/%.R : topics/%.Rmd purify.py
	python purify.py <$< >$@

r-more-files.zip : r-more-files/* r-more-files/fastqc-output/* $(RS)
	zip -FSr r-more-files.zip r-more-files

clean :
	rm $(HTMLS) $(RS) r-more-files.zip
