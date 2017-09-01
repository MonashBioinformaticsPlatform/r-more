
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

# Create unevaluated versions (compact teacher's notes)
UNEVALS=topics/programming_uneval.html \
        topics/tidyverse_uneval.html \
        topics/sequences_and_features_uneval.html

all : $(RS) $(HTMLS) $(UNEVALS) r-more-files.zip

%.html : %.Rmd
	Rscript -e 'rmarkdown::render("$<", "all")'

%_uneval.html : %.Rmd Makefile
	python unevalify.py <$< >topics/temp.Rmd
	Rscript -e 'rmarkdown::render("topics/temp.Rmd", "all")'
	mv topics/temp.html $@
	rm topics/temp.Rmd

r-more-files/%.R : topics/%.Rmd purify.py
	python purify.py <$< >$@

r-more-files.zip : r-more-files/* r-more-files/fastqc-output/* $(RS)
	zip -FSr r-more-files.zip r-more-files

clean :
	rm -f $(HTMLS) $(RS) $(UNEVALS) r-more-files.zip
	rm -rf topics/sequences_and_features_cache
	rm -rf topics/programming_cache
