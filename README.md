# r-more

Further topics in R.

* [Go to website.](https://monashbioinformaticsplatform.github.io/r-more/)

## Building

To completely rebuild, type:

```
make clean ; make
```

On an older Ubuntu (eg 14), you might need to use the RStudio pandoc:

```
PATH=/usr/lib/rstudio-server/bin/pandoc:$PATH make
```

## Dependencies

Command line software:

* R
* pandoc
* java (Ubuntu: apt-get install default-jre)
* fastqc

R packages:

```
install.packages(c(
    "readr",
    "tidyr",
    "dplyr",
    "ggplot2",
    "viridis",
    "broom",
    "shiny",
    "rmarkdown"
))

source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite(c(
    "Biostrings",
    "GenomicRanges",
    "BSgenome",
    "rtracklayer",
    "motifRG"
))
```

## License

This material is made available under the Creative Commons Attribution 4.0 International License. See index.html for details.
