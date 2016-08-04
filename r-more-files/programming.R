# This file is generated from the corresponding .Rmd file



# Don't run this if you are using our biotraining server, the packages are already installed!
install.packages(c(
    "readr",    # read tabular data
    "dplyr"     # general data frame manipulation
))


# =========
# Functions
# =========

fahr_to_kelvin <- function(temp) {
    (temp-32) * (5/9) + 273.15
}


# freezing point of water
fahr_to_kelvin(32)
# boiling point of water
fahr_to_kelvin(212)


## ----------
## Variations
## ----------

fahr_to_kelvin <- function(temp) {
    kelvin <- (temp-32) * (5/9) + 273.15
    kelvin
}


fahr_to_kelvin <- function(temp) (temp-32) * (5/9) + 273.15


fahr_to_kelvin <- function(temp) {
    kelvin <- (temp-32) * (5/9) + 273.15
    return(kelvin)
    plot(1:10)
}


fahr_to_kelvin <- function(temp) {
    kelvin <-
        (temp-32) * (5/9) +
        273.15
    return(kelvin)
}


   fahr_to_kelvin<-
function(
                   temp){(temp
  -32  )*(    5
            /9 )+

273.15}


fahr_to_kelvin_broken <- function(temp) {
    (temp-32) * (5/9)
        + 273.15
}

fahr_to_kelvin_broken(212)


charmander <- function(bulbasaur) {
    mew <- (bulbasaur-32) * (5/9) + 273.15
    mew
}


## -------------------
## Composing Functions
## -------------------

kelvin_to_celsius <- function(temp) {
    temp - 273.15
}

#absolute zero in Celsius
kelvin_to_celsius(0)


fahr_to_celsius <- function(temp) {
    temp_k <- fahr_to_kelvin(temp)
    temp_c <- kelvin_to_celsius(temp_k)
    temp_c
}

# freezing point of water in Celsius
fahr_to_celsius(32.0)


debugonce(fahr_to_celsius)
fahr_to_celsius(212)


## ---------
## Challenge
## ---------
# 
# Write a function to calculate the length of the hypotenuse of a right
# angled triangle using Pythagorus's rule, given the lengths of the
# other sides.
# 
# Hint: `sqrt` calculates the square root of a number.
# 
# Testing your code is important. Invent a test case for your code
# consisting of:
# 
# * The input arguments to your function.
# * The return value you expect.
# 
# Confirm that your function works as expected.
# 
# 
#
# ================================
# A practical programming exercise
# ================================

## -------------------------
## Running external software
## -------------------------

system("uptime")


system("fastqc --extract --outdir . r-more-files/Day0.fastq")


## -------------------
## For-loops revisited
## -------------------

for(i in c("a","b","c","d","e")) {  #    1
    cat("i is",i,"\n")              #      2   4   6   8   10
    cat("potato\n")                 #        3   5   7   9    11
}                                   #
cat("done\n")                       #                            12
                                    #   --order-of-execution-->
                                    #
                                    # i= - a a b b c c d d e  e  e ...


day <- 0
command <- paste0("fastqc --extract --outdir . r-more-files/Day", day, ".fastq")
command


days <- c(0,4,7,10,15,20)

for(day in days) {
    command <- paste0("fastqc --extract --outdir . r-more-files/Day", day, ".fastq")
    cat("Running the command:", command, "\n")
    system(command)
}


## -----------------------------
## Loading the summary.txt files
## -----------------------------

install.packages("readr")


library(readr)


read_tsv("Day0_fastqc/summary.txt")


read_tsv("Day0_fastqc/summary.txt", col_names=FALSE)


filename <- "Day0_fastqc/summary.txt"
sumtab <- read_tsv(filename, col_names=FALSE)
colnames(sumtab) <- c("grade", "test", "file")
sumtab$grade <- factor(sumtab$grade, c("FAIL","WARN","PASS"))
sumtab


load_fastqc <- function(filename) {
    sumtab <- read_tsv(filename, col_names=FALSE)
    colnames(sumtab) <- c("grade", "test", "file")
    sumtab$grade <- factor(sumtab$grade, c("FAIL","WARN","PASS"))
    sumtab
}

load_fastqc("Day0_fastqc/summary.txt")


## ---------------------
## Applying the function
## ---------------------

days <- c(0,4,7,10,15,20)
filenames <- paste0("Day", days, "_fastqc/summary.txt")
filenames


sumtabs <- lapply(filenames, load_fastqc)


# Base R way
bigtab <- do.call(rbind, sumtabs)

# The dplyr way
library(dplyr)
bigtab <- bind_rows(sumtabs)

bigtab
table(bigtab$test, bigtab$grade)


# =============
# If statements
# =============

num <- 37                   # 1
if (num > 100) {            #   2
  cat("greater\n")          #
} else {                    #
  cat("not greater\n")      #     3
}                           #
cat("done\n")               #       4
                            # --time-->


num <- 53                            # 1
if (num > 100) {                     #   2
  cat("num is greater than 100\n")   #
}                                    #
cat("done\n")                        #     3
                                     # --time-->


sign <- function(num) {
    if (num > 0) {            # line 1
        return(1)             # line 2
    } else if (num == 0) {    # line 3
        return(0)             # line 4
    } else {
        return(-1)            # line 5
    }
}

sign(-3)
sign(0)
sign(2/3)


## ----
## Quiz
## ----
# 
# Which lines of the function `sign` executed when it was called above,
# and in what order?
# 
# 
#
## ---------------------
## Improving load_fastqc
## ---------------------

load_fastqc("nosuchfile.txt")


load_fastqc <- function(filename) {
    # Check arguments are sane
    if (!file.exists(filename)) {
        warning("No such file: ", filename)
        return(NULL)
    }

    # Load and tidy data
    sumtab <- read_tsv(filename, col_names=FALSE)
    colnames(sumtab) <- c("grade", "test", "file")
    sumtab$grade <- factor(sumtab$grade, c("FAIL","WARN","PASS"))
    sumtab
}

load_fastqc("nosuchfile.txt")


# =================
# Sourcing .R files
# =================

# fastqc.R file should contain:

library(readr)

load_fastqc <- function(filename) {
    sumtab <- read_tsv(filename, col_names=FALSE)
    colnames(sumtab) <- c("grade", "test", "file")
    sumtab$grade <- factor(sumtab$grade, c("FAIL","WARN","PASS"))
    sumtab
}


# From the console:

source("fastqc.R")


## ----------
## Discussion
## ----------
# 
# What other R code from this lesson could we put in a .R file? Or a
# .Rmd file?
# 
# How should we break up a large project (paper/thesis/software package)
# into files, for our own workflow?
# 
# What about managing data files?
# 
# How do we ensure the code for a project can be run by others, ensuring
# it is reproducible?
# 
# How should we share a project with others?
# 
# * [Software Carpentry's list of best practices in
# R](http://swcarpentry.github.io/r-novice-inflammation/06-best-
# practices-R/)
# 
# 
#
# ========
# Packages
# ========

# Create an empty package template
devtools::create("mylibrary")

# ... Edit mylibrary/DESCRIPTION file
# ... Write .R files in mylibrary/R folder

# Load package. Use this during development.
devtools::load_all("mylibrary")

# Build package, including converting inline documentation to .Rds files using roxygen2.
# Check for common problems and missing documentation.
# A CRAN package must pass all checks.
devtools::check("mylibrary")


# To install from GitHub:
devtools::install_github("myusername/mylibrary")


sessionInfo()

