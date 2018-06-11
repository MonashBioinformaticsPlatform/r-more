# This file is generated from the corresponding .Rmd file



# Don't run this if you are using our biotraining server, the packages are already installed!
#
# Install the entire Tidyverse collection of packages with:
#
#   install.packages("tidyverse")
#
# or just install packages needed for this section with:
#
#   install.packages(c(
#       "readr",    # read tabular data
#       "dplyr"     # general data frame manipulation
#   ))


# =========
# Functions
# =========

FUNCTION_NAME <- function(ARGUMENT_NAME1, ARGUMENT_NAME2, ...) {
    FUNCTION_BODY
}


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


# debugging a function
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
# 
#
# =========
# For-loops
# =========

for(VARIABLE_NAME in VECTOR) {
  FOR_LOOP_BODY
}


i <- 10
cat("i is",i,"\n")
cat("i squared is",i*i,"\n")
i <- 20
cat("i is",i,"\n")
cat("i squared is",i*i,"\n")
i <- 30
cat("i is",i,"\n")
cat("i squared is",i*i,"\n")
i <- 40
cat("i is",i,"\n")
cat("i squared is",i*i,"\n")
i <- 50
cat("i is",i,"\n")
cat("i squared is",i*i,"\n")


for(i in c(10,20,30,40,50)) {       #    1
    cat("i is",i,"\n")              #      2   4   6   8   10
    cat("i squared is",i*i,"\n")    #        3   5   7   9    11
}                                   #
                                    #
cat("done\n")                       #                            12
                                    #   --order-of-execution-->
                                    #
                                    # i= - a a b b c c d d e  e  e ...


## ---------
## Challenge
## ---------
# 
# 1. What do you think this R code will do? Read the code and make a
# guess, then try it in R. Try adding `cat`s or `print`s to the loop
# body to check what is going on.
# 

myvec <- c(10,20,30,40)

total <- 0
for(item in myvec) {
    total <- total + item
}

total

# 
# 2. What do you think this R code will do? How could it be changed to
# work with any length of `myvec`?
# 

myvec <- c(10,20,30,40)

for(index in 1:4) {
    myvec[index] <- myvec[index] * 2
}

# 
# 3. Write a for-loop to calculate 10 factorial, ie `1*2*3*4* ... *10`.
# 

numbers <- 1:10

... your code here ...

# 
# 
# 
#
# ================================
# A practical programming exercise
# ================================

## -------------------------
## Running external software
## -------------------------

uptime


system("uptime")


system("fastqc --extract --outdir . r-more-files/Day0.fastq")


## ----------------
## Using a for-loop
## ----------------

# construct a command to run
day <- 0
command <- paste0("fastqc --extract --outdir . r-more-files/Day", day, ".fastq")
command


# run fastqc on each of the files

days <- c(0,4,7,10,15,20)

for(day in days) {
    command <- paste0("fastqc --extract --outdir . r-more-files/Day", day, ".fastq")
    cat("Running the command:", command, "\n")
    system(command)
}


## -----------------------------
## Loading the summary.txt files
## -----------------------------

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


sumtabs
class(sumtabs)
length(sumtabs)
str(sumtabs)


sumtabs[[1]]
sumtabs[[2]]


library(dplyr)
bigtab <- bind_rows(sumtabs)

bigtab
table(bigtab$test, bigtab$grade)


# =============
# If-statements
# =============

if (LOGICAL_VALUE) {
    THING_TO_DO_IF_TRUE
} else {
    THING_TO_DO_IF_FALSE
}


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
# into files?
# 
# What about managing data files?
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
library(devtools)
create("mypack")

# ... Edit mylibrary/DESCRIPTION file
# ... Write .R files in mylibrary/R folder

# Load package. Use this during development.
load_all("mypack")

# Build package, including converting inline documentation to .Rd files using roxygen2.
# Check for common problems and missing documentation.
# A CRAN package must pass all checks.
check("mypack")


# To install from GitHub:
devtools::install_github("myusername/mypack")


sessionInfo()

