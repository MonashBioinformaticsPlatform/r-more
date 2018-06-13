# This file is generated from the corresponding .Rmd file



# install.packages(c(
#     "tidyverse",
#     "viridis",
#     "broom"
# ))


library(tidyverse) # Load all "tidyverse" libraries.
# OR
# library(readr)   # Read tabular data.
# library(tidyr)   # Data frame tidying functions.
# library(dplyr)   # General data frame manipulation.
# library(ggplot2) # Flexible plotting.

library(viridis)   # Viridis color scale.


vignette()
vignette(package="dplyr")
vignette("dplyr", package="dplyr")


bigtab <- read_csv("r-more-files/fastqc.csv")


# =================
# ggplot2 revisited
# =================

ggplot(bigtab, aes(x=file,y=test,color=grade)) +
    geom_point()


ggplot(bigtab, aes(x=file,y=test,fill=grade)) +
    geom_tile()


## --------------------------
## Publication quality images
## --------------------------

y_order <- sort(unique(bigtab$test), decreasing=T)  # y axis plots from bottom to top, so reverse
bigtab$test <- factor(bigtab$test, levels=y_order)

x_order <- unique(bigtab$file)
bigtab$file <- factor(bigtab$file, levels=x_order)

# Only necessary if not continuing from previous lesson on programming!
color_order <- c("FAIL", "WARN", "PASS")
bigtab$grade <- factor(bigtab$grade, levels=color_order)

myplot <- ggplot(bigtab, aes(x=file, y=test, fill=grade)) +
    geom_tile(color="black", size=0.5) +           # Black border on tiles
    scale_fill_manual(                             # Colors, as color hex codes
        values=c("#ee0000","#ffee00","#00aa00")) +
    labs(x="", y="", fill="") +                    # Remove axis labels
    coord_fixed() +                                # Square tiles
    theme_minimal() +                              # Minimal theme, no grey background
    theme(panel.grid=element_blank(),              # No underlying grid lines
          axis.text.x=element_text(                # Vertical text on x axis
              angle=90,vjust=0.5,hjust=0))
myplot


ggsave("plot1.png", myplot, width=5,  height=5,  dpi=600)
ggsave("plot2.png", myplot, width=10, height=10, dpi=300)


# =====
# dplyr
# =====

# input         +--------+        +--------+        +--------+      result
#  data   %>%   |  verb  |  %>%   |  verb  |  %>%   |  verb  |  ->   data
#   frame       +--------+        +--------+        +--------+        frame


## -------
## tibbles
## -------

bigtab


as.data.frame(bigtab)
View(bigtab)


## ------
## filter
## ------

filter(bigtab, grade == "FAIL")


## -------
## arrange
## -------

arrange(bigtab, grade)

# desc( ) can be used to reverse the sort order
arrange(bigtab, desc(grade))


## ------
## select
## ------

select(bigtab, test,grade)
select(bigtab, 2,1)
select(bigtab, foo=file, bar=test, baz=grade)


select(bigtab, -file)


## -----
## Joins
## -----

fwp <- c("FAIL","WARN","PASS")
scoring <- tibble(grade=factor(fwp,levels=fwp), score=c(0,0.5,1))

# Or:
# scoring <- data.frame(grade=factor(fwp,levels=fwp), score=c(0,0.5,1))

scoring


scoretab <- left_join(bigtab, scoring, by="grade")
scoretab


### ---------
### Challenge
### ---------
# 
# 1. Filter `scoretab` to get only "WARN" or "FAIL" grades.
# 
# 2. Sort the result so that "FAIL"s are at the top.
# 
# 
#
## ------
## mutate
## ------

mutate(scoretab, doublescore = score*2)


scoretab2 <- scoretab
scoretab2$doublescore <- scoretab2$score * 2


## ---------
## summarize
## ---------

summarize(scoretab, total=sum(score))


group_by(scoretab, file)

summarize(group_by(scoretab, file), average_score=mean(score))


summarize(group_by(scoretab, grade), count=n())


# ============
# The pipe %>%
# ============

scoretab %>% group_by(grade) %>% summarize(count=n())


rep(paste("hello", "world"), 5)

"hello" %>% paste("world") %>% rep(5)


# =====
# tidyr
# =====

untidy <- read_csv(
    "country,     male-young, male-old, female-young, female-old
     Australia,            1,        2,            3,          4
     New Zealand,          5,        6,            7,          8")
untidy


gathered <- gather(untidy, key=group, value=cases, -country)
gathered


spread(bigtab, key=file, value=grade)


separate(gathered, col=group, into=c("gender","age"))


tidied <- untidy %>%
    gather(key=group, value=cases, -country) %>%
    separate(group, into=c("gender","age"))


# Advanced
nested <- nest(gathered, -country)
nested
nested$data
unnest(nested)


### ---------
### Challenge
### ---------
# 
# You receive data on a set of points. The points are in two dimensions
# (`dim`), and each point has x and y coordinates. Unfortunately it
# looks like this:
# 

df <- read_csv(
    "dim, A_1, A_2, B_1, B_2, B_3, B_4, B_5
     x,   2,   4,   1,   2,   3,   4,   5
     y,   4,   4,   2,   1,   1,   1,   2")

# 
# 1. Tidy the data by *gathering* all of the columns except `dim`. What
# what does each row now represent?
# 
# 2. We want to plot the points as a scatter-plot, using either `plot`
# or `ggplot`. *Spread* the gathered data so that this is possible. Now
# what do the rows represent?
# 
# 3. What other tidying operation could be applied to this data?
# 
# 
#
# ==================
# An RNA-Seq example
# ==================

## -------
## Tidying
## -------

# Use readr's read_csv function to load the file
counts_untidy <- read_csv("r-more-files/read-counts.csv")

counts <- counts_untidy %>%
    gather(sample, count, -Feature, factor_key=TRUE) %>%
    separate(sample, sep=":", into=c("strain","time"), convert=TRUE, remove=FALSE) %>%
    mutate(
        strain = factor(strain, levels=c("WT","SET1","RRP6","SET1-RRP6")),
        set1 = (strain == "SET1" | strain == "SET1-RRP6"),
        rrp6 = (strain == "RRP6" | strain == "SET1-RRP6")
    ) %>%
    filter(time >= 2) %>%
    select(sample, strain, set1, rrp6, time, gene=Feature, count)

summary(counts)


## --------------------------------
## Transformation and normalization
## --------------------------------

ggplot(counts, aes(x=sample, y=count)) +
    geom_boxplot() +
    coord_flip()


ggplot(counts, aes(x=sample, y=log2(count))) +
    geom_boxplot() +
    coord_flip()

ggplot(counts, aes(x=log2(count), group=sample)) +
    geom_line(stat="density")


normalizer <- counts %>%
    filter(gene == "SRP68") %>%
    select(sample, norm=count)

moderation <- 0.5/mean(normalizer$norm)

counts_norm <- counts %>%
    left_join(normalizer, by="sample") %>%
    mutate(
        norm_count = count/norm,
        log_norm_count = log2(norm_count+moderation)
    )

ggplot(counts_norm, aes(x=sample, y=log_norm_count)) +
    geom_boxplot() +
    coord_flip()


# For a full sized RNA-Seq dataset:
library(edgeR)
mat <- counts_untidy %>% select(-Feature) %>% as.matrix
adjusted_lib_sizes <- colSums(mat) * calcNormFactors(mat)
normalizer_by_tmm <- tibble(sample=names(adjusted_lib_sizes), norm=adjusted_lib_sizes)


## -------------
## Visualization
## -------------

### ---------
### Challenge
### ---------
# 
# 1. Get all the rows in `counts_norm` relating to the histone gene
# "HHT1".
# 
# 2. Plot this data with ggplot2. Use `time` as the x-axis,
# `log_norm_count` as the y-axis, and color the data by `strain`. Try
# using geoms: `geom_point()`, `geom_line()`.
# 

ggplot( ... , aes(x= ... , y= ... , color= ... )) + ...

# 
# Extensions:
# 
# Compare plots of `log_norm_count`, `norm_count`, and `count`.
# 
# Experiment with other geoms and other ways to assign columns to
# aesthetics.
# 
# 
#
### ---------------------------
### Whole dataset visualization
### ---------------------------

ggplot(counts_norm, aes(x=sample, y=gene, fill=log_norm_count)) +
    geom_tile() +
    scale_fill_viridis() +
    theme_minimal() +
    theme(axis.text.x=element_text(           # Vertical text on x axis
              angle=90,vjust=0.5,hjust=1))


ggplot(counts_norm, aes(x=time, y=gene, fill=log_norm_count)) +
    geom_tile() +
    facet_grid(~ strain) +
    scale_fill_viridis() +
    theme_minimal()


ggplot(counts_norm, aes(x=time, y=strain, fill=log_norm_count)) +
    geom_tile() +
    facet_wrap(~ gene) +
    scale_fill_viridis() +
    theme_minimal()


ggplot(counts_norm, aes(x=time, y=log_norm_count, color=strain, group=strain)) +
    geom_line() +
    facet_wrap(~ gene, scale="free")


## ---------
## Exercises
## ---------
# 
# 1. Which are the three most variable genes?
# 
# Hint:
# `intToUtf8(utf8ToInt("xvh#jurxsbe|/#vxppdul}h/#vg/#dqg#duudqjh")-3)`
# 
# 2. Different genes have different average expression levels, but what
# we are interested in is how they change over time. Further normalize
# the data by subtracting the average for each gene from
# `log_norm_count`.
# 
# 
#
# ===============================
# Appendix: Tidy linear modelling
# ===============================

sut476 <- counts_norm %>% filter(gene=="SUT476")
sut476_wt <- sut476 %>% filter(strain=="WT")

ggplot(sut476_wt,aes(x=time,y=log_norm_count)) +
    geom_point() +
    geom_smooth(method="lm")


model <- lm(log_norm_count ~ time, data=sut476_wt)
model
# Note: The p-values listed by summary aren't always meaningful.
# It makes no sense to remove the intercept term but retain the time term.
# Significance testing can better be done with anova() and the multcomp package.
summary(model)
model.matrix(log_norm_count ~ time, data=sut476_wt)


null_model <- lm(log_norm_count ~ 1, data=sut476_wt)
anova(null_model, model)


library(broom)

augment(model, sut476_wt) %>% head

# augment() can also be used to produce predictions for new inputs
augment(model, newdata=sut476_wt[1:5,])

# Let's look at the model fit
augment(model, sut476_wt) %>%
ggplot(aes(x=time)) +
    geom_point(aes(y=log_norm_count)) +
    geom_line(aes(y=.fitted))


model2 <- lm(log_norm_count ~ strain*poly(time,3), data=sut476)
model2

augment(model2, sut476) %>%
ggplot(aes(x=time, color=strain, group=strain)) +
    geom_point(aes(y=log_norm_count), alpha=0.5) +
    geom_line(aes(y=.fitted))


sessionInfo()

