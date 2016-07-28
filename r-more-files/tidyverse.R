# This file is generated from the corresponding .Rmd file



# Don't run this if you are using our biotraining server, the packages are already installed!
install.packages(c(
    # Hadley Wickham packages
    "readr",    # read tabular data
    "tidyr",    # data frame tidying functions
    "dplyr",    # general data frame manipulation
    "ggplot2",  # flexible plotting

    # Viridis color scale
    "viridis",

    # (Advanced!) Tidy linear modelling
    "broom"
))


library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(broom)


vignette()
vignette(package="dplyr")
vignette("introduction", package="dplyr")


# =====
# dplyr
# =====

bigtab <- read_csv("r-more-files/fastqc.csv")
bigtab$grade <- factor(bigtab$grade, c("FAIL","WARN","PASS"))

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

summarize(group_by(scoretab, file), total=sum(score))


summarize(group_by(bigtab, grade), count=n())


# ============
# The pipe %>%
# ============

bigtab %>% group_by(grade) %>% summarize(count=n())


rep("hello", 5)
"hello" %>% rep(5)


# =================
# ggplot2 revisited
# =================

ggplot(bigtab,aes(x=file,y=test,color=grade)) + geom_point()


ggplot(bigtab,aes(x=file,y=test,fill=grade)) + geom_tile()


## --------------------------
## Publication quality images
## --------------------------

plot <- ggplot(bigtab,aes(x=file,y=test,fill=grade)) +
    geom_tile(color="black",size=0.5) +
    labs(x="",y="",fill="") +                 # Remove axis labels
    coord_fixed() +                           # Square tiles
    theme_minimal() +                         # Minimal theme, no grey background
    theme(panel.grid=element_blank(),         # No underlying grid lines
          axis.text.x=element_text(           # Vertical text on x axis
              angle=90,vjust=0.5,hjust=1))
plot


ggsave("plot1.png", plot, width=5,  height=5,  dpi=600)
ggsave("plot2.png", plot, width=10, height=10, dpi=300)


# ========================
# A simple RNA-Seq example
# ========================

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

ggplot(counts, aes(x=sample, y=count)) + geom_boxplot() + coord_flip()


ggplot(counts, aes(x=sample, y=log2(count))) + geom_boxplot() + coord_flip()

ggplot(counts, aes(x=log2(count), group=sample)) + geom_line(stat="density")


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

ggplot(counts_norm, aes(x=sample, y=log_norm_count)) + geom_boxplot() + coord_flip()


## -------------
## Visualization
## -------------

ggplot(counts_norm, aes(x=time, y=gene, fill=log_norm_count)) +
    geom_tile() + facet_grid(~ strain) +
    scale_fill_viridis() + theme_minimal()


ggplot(counts_norm, aes(x=time, y=strain, fill=log_norm_count)) +
    geom_tile() + facet_wrap(~ gene) +
    scale_fill_viridis() + theme_minimal()


ggplot(counts_norm, aes(x=time, y=log_norm_count, color=strain, group=strain)) +
    geom_line() + facet_wrap(~ gene, scale="free")


## ---------
## Challenge
## ---------
# 
# 1. Make a line plot as above but just showing the gene HHF1.
# 
# Hint: intToUtf8(utf8ToInt("vtf!gjmufs")-1)
# 
# 2. Which are the three most variable genes?
# 
# Hint:
# intToUtf8(utf8ToInt("xvh#jurxsbe|/#vxppdul}h/#vg/#dqg#duudqjh")-3)
# 
# 3. Different genes have different average expression levels, but what
# we are interested in is how they change over time. Further normalize
# the data by subtracting the average for each gene from
# `log_norm_count`.
# 
# 
#
# ================
# Linear modelling
# ================

sut476 <- counts_norm %>% filter(gene=="SUT476")
sut476_wt <- sut476 %>% filter(strain=="WT")

ggplot(sut476_wt,aes(x=time,y=log_norm_count)) +
    geom_point() +
    geom_smooth(method="lm")


model <- lm(log_norm_count ~ time, data=sut476_wt)
model
# Note: The p-values listed by summary aren't always meaningful.
# It makes no sense to remove the intercept term but retain the time term.
# Significance testing can better be done with anova() and the glht package.
summary(model)
model.matrix(log_norm_count ~ time, data=sut476_wt)


null_model <- lm(log_norm_count ~ 1, data=sut476_wt)
anova(null_model, model)


augment(model, sut476_wt) %>% head

# augment() can also be used to produce predictions for new inputs
augment(model, newdata=sut476_wt[1:5,])

# Let's look at the model fit, and how influential each observation was
augment(model, sut476_wt) %>%
ggplot(aes(x=time, y=log_norm_count)) +
    geom_point(aes(size=.cooksd), alpha=0.5) +
    geom_line(aes(y=.fitted))


model2 <- lm(log_norm_count ~ strain*poly(time,3), data=sut476)
model2

augment(model2, sut476) %>%
ggplot(aes(x=time, y=log_norm_count, color=strain, group=strain)) +
    geom_point(aes(size=.cooksd), alpha=0.5) +
    geom_line(aes(y=.fitted))


sessionInfo()

