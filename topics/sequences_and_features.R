# This file is generated from the corresponding .Rmd file



options(width=90)


# This loads the Bioconductor installer
source("https://bioconductor.org/biocLite.R")

# Install a basic set of packages
biocLite()

# Install further packages used in this tutorial
biocLite(c(
    "Biostrings",
    "GenomicRanges",
    "BSgenome",
    "rtracklayer",
    "motifRG"
))

# If you want to install further packages in future, you can use
#   library(BiocInstaller)
#   biocLite( ... )


vignette()
vignette(package="Biostrings")
vignette("BiostringsQuickOverview", package="Biostrings")


# ================================
# DNA sequences and genomic ranges
# ================================


library(Biostrings)     # Provides DNAString, DNAStringSet, etc
library(BSgenome)       # Provides getSeq()
library(GenomicRanges)  # Provides GRanges, etc
library(rtracklayer)    # Provides import() and export()

library(dplyr)          # Pipe %>%


## ---------
## DNAString
## ---------


myseq <- DNAString("ACCATTGATTAT")
myseq
class(myseq)
reverseComplement(myseq)
translate(myseq)
subseq(myseq, 3,5)
as.character(myseq)


methods(class="DNAString")


## ------------
## DNAStringSet
## ------------


myset <- DNAStringSet(list(chrI=myseq, chrII=DNAString("ACGTACGT")))

# A DNAStringSet is list-like
myset$chrII
# or myset[["chrII"]]
# or myset[[2]]


## -------
## GRanges
## -------


range1 <- GRanges("chrI", IRanges(start=3,end=5), strand="+")
getSeq(myset, range1)

range2 <- GRanges("chrI", IRanges(start=3,end=5), strand="-")
getSeq(myset, range2)


seqnames(range1)
start(range1)
end(range1)
strand(range1)
as.data.frame(range1)


# GRanges are sometimes like vectors:
c(range1, range2)

# GRanges can have metadata columns, so they are also like data frames:
mcols(range1)$wobble <- 10
range1
mcols(range1)$wobble
range1$wobble

# A handy way to create a GRanges
as("chrI:3-5:+", "GRanges")


## ----------------------
## Challenge {.challenge}
## ----------------------


Reverse complement the following DNA sequence and then translate to an amino acid sequence:


TTCCATTTCCAT




# =============
# Loading files
# =============


## -----------------
## Loading sequences
## -----------------


seqs <- readDNAStringSet("r-more-files/Escherichia_coli_k_12.GCA_000800765.1.29.dna.genome.fa")
seqs

# Our chromosome name is too verbose.
# Remove everything from the name after the first space.
names(seqs)
names(seqs) <- sub(" .*","",names(seqs))
names(seqs)


## ----------------
## Loading features
## ----------------


features <- import("r-more-files/Escherichia_coli_k_12.GCA_000800765.1.29.gtf")

# Optional: just retain the columns of metadata we need
mcols(features) <- mcols(features)[,c("type","gene_name","gene_id")]

features


feat <- features[4,]
feat
feat_seq <- getSeq(seqs, feat)
feat_seq
translate(feat_seq)


subset(features, gene_name == "lacA")
# Equivalently:
#   features[features$gene_name == "lacA" & !is.na(features$gene_name),]


cds <- subset(features, type == "CDS")
cds
# Equivalently:
#   features[features$type == "CDS",]


# =============================
# Further operations on GRanges
# =============================


## -----------
## Intra-range
## -----------


?"intra-range-methods"


feat <- features[4,]
feat_stop <- resize(feat, width(feat)+3)
seq_stop <- getSeq(seqs, feat_stop)
translate(seq_stop)


input
                        (----)
                        .    .
resize                  (--------)
resize(fix="end")   (--------)
                        .    .
flank              (---).    .
flank(start=F)          .    .(---)
                        .    .
promoters          (------)  .
                        .    .
narrow                  .(--).


## -----------
## Inter-range
## -----------


?"inter-range-methods"
?"setops-methods"


query <- as("Chromosome:9500-10000:+", "GRanges")
hits <- findOverlaps(query, features, ignore.strand=TRUE)
hits
subjectHits(hits)
features[subjectHits(hits),]

findOverlaps(query, features, ignore.strand=FALSE)


input
        (--------)
             (----------------)
                     (-----)
                                    (---)

range   (-------------------------------)

reduce  (----------------------)    (---)

disjoin (---)(---)(-)(-----)(--)    (---)

GenomicRanges::setdiff(range(input),input)
                                (--)


## ----------------------
## Challenge {.challenge}
## ----------------------


What are E. coli's most common start and stop codons?

The start codon is the first three bases of the CDS, and the stop codon is the three bases following the end of the CDS.

Hint: Recall that we could get all CDS ranges with:


cds <- subset(features, type == "CDS")


Hint: Use `flank()` and `resize()` to manipulate these ranges.



# =====================
# Finding a known motif
# =====================


vmatchPattern("AGGAGGT", seqs)


vmatchPattern(reverseComplement(DNAString("AGGAGGT")), seqs)


query <- DNAString("AGGAGGT")
max.mismatch <- 1

fwd <- vmatchPattern(query, seqs, max.mismatch=max.mismatch)
fwd <- as(fwd, "GRanges")
strand(fwd) <- "+"
rev <- vmatchPattern(reverseComplement(query), seqs, max.mismatch=max.mismatch)
rev <- as(rev, "GRanges")
strand(rev) <- "-"

complete <- c(fwd, rev)
complete

# Write to GFF file
export(complete, "motif-matches.gff")


# =====================
# De novo motif finding
# =====================


# Note: bacteria do not have introns
# In a eukaryote, you would need to merge CDS by transcript

size <- 20

initiation_regions <- flank(cds, size, start=TRUE)
initiation_seqs <- getSeq(seqs, initiation_regions)
names(initiation_seqs) <- initiation_regions$gene_id

# Look for any composition bias
library(seqLogo)
letter_counts <- consensusMatrix(initiation_seqs)
probs <- prop.table(letter_counts[1:4,], 2)
seqLogo(probs, ic.scale=FALSE)
seqLogo(probs)

# Generate a background set of sequences by shuffling
shuffle <- function(dna) {
    strsplit(as.character(dna),"")[[1]] %>%
        sample %>%
        paste(collapse="") %>%
        DNAString
}

background_seqs <- lapply(initiation_seqs, shuffle) %>% DNAStringSet
names(background_seqs) <- paste0(names(background_seqs), "-shuffled")


library(motifRG)

results <- findMotifFgBg(
    initiation_seqs, background_seqs,
    both.strand=FALSE, start.width=4)

summaryMotif(results$motif, results$category)

motifHtmlTable(results)

refined <- refinePWMMotif(results$motifs[[1]]@match$pattern, initiation_seqs)
seqLogo(refined$model$prob)


writeXStringSet(initiation_seqs, "fg.fa")
writeXStringSet(background_seqs, "bg.fa")


system("meme -dna -maxsize 1000000 fg.fa")


# ==========
# Next steps
# ==========


### ---------------------------
### BSgenome.Hsapiens.UCSC.hg38
### ---------------------------


### ---------------------------------
### TxDb.Hsapiens.UCSC.hg38.knownGene
### ---------------------------------


### ------------
### org.Hs.eg.db
### ------------


### -------
### biomaRt
### -------


### -------------
### AnnotationHub
### -------------


sessionInfo()

