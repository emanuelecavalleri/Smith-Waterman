## ----style, include = FALSE---------------------------------------------------
library(BiocStyle)
library(knitr)

## ---- message=FALSE-----------------------------------------------------------
library(Biostrings)
library(benchmarkme)
library(ShortRead)

library(smithWaterman)

## -----------------------------------------------------------------------------
X <- "TATCT"
Y <- "TAGTCT"
match <- 7
mismatch <- -10
gap <- -7
# Same with: mismatch <- 10 and/or gap <- 7

sw <- smith_waterman_string(X, Y, match, mismatch, gap)
sw

## -----------------------------------------------------------------------------
sw$show()

## -----------------------------------------------------------------------------
X <- DNAString(X)
Y <- DNAString(Y)

sw <- smith_waterman_DNAString(X, Y, match, mismatch, gap)
sw$show()

## -----------------------------------------------------------------------------
data(phiX174Phage)
X <- subseq(phiX174Phage[1], 1, 500) # First 500 nucleotides
Y <- subseq(phiX174Phage[2], 1, 500)

start_time <- Sys.time() # Start the stopwatch
# Execute it for ten times
for (i in 1:10)
  sw <- smith_waterman_DNAString(X, Y, match, mismatch, gap)
end_time <- Sys.time() # Stop the stopwatch

cat("Best local alignment's score according to Smith-Waterman algorithm:",
    sw$get_alignment_score())

## ----out.width="20%", echo = FALSE--------------------------------------------
url <- "https://www.humboldtmfg.com/product-originals/_lrg@1x/H-2260.png"
knitr::include_graphics(url)

## -----------------------------------------------------------------------------
cat("Mean execution time", round(end_time - start_time, 1) / 10, "second(s) on a(n)",
    get_cpu()$model_name, "having", round(get_ram()/1000000000, 2), "GB of RAM")

## -----------------------------------------------------------------------------
fastqfile <- system.file("extdata", "SP1.fq", package="smithWaterman")
seqset <- readDNAStringSet(fastqfile, format="fastq", with.qualities=TRUE)
X <-  seqset[1]
Y <-  seqset[246]

sw <- smith_waterman_DNAString(X, Y, match, mismatch, gap)
cat("Best local alignment's score according to Smith-Waterman algorithm:",
    sw$alignment_score)

## -----------------------------------------------------------------------------
data(HNF4alpha)
X <-  HNF4alpha[1]
Y <-  HNF4alpha[6]

sw <- smith_waterman_DNAString(X, Y, match, mismatch, gap)
sw$show()

## -----------------------------------------------------------------------------
sw <- smith_waterman(HNF4alpha[1], "AGGCCAAAGTCCT", match, mismatch, gap)
sw$show()

## -----------------------------------------------------------------------------
sessionInfo()

