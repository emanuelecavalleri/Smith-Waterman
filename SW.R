# _Short description_: Implement a pairwise local sequence alignment method,
# to find the optimal global sequence alignment between two nucleotide sequences,
# using the Smith-Waterman algorithm.

# _Expected outcome_:
# Develop a package which consist of:
# 1. Appropriate data structures;
# 2. Based on the design in the point (1), implement the Smith-Waterman algorithm
    # for optimal local alignment.
# 3. The expected outcome of the project is a function that takes as parameters:
    # (a1) the cost of a match, (a2) the cost of a mismatch, (a3) the cost of a
    # gap/indel, (b) the first nucleotide sequence, (c) the second nucleotide sequence.
    # As result, returns one optimal alignment between the two input sequences (b)
    # and (c).
# 4. Documentation of the design and development.
# 5. The implementation should work smoothly for sequences of length equals to ~500
    # nucelotides.
# 6. _Optionally_: Unit test implementation

# _Hints/Notes_:
# For the Smith-Waterman algorithm the student may want to define two matrices;
# one to store the scores of the alignment and a second one to store the directions
# of the alignment, which will be used in the trace back phase.

# _Literature/Resources_:
# Smith-Waterman algorithm:
# Smith, Temple F. & Waterman, Michael S. (1981). "Identification of Common
# Molecular Subsequences". Journal of Molecular Biology. 147 (1): 195–197.

rm(list=ls())

# Get sequences
# X <- readline(prompt="Enter the first sequence: ")
# Y <- readline(prompt="Enter the second sequence: ")
X = "TATCT"
Y = "TTCGG"
match <- 7
mismatch <- -10
d <- -7 # gap/indel score

# Add a zero before the sequence and split it into single characters
seq.x <- c(0, unlist(strsplit(X, '')))
seq.y <- c(0, unlist(strsplit(Y, '')))

# Initialization of the score matrix (F)
xlen <- length(seq.x)
ylen <- length(seq.y)
F <- matrix(NA, xlen, ylen)
# F(i,0) = 0
F[,1] <- integer(xlen)
# F(0,j) = 0
F[1,] <- integer(ylen)
dimnames(F) <- list(seq.x, seq.y)

# Local alignment with dynamic programming, alignment equation:
# F(i,j) = max{ F(i-1, j-1) + submatrix(x_i, y_i),
#               F(i-1, j) + d, F(i, j-1) + d, 0 }
for (i in 2:xlen) {
  for (j in 2:ylen) {
    # F(i,j) = max{ F(i-1, j-1) + submatrix(x[i], y[i]),
    if ( seq.x[i] == seq.y[j] ) aligned = F[i-1, j-1] + match # x[i] aligned to y[j]
    else aligned = F[i-1, j-1] + mismatch # x[i] aligned to y[j]
    # F(i-1, j) + d,
    delete <- F[i-1, j] + d # x[i] aligned to -
    # F(i, j-1) + d,
    insert <- F[i, j-1] + d # y[j] aligned to -
    # 0 }
    F[i, j] = max(aligned, delete, insert, 0)
  }
}

# Traceback, starting at the highest score in the scoring matrix F
alignment_score <- max(F)
max_match <- which(F == alignment_score, arr.ind = TRUE, useNames = FALSE)
i <- max_match[1,1]
j <- max_match[1,2]
ax <- character()
ay <- character()
while (i > 1 && j > 1) {
  # Ending at a matrix cell that has a score of 0
  if (F[i, j] == 0) break
  # if x[i] was aligned to y[j]
  if ((seq.x[i] == seq.y[j] && F[i, j] == F[i-1, j-1] + match) ||
      (seq.x[i] != seq.y[j] && F[i, j] == F[i-1, j-1] + mismatch)) { 
    ax <- c(seq.x[i], ax)
    ay <- c(seq.y[j], ay)
    i <- i -1
    j <- j - 1
    print(F[i,j])
  }
  # if x[i] was aligned to -
  else if (F[i, j] == F[i-1, j] + d) {
    ax <- c(seq.x[i], ax)
    ay <- c("-", ay)
    i <- i - 1
    print(F[i,j])
  }
  # if y[j] was aligned to -
  else { # if (j > 1 && F[i, j] == F[i, j-1] + d)
    ax <- c("-", ax)
    ay <- c(seq.y[j], ay)
    j <- j - 1
    print(F[i,j])
  }
}

cat("First Sequence: ", X, "\n")
cat("Second Sequence: ", Y, "\n")
cat("Scoring System: ", match, " for match; ", mismatch, " for mismatch; ", d, " for gap", "\n\n")

cat("Dynamic Programming Matrix:\n")
print(F)

cat("\nAlignment:\n")
cat(paste(ax, collapse=''), "\n")
cat(paste(ay, collapse=''),"\n\n")
cat("Optimum alignment score: ", alignment_score,"\n")

# https://gtuckerkellogg.github.io/pairwise/demo/