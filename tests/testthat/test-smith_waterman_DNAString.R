#context("Apply Smith-Waterman algorithm on DNAString(Set)")
library(Biostrings)
test_that("Apply Smith-Waterman algorithm on strings", {
  ax <- c("T", "A", "-", "T", "C", "T")
  ay <- c("T", "A", "G", "T", "C", "T")
  F <- matrix()
  # Check at: https://gtuckerkellogg.github.io/pairwise/demo/
  F <- rbind(c(0,0,0,0,0,0,0),c(0,7,0,0,7,0,7),c(0,0,14,7,0,0,0),c(0,7,7,4,14,7,7),c(0,0,0,0,7,21,14),c(0,7,0,0,7,14,28))
  rownames(F) <- c(0, unlist(strsplit("TATCT", '')))
  colnames(F) <- c(0, unlist(strsplit("TAGTCT", '')))
  expected <- list(ax, ay, 28, F)
  expect_equal(smith_waterman_DNAString(DNAString("TATCT"), DNAString("TAGTCT"), 7, -10, -7), expected)
  expect_equal(smith_waterman_DNAString(DNAStringSet("TATCT")[1], DNAStringSet("TAGTCT")[1], 7, -10, -7), expected)
})

test_that("returns error if not string", {
  expect_error(smith_waterman_DNAString(DNAString("TATCT"), "TATCT", 7, -10, -7))
})


