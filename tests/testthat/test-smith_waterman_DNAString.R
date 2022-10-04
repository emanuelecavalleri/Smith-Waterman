#context("Apply Smith-Waterman algorithm on DNAString(Set) objects")
test_that("Apply Smith-Waterman algorithm on DNAString objects", {
    X <- "TATCT"
    Y <- "TAGTCT"
    match <- 7
    mismatch <- -10
    d <- -7
    
    
    aligned_X <- c("T", "A", "-", "T", "C", "T")
    aligned_Y <- c("T", "A", "G", "T", "C", "T")
    alignment_score <- 28
    F <- matrix()
    # Check at: https://gtuckerkellogg.github.io/pairwise/demo/
    F <- rbind(c(0,0,0,0,0,0,0),c(0,7,0,0,7,0,7),c(0,0,14,7,0,0,0),
               c(0,7,7,4,14,7,7),c(0,0,0,0,7,21,14),c(0,7,0,0,7,14,28))
    rownames(F) <- c(0, unlist(strsplit(X, '')))
    colnames(F) <- c(0, unlist(strsplit(Y, '')))
    expected <- list(aligned_X, aligned_Y, alignment_score, F)
    expect_equal(smith_waterman_DNAString(DNAString(X), DNAString(Y), match,
                                          mismatch, d), expected)
    expect_equal(smith_waterman_DNAString(DNAStringSet(X)[1], DNAStringSet(Y)[1],
                                          match, mismatch, d), expected)
})

test_that("returns error if not DNAString object", {
  expect_error(smith_waterman_DNAString(DNAString("TATCT"), "TATCT", 7, -10, -7))
})


