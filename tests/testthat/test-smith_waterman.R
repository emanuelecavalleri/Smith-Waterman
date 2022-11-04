#context("Apply Smith-Waterman algorithm on DNAString(Set) objects and/or strings")
test_that("Apply Smith-Waterman algorithm applied on DNAString(Set) objects and/or strings", {
    X <- "TATCT"
    Y <- "TAGTCT"
    match <- 7
    mismatch <- -10
    d <- -7
    
    expect_equal(smith_waterman(X, Y, match, mismatch, d),
                 smith_waterman(DNAString(X), DNAString(Y), match, mismatch, d))
    expect_equal(smith_waterman(DNAStringSet(X)[1], DNAStringSet(Y)[1], match, mismatch, d),
                 smith_waterman(DNAString(X), Y, match, mismatch, d))
})

test_that("returns error if not string/DNAString object", {
    expect_error(smith_waterman("TATCT", 98, 7, -10, -7))
    expect_error(smith_waterman(67, DNAString("TAGTCT"), 7, -10, -7))
})