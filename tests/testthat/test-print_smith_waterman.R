#context("Print Smith-Waterman algorithm on DNAString(Set) or canonical strings")
library(Biostrings)
test_that("Print Smith-Waterman algorithm on DNAString(Set) or canonical strings", {
  expect_equal(print_smith_waterman("TATCT", "TAGTCT", 7, -10, -7), print_smith_waterman(DNAString("TATCT"), DNAString("TAGTCT"), 7, -10, -7))
  expect_equal(print_smith_waterman(DNAStringSet("TATCT")[1], DNAStringSet("TAGTCT")[1], 7, -10, -7), print_smith_waterman(DNAString("TATCT"), "TAGTCT", 7, -10, -7))
})

test_that("returns error if not string", {
  expect_error(print_smith_waterman("TATCT", 98, 7, -10, -7))
})