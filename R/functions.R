#' Apply Smith-Waterman algorithm
#'
#' This function take as input two strings and returns
#' one of the possible best local alignment in terms of score
#' according to Smith-Waterman greedy algorithm.
#'
#' @usage smith_waterman(X, Y, match, mismatch, d)
#' @param X A nucleotide sequence
#' @param Y Another nucleotide sequence
#' @param match Score associated with a match
#' @param mismatch Penalty associated with a mismatch
#' @param d Penalty associated with a gap/indel
#' @return A list containing the local alignment with scoring matrix and best score
#' @author Emanuele Cavalleri\cr Politecnico di Milano\cr Maintainer: Emanuele
#' Cavalleri\cr E-Mail: <emanuele.cavalleri@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm}\cr
#' @seealso \code{\link{smith_waterman_DNAString}}\cr
#' \code{\link{print_smith_waterman}}\cr
#' @importFrom "methods" "is"
#' @examples
#' smith_waterman("TATCT", "TTCGG", 7, -10, -7)
#' 
#' @export
smith_waterman <- function(X, Y, match, mismatch, d) {

    # Check X and Y sequences are both strings
    if(!is.character(X)) stop("X is not a string!")
    if(!is.character(Y)) stop("Y is not a string!")

    # Check match, mismatch and d are all numeric
    if (!(is.numeric(match) && is.numeric(mismatch) && is.numeric(d)))
        stop("Argument is not numeric!")
    
    # Penalties must be 0 or negative
    if (mismatch > 0) mismatch <- -mismatch
    if (d > 0) d <- -d
    
    # '0' followed by single characters of X (needed to build the scoring matrix)
    seq_X <- c(0, unlist(strsplit(X, '')))
    seq_Y <- c(0, unlist(strsplit(Y, '')))

    # Initialization of the scoring matrix F
    xlen <- length(seq_X)
    ylen <- length(seq_Y)
    F <- matrix(NA, xlen, ylen)
    # F(i,0) = 0
    F[,1] <- integer(xlen)
    # F(0,j) = 0
    F[1,] <- integer(ylen)
    dimnames(F) <- list(seq_X, seq_Y)

    # Local alignment with dynamic programming, alignment equation:
    # F(i,j) = max{ F(i-1, j-1) + submatrix(x[i], y[j]),
    #               F(i-1, j) + d, F(i, j-1) + d, 0 }
    for (i in 2:xlen) {
      for (j in 2:ylen) {
        # F(i,j) = max{ F(i-1, j-1) + submatrix(x[i], y[j]),
        if ( seq_X[i] == seq_Y[j] ) aligned <- F[i-1, j-1] + match # x[i] aligned to y[j]
        else aligned <- F[i-1, j-1] + mismatch # x[i] aligned to y[j]
        # F(i-1, j) + d,
        delete <- F[i-1, j] + d # x[i] aligned to -
        # F(i, j-1) + d,
        insert <- F[i, j-1] + d # y[j] aligned to -
        # 0 }
        F[i, j] <- max(aligned, delete, insert, 0)
      }
    }

    # Traceback, starting at the highest score in the scoring matrix F
    alignment_score <- max(F)
    max_match <- which(F == alignment_score, arr.ind = TRUE, useNames = FALSE)
    i <- max_match[1,1]
    j <- max_match[1,2]
    
    # Initialize aligned subsequences
    aligned_X <- character()
    aligned_Y <- character()
    
    while (i > 1 && j > 1) {
      # Ending at a matrix cell that has a score of 0
      if (F[i, j] == 0) break
      # if x[i] was aligned to y[j] (could've been be a match or a mismatch)
      if ((seq_X[i] == seq_Y[j] && F[i, j] == F[i-1, j-1] + match) ||
          (seq_X[i] != seq_Y[j] && F[i, j] == F[i-1, j-1] + mismatch)) { 
            aligned_X <- c(seq_X[i], aligned_X)
            aligned_Y <- c(seq_Y[j], aligned_Y)
            i <- i - 1
            j <- j - 1
      }
      # if x[i] was aligned to -
      else if (F[i, j] == F[i-1, j] + d) {
          aligned_X <- c(seq_X[i], aligned_X)
          aligned_Y <- c("-", aligned_Y)
          i <- i - 1
      }
      # if y[j] was aligned to -
      else { # if (F[i, j] == F[i, j-1] + d)
          aligned_X <- c("-", aligned_X)
          aligned_Y <- c(seq_Y[j], aligned_Y)
          j <- j - 1
      }
    }
    
    result <- list(aligned_X, aligned_Y, alignment_score, F)
    return(result)
}

#' Apply Smith-Waterman algorithm to DNAString objects
#'
#' This function take as input two DNAString objects
#' and returns one of the possible best local alignment in terms
#' of score according to Smith-Waterman greedy algorithm.
#'
#' @usage smith_waterman_DNAString(X, Y, match, mismatch, d)
#' @param X A DNAString object
#' @param Y Another DNAString object
#' @param match Score associated with a match
#' @param mismatch Penalty associated with a mismatch
#' @param d Penalty associated with a gap/indel
#' @return A list containing the local alignment with scoring matrix and best score
#' @author Emanuele Cavalleri\cr Politecnico di Milano\cr Maintainer: Emanuele
#' Cavalleri\cr E-Mail: <emanuele.cavalleri@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm}\cr
#' @seealso \code{\link{smith_waterman}}\cr
#' \code{\link{print_smith_waterman}}\cr
#' @import Biostrings
#' @importFrom "methods" "is"
#' @examples
#' library(Biostrings)
#' smith_waterman_DNAString(DNAString("TATCT"), DNAString("TTCGG"), 7, -10, -7)
#' 
#' @export
smith_waterman_DNAString <- function(X, Y, match, mismatch, d) {
    
    # Check X and Y sequences are both DNAString objects
    if(!(is(X,"DNAString") || is(X,"DNAStringSet"))) stop("X is not a DNAString!")
    if(!(is(Y,"DNAString") || is(Y,"DNAStringSet"))) stop("Y is not a DNAString!")
  
    return(smith_waterman(as.character(X), as.character(Y), match, mismatch, d))
}

#' Apply Smith-Waterman algorithm and show its output
#'
#' This function takes as input two strings/DNAString objects and returns
#' a visualization of one of the possible best local alignment in terms of 
#' score according to Smith-Waterman greedy algorithm.
#'
#' @param X A nucleotide sequence
#' @param Y Another nucleotide sequence
#' @param match Score associated with a match
#' @param mismatch Penalty associated with a mismatch
#' @param d Penalty associated with a gap/indel
#' @return A visualization of the local alignment with scoring matrix and best score
#' @author Emanuele Cavalleri\cr Politecnico di Milano\cr Maintainer: Emanuele
#' Cavalleri\cr E-Mail: <emanuele.cavalleri@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm}\cr
#' @seealso \code{\link{smith_waterman}}\cr
#' \code{\link{smith_waterman_DNAString}}\cr
#' @import Biostrings
#' @importFrom "methods" "is"
#' @examples
#' library(Biostrings)
#' print_smith_waterman(DNAString("TATCT"), "TTCGG", 7, -10, -7)
#' 
#' @export
print_smith_waterman <- function(X, Y, match, mismatch, d) {
    # Check X and Y sequences are both Strings or DNAString objects
    if(!(is.character(X) || is(X,"DNAString") || is(X,"DNAStringSet")))
        stop("X is not a string or DNAString!")
    if(!(is.character(Y) || is(Y,"DNAString") || is(Y,"DNAStringSet")))
        stop("Y is not a string or DNAString!")
  
    X <- as.character(X)
    Y <- as.character(Y)
    sw <- smith_waterman(X, Y, match, mismatch, d)
    
    #Reconstruct output
    aligned_X <- sw[[1]]
    aligned_Y <- sw[[2]]
    alignment_score <- sw[[3]]
    F <- sw[[4]]

    cat("First Sequence: ", X, "\n")
    cat("Second Sequence: ", Y, "\n")
    cat("Scoring System: ", match, " for match; ", mismatch, " for mismatch; ", d, " for gap/indel\n")

    cat("\nDynamic Programming Matrix:\n")
    print(F)

    cat("\nAlignment:\n")
    cat(paste(aligned_X, collapse=''), "\n")
    cat(paste(aligned_Y, collapse=''), "\n")
    cat("\nOptimum alignment score: ", alignment_score)
}
