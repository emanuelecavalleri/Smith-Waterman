#' Apply Smith-Waterman algorithm
#'
#' This function take as input two strings and returns
#' one of the possible best local alignment in terms of score
#' according to Smith-Waterman greedy algorithm.
#'
#' @usage smith_waterman(X, Y, match, mismatch, d)
#' @param X A nucleotidic sequence
#' @param Y Another nucleotidic sequence
#' @param match Score of a match
#' @param mismatch Penalty of a mismatch
#' @param d Penalty of a gap/indel
#' @return A list containing a possible local alignment with cost matrix and best score
#' @author Emanuele Cavalleri\cr Politecnico di Milano\cr Maintainer: Emanuele
#' Cavalleri\cr E-Mail: <emanuele.cavalleri@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm}\cr
#' @seealso \code{\link{smith_waterman_DNAString}}\cr
#' \code{\link{print_smith_waterman}}\cr
#' @import Biostrings
#' @import ShortRead
#' @import benchmarkme
#' @importFrom "methods" "is"
#' @examples
#' smith_waterman("TATCT", "TTCGG", 7, -10, -7)
#' 
#' @export
smith_waterman <- function(X, Y, match, mismatch, d) {

    # Check X and Y sequences are both strings
    if(!is.character(X)) stop("X is not a String!")
    if(!is.character(Y)) stop("Y is not a String!")

    # Check match, mismatch and d are all numeric
    if (!is.numeric(match) || !is.numeric(mismatch) || !is.numeric(d)) {
        stop("Argument is not numeric!")
    }
    
    # Penalties must be 0 or negative
    if (mismatch > 0) mismatch = -mismatch
    if (d > 0) d = -d
    
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
      }
      # if x[i] was aligned to -
      else if (F[i, j] == F[i-1, j] + d) {
        ax <- c(seq.x[i], ax)
        ay <- c("-", ay)
        i <- i - 1
      }
      # if y[j] was aligned to -
      else { # if (j > 1 && F[i, j] == F[i, j-1] + d)
        ax <- c("-", ax)
        ay <- c(seq.y[j], ay)
        j <- j - 1
      }
    }
    
    result = list(ax, ay, alignment_score, F)
    return(result)
}

#' Apply Smith-Waterman algorithm with DNAStrings
#'
#' This function take as input two DNAString(Set) (Biostrings) objects
#' and returns one of the possible best local alignment in terms
#' of score according to Smith-Waterman greedy algorithm.
#'
#' @usage smith_waterman_DNAString(X, Y, match, mismatch, d)
#' @param X A DNAString object
#' @param Y Another DNAString object
#' @param match Score of a match
#' @param mismatch Penalty of a mismatch
#' @param d Penalty of a gap/indel
#' @return A list containing a possible local alignment with cost matrix and best score
#' @author Emanuele Cavalleri\cr Politecnico di Milano\cr Maintainer: Emanuele
#' Cavalleri\cr E-Mail: <emanuele.cavalleri@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm}\cr
#' @seealso \code{\link{smith_waterman}}\cr
#' \code{\link{print_smith_waterman}}\cr
#' @import Biostrings
#' @import ShortRead
#' @import benchmarkme
#' @importFrom "methods" "is"
#' @examples
#' library(Biostrings)
#' smith_waterman_DNAString(DNAString("TATCT"), DNAString("TTCGG"), 7, -10, -7)
#' 
#' @export
smith_waterman_DNAString <- function(X, Y, match, mismatch, d) {
  
  
  # Check X and Y sequences are both DNAString object
  if(!is(X,"DNAString") && !is(X,"DNAStringSet")) stop("X is not a DNAString!")
  if(!is(Y,"DNAString") && !is(Y,"DNAStringSet")) stop("Y is not a DNAString!")
  
  return(smith_waterman(as.character(X), as.character(Y), match, mismatch, d))
}

#' Show Smith-Waterman algorithm output
#'
#' This function take as input two strings and returns a visualization
#' of one of the possible best local alignment in terms of score
#' according to Smith-Waterman greedy algorithm.
#'
#' @param X A nucleotidic sequence
#' @param Y Another nucleotidic sequence
#' @param match Score of a match
#' @param mismatch Penalty of a mismatch
#' @param d Penalty of a gap/indel
#' @return A visualization of the two aligned sequences with the best score
#' @author Emanuele Cavalleri\cr Politecnico di Milano\cr Maintainer: Emanuele
#' Cavalleri\cr E-Mail: <emanuele.cavalleri@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm}\cr
#' @seealso \code{\link{smith_waterman}}\cr
#' \code{\link{smith_waterman_DNAString}}\cr
#' @import Biostrings
#' @import ShortRead
#' @import benchmarkme
#' @importFrom "methods" "is"
#' @examples
#' library(Biostrings)
#' print_smith_waterman(DNAString("TATCT"), "TTCGG", 7, -10, -7)
#' 
#' @export
print_smith_waterman <- function(X, Y, match, mismatch, d) {
  
    # Check X and Y sequences are both Strings or DNAString objects
    if(!is.character(X) && !is(X,"DNAString") && !is(X,"DNAStringSet")) stop("X is not a String or DNAString!")
    if(!is.character(Y) && !is(Y,"DNAString") && !is(Y,"DNAStringSet")) stop("Y is not a String or DNAString!")
  
    X <- as.character(X)
    Y <- as.character(Y)
    sw = smith_waterman(X, Y, match, mismatch, d)
    ax = sw[[1]]
    ay = sw[[2]]
    alignment_score = sw[[3]]
    F = sw[[4]]

    cat("First Sequence: ", X, "\n")
    cat("Second Sequence: ", Y, "\n")
    cat("Scoring System: ", match, " for match; ", mismatch, " for mismatch; ", d, " for gap/indel", "\n\n")

    cat("Dynamic Programming Matrix:\n")
    print(F)

    cat("\nAlignment:\n")
    cat(paste(ax, collapse=''), "\n")
    cat(paste(ay, collapse=''),"\n\n")
    cat("Optimum alignment score: ", alignment_score)
}