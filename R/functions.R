#' Create Smith_Waterman (RefClass) object
#'
#' @param X A nucleotide sequence
#' @param Y Another nucleotide sequence
#' @param match Score associated with a match
#' @param mismatch Penalty associated with a mismatch
#' @param d Penalty associated with a gap/indel
#' @param aligned_X The aligned nucleotide sequence obtained according to Smith-Waterman algorithm
#' @param aligned_Y The other aligned nucleotide sequence obtained according to Smith-Waterman algorithm
#' @param alignment_score Best alignment score obtained according to Smith-Waterman algorithm
#' @param F Scoring matrix obtained according to Smith-Waterman algorithm
#' @return A new Smith_Waterman object
#' @author Emanuele Cavalleri\cr Politecnico di Milano\cr Maintainer: Emanuele
#' Cavalleri\cr E-Mail: <emanuele.cavalleri@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm}\cr
Smith_Waterman <- setRefClass("Smith_Waterman",
                              fields=list(X = "character", Y = "character",
                                          match = "numeric", mismatch = "numeric",
                                          d = "numeric",
                                          aligned_X = "character", aligned_Y = "character",
                                          alignment_score = "numeric", F = "matrix"),
                              
                              methods=list(
                                show = function(){
                                  cat("First Sequence: ", X, "\n")
                                  cat("Second Sequence: ", Y, "\n")
                                  cat("Scoring System: ", match, " for match; ",
                                      mismatch, " for mismatch; ", d, " for gap/indel\n")
    
                                  cat("\nDynamic Programming Matrix:\n")
                                  print(F)
    
                                  cat("\nBest local alignment found:\n")
                                  cat(paste(aligned_X, collapse=''), "\n")
                                  cat(paste(aligned_Y, collapse=''), "\n")
                                  cat("\nOptimum alignment score: ", alignment_score)
                                },
                              
                                get_alignment_score = function(){
                                  return(alignment_score)
                                  # N.B.: we could've used Smith_Waterman$accessors(c("alignment_score"))
                                  # but it returns getter+setter, and we don't want a setter
                                }
                              )
                          ) 

#' Apply Smith-Waterman algorithm to strings
#'
#' This function takes as input two strings and returns
#' one of the possible best local alignment in terms of score
#' according to Smith-Waterman greedy algorithm.
#'
#' @usage smith_waterman_string(X, Y, match, mismatch, d)
#' @param X A nucleotide sequence
#' @param Y Another nucleotide sequence
#' @param match Score associated with a match
#' @param mismatch Penalty associated with a mismatch
#' @param d Penalty associated with a gap/indel
#' @return A Smith_Waterman object
#' @author Emanuele Cavalleri\cr Politecnico di Milano\cr Maintainer: Emanuele
#' Cavalleri\cr E-Mail: <emanuele.cavalleri@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm}\cr
#' @seealso \code{\link{smith_waterman_DNAString}}\cr
#' \code{\link{smith_waterman}}\cr
#' @importFrom "methods" "is"
#' @importFrom "methods" "new"
#' @examples
#' smith_waterman_string("TATCT", "TTCGG", 7, -10, -7)
#' 
#' @export
smith_waterman_string <- function(X, Y, match, mismatch, d) {

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

    # Local alignment with dynamic programming: build scoring matrix F
    F <- build_F(F, seq_X, seq_Y, match, mismatch, d)  

    # Traceback
    aligned_sequences <- traceback(F, seq_X, seq_Y, match, mismatch, d)
    
    result <- Smith_Waterman(X=X, Y=Y, match=match, mismatch=mismatch, d=d,
                             aligned_X=aligned_sequences$aligned_X,
                             aligned_Y=aligned_sequences$aligned_Y,
                             alignment_score=max(F), F=F)
    return(result)
}

#' Build scoring matrix F
#' 
#' This function builds and returns the scoring matrix
#' F according to Smith-Waterman greedy algorithm 
#'
#' @param F Pre-initialized scoring matrix
#' @param seq_X A nucleotide sequence preceeded by '0'
#' @param seq_Y Another nucleotide sequence preceeded by '0'
#' @param match Score associated with a match
#' @param mismatch Penalty associated with a mismatch
#' @param d Penalty associated with a gap/indel
#' @return Scoring matrix F obtained according to Smith-Waterman algorithm
#' @author Emanuele Cavalleri\cr Politecnico di Milano\cr Maintainer: Emanuele
#' Cavalleri\cr E-Mail: <emanuele.cavalleri@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm}\cr
#' @seealso \code{\link{smith_waterman_string}}\cr
build_F <- function(F, seq_X, seq_Y, match, mismatch, d) {
  xlen <- length(seq_X)
  ylen <- length(seq_Y)
  
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
  
  return(F)
}

#' Traceback
#' 
#' Traceback implemented according to an efficient modified version
#' of the Needleman-Wunsch algorithm (see Wikipedia reference link)
#'
#' @param F Scoring matrix obtained according to Smith-Waterman algorithm
#' @param seq_X A nucleotide sequence preceeded by '0'
#' @param seq_Y Another nucleotide sequence preceeded by '0'
#' @param match Score associated with a match
#' @param mismatch Penalty associated with a mismatch
#' @param d Penalty associated with a gap/indel
#' @return A list containing the two locally aligned strings
#' @author Emanuele Cavalleri\cr Politecnico di Milano\cr Maintainer: Emanuele
#' Cavalleri\cr E-Mail: <emanuele.cavalleri@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm#Advanced_presentation_of_algorithm}\cr
#' @seealso \code{\link{smith_waterman_string}}\cr
traceback <- function(F, seq_X, seq_Y, match, mismatch, d) {
  # Initialize aligned subsequences for traceback
  aligned_X <- character()
  aligned_Y <- character()
  
  # Traceback, starting at the highest score in the scoring matrix F
  alignment_score <- max(F)
  max_match <- which(F == alignment_score, arr.ind = TRUE, useNames = FALSE)
  i <- max_match[1,1]
  j <- max_match[1,2]
  
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
  
  return(list(aligned_X=aligned_X, aligned_Y=aligned_Y))
}

#' Apply Smith-Waterman algorithm to DNAString objects
#'
#' This function takes as input two DNAString objects
#' and returns one of the possible best local alignment in terms
#' of score according to Smith-Waterman greedy algorithm.
#'
#' @usage smith_waterman_DNAString(X, Y, match, mismatch, d)
#' @param X A DNAString object
#' @param Y Another DNAString object
#' @param match Score associated with a match
#' @param mismatch Penalty associated with a mismatch
#' @param d Penalty associated with a gap/indel
#' @return A Smith_Waterman object
#' @author Emanuele Cavalleri\cr Politecnico di Milano\cr Maintainer: Emanuele
#' Cavalleri\cr E-Mail: <emanuele.cavalleri@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm}\cr
#' @seealso \code{\link{smith_waterman_string}}\cr
#' \code{\link{smith_waterman}}\cr
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
  
    result <- smith_waterman_string(as.character(X), as.character(Y),
                                    match, mismatch, d)
    return(result)
}

#' Apply Smith-Waterman algorithm and show its output
#'
#' This function takes as input two strings/DNAString objects and returns
#' one of the possible best local alignment in terms
#' of score according to Smith-Waterman greedy algorithm.
#'
#' @param X A nucleotide sequence or a DNAString object
#' @param Y Another nucleotide sequence or another DNAString object
#' @param match Score associated with a match
#' @param mismatch Penalty associated with a mismatch
#' @param d Penalty associated with a gap/indel
#' @return A Smith-Waterman object
#' @author Emanuele Cavalleri\cr Politecnico di Milano\cr Maintainer: Emanuele
#' Cavalleri\cr E-Mail: <emanuele.cavalleri@@mail.polimi.it>
#' @references \url{https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm}\cr
#' @seealso \code{\link{smith_waterman_string}}\cr
#' \code{\link{smith_waterman_DNAString}}\cr
#' @import Biostrings
#' @importFrom "methods" "is"
#' @examples
#' library(Biostrings)
#' smith_waterman(DNAString("TATCT"), "TTCGG", 7, -10, -7)
#' 
#' @export
smith_waterman <- function(X, Y, match, mismatch, d) {
    # Check X and Y sequences are strings or DNAString objects
    if(!(is.character(X) || is(X,"DNAString") || is(X,"DNAStringSet")))
        stop("X is not a string or DNAString!")
    if(!(is.character(Y) || is(Y,"DNAString") || is(Y,"DNAStringSet")))
        stop("Y is not a string or DNAString!")
  
    X <- as.character(X)
    Y <- as.character(Y)
    result <- smith_waterman_string(X, Y, match, mismatch, d)
    
    return(result)
}
