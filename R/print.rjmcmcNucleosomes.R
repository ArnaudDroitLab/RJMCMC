#' @title Formated output of predicted nucleosomes
#'
#' @description Generated a formated output of a list marked as
#' an \code{rjmcmcNucleosomes} class
#'
#' @method print rjmcmcNucleosomes
#'
#' @param x the output object from \code{rjmcmc} function to be printed
#'
#' @param \ldots arguments passed to or from other methods
#'
#' @examples
#'
#' ## Loading dataset
#' data(RJMCMC_result)
#'
#' print(RJMCMC_result)
#'
#' @author Astrid Deschenes
#' @export
print.rjmcmcNucleosomes <- function(x, ...) {
    # Print title before printing the content
    cat("RJMCMC - Predicted nucleosomes\n")
    cat("\nCall:\n")
    print(x$call, ...)
    cat("\nNumber of nucleosomes:\n")
    print(x$k, ...)
    cat("\nNucleosomes positions:\n")
    print(x$mu, ...)
    invisible(x)
}
