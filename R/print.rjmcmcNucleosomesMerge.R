#' @title Formated output of predicted nucleosomes
#'
#' @description Generated a formated output of a list marked as
#' an \code{rjmcmcNucleosomesMerge} class
#'
#' @method print rjmcmcNucleosomesMerge
#'
#' @param x the output object from \code{mergeAllRDSFilesFromDirectory}
#' function to be printed
#'
#' @param \ldots arguments passed to or from other methods
#'
#' @examples
#'
#' ## Use a directory present in the RJMCMC package
#' directoryWithRDSFiles <- system.file("extdata", package = "RJMCMC")
#'
#' ## Merge nucleosomes info from RDS files present in directory
#' ## It is assumed that all files present in the directory are nucleosomes
#' ## result for the same chromosome
#' result <- mergeAllRDSFilesFromDirectory(directoryWithRDSFiles)
#'
#' ## Show resulting nucleosomes
#' print(result)
#'
#' ## or simply
#' result
#'
#' @author Astrid Deschenes
#' @export
print.rjmcmcNucleosomesMerge <- function(x, ...) {
    # Print title before printing the content of the regression object
    cat("\nCall:\n")
    print(x$call, ...)
    cat("\nNumber of nucleosomes:\n")
    print(x$k, ...)
    cat("\nNucleosomes positions:\n")
    print(x$mu, ...)
    invisible(x)
}
