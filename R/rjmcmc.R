#' rjmcmc: Nucleosome Positioning
#'
#' This package does nucleosome positioning using informative
#' Multinomial-Dirichlet prior in a t-mixture with reversible jump
#' estimation of nucleosome positions for genome-wide profiling.
#'
#' @docType package
#'
#' @name rjmcmc-package
#'
#' @aliases rjmcmc-package rjmcmc
#'
#' @author  Rawane Samb,
#' Astrid Louise Deschenes,
#' Khader Khadraoui,
#' Lajmi Lakhal,
#' Pascal Belleau and
#' Arnaud Droit
#'
#' Maintainer:
#' Astrid Louise Deschenes <astrid-louise.deschenes@@crchudequebec.ulaval.ca>
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{RJMCMC}} {for doing nucleosomes positioning.}
#' }
#'
#' @keywords package
NULL

#' Forward reads and reverse reads (for demo purpose).
#'
#' A group of froward and reverse reads that can be used to test the
#' \code{RJMCMC} function.
#'
#' @name reads_demo
#'
#' @docType data
#'
#' @aliases reads_demo
#'
#' @format A \code{list} containing:
#' \itemize{
#'     \item \code{readsFoward} a \code{vector} of non-negative \code{numeric},
#'     the start positions of the forward reads.
#'     \item \code{readsReverse} a \code{vector} of non-negative
#'     \code{numeric}, the start
#'     positions of the reverse reads. Beware that the start position of
#'     a reverse read is always higher that the end positition.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{RJMCMC}} {for doing nucleosomes positioning.}
#' }
#'
#' @usage data(reads_demo)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset
#' data(reads_demo)
#'
#' ## Nucleosome positioning
#' RJMCMC(startPosForwardReads = reads_demo$readsForward,
#'          startPosReverseReads = reads_demo$readsReverse,
#'          nbrIterations = 100, lambda = 3, kmax = 30,
#'          minInterval = 146, maxInterval = 292, minReads = 5)
#'
NULL

