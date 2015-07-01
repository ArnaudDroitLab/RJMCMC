#' rjmcmc: Nucleosome Positioning
#'
#' This package does nucleosome positioning using informative
#' Multinomial-Dirichlet prior in a t-mixture with reversible jump
#' estimation of nucleosome positions for genome-wide.
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
#' @name reads_demo
#'
#' @docType data
#'
#' @aliases reads_demo
#'
#' @format A \code{list} containing:
#' \itemize{
#'     \item \code{yf} a \code{vector} of \code{numeric}, the positions of the
#' forward reads.
#'     \item \code{yr} a \code{vector} of \code{numeric}, the positions of the
#' reverse reads.
#' }
#'
#' @source The Encyclopedia of DNA Elements (ENCODE) (DCC accession:
#' ENCFF002CFN)
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
#' RJMCMC(yf = reads_demo$readsForward, yr = reads_demo$readsReverse,
#'      nbrIterations = 100, lambda = 3, kmax = 30,
#'      minInterval = 146, maxInterval = 292, minReads = 5)
#'
NULL

