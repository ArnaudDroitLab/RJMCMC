#' RJMCMC: Bayesian hierarchical model for genome-wide profiling of nucleosome positions based on high-throughput short-read data (MNase-Seq data)
#'
#' This package does nucleosome positioning using informative
#' Multinomial-Dirichlet prior in a t-mixture with reversible jump
#' estimation of nucleosome positions for genome-wide profiling.
#'
#' @docType package
#'
#' @name RJMCMC-package
#'
#' @aliases RJMCMC-package RJMCMC
#'
#' @author  Rawane Samb,
#' Khader Khadraoui,
#' Pascal Belleau,
#' Astrid Deschenes,
#' Lajmi Lakhal and
#' Arnaud Droit
#'
#' Maintainer:
#' Astrid Deschenes <astrid-louise.deschenes@@crchudequebec.ulaval.ca>
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{rjmcmc}} { for profiling of nucleosome positions}
#' }
#'
#' @keywords package
NULL

#' Forward reads and reverse reads in \code{numeric} format (for demo purpose).
#'
#' A group of froward and reverse reads that can be used to test the
#' \code{rjmcmc} function.
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
#'     \item \code{\link{rjmcmc}} {for profiling of nucleosome positions}
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
#' rjmcmc(startPosForwardReads = reads_demo$readsForward,
#'          startPosReverseReads = reads_demo$readsReverse,
#'          nbrIterations = 100, lambda = 3, kMax = 30,
#'          minInterval = 146, maxInterval = 292, minReads = 5)
#'
NULL

#' Forward reads and reverse reads in \code{integer} format
#' (for demo purpose).
#'
#' A group of froward and reverse reads that can be used to test the
#' \code{rjmcmc} function.
#'
#' @name reads_demo_02
#'
#' @docType data
#'
#' @aliases reads_demo_02
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
#'     \item \code{\link{rjmcmc}} {for profiling of nucleosome positions}
#' }
#'
#' @usage data(reads_demo_02)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset
#' data(reads_demo_02)
#'
#' ## Nucleosome positioning
#' rjmcmc(startPosForwardReads = reads_demo_02$readsForward,
#'          startPosReverseReads = reads_demo_02$readsReverse,
#'          nbrIterations = 150, lambda = 3, kMax = 30,
#'          minInterval = 144, maxInterval = 290, minReads = 6)
#'
NULL
