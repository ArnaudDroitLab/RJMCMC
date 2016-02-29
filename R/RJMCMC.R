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
#' Astrid Deschênes,
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

#' Nucleosomes obtained by running RJMCMC function
#' (for demo purpose).
#'
#' A \code{list} of \code{class}
#' "rjmcmcNucleosomes" which contains the information about the
#' detected nucleosomes.
#'
#'
#' @name RJMCMC_result
#'
#' @docType data
#'
#' @aliases RJMCMC_result
#'
#' @format A \code{list} of \code{class} "rjmcmcNucleosomes" containing:
#' \itemize{
#' \item \code{call} the matched call.
#' \item \code{K} a \code{vector} of \code{integer}, the estimation of the
#' number of the nucleosomes for each iteration.
#' \item \code{k} a \code{integer}, the final estimation of the number
#' of nucleosomes.
#' \item \code{mu} a \code{vector} of \code{numeric} of length
#' \code{k}, the positions of the nucleosomes.
#' \item \code{sigmaf} a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the forward reads for each nucleosome.
#' \item \code{sigmar} a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the reverse reads for each nucleosome.
#' \item \code{delta} a \code{vector} of \code{numeric} of length
#' \code{k}, the distance between the maxima of the forward and reverse reads
#' position densities for each nucleosome.
#' \item \code{df} a \code{vector} of \code{numeric} of length
#' \code{k}, the degrees of freedom for each nucleosome.
#' \item \code{w} a \code{vector} of positive \code{numerical} of length
#' \code{k}, the weight for each nucleosome. The sum of all \code{w} values
#' must be equal to \code{1}.
#' \item \code{qmu} a \code{matrix} of \code{numerical} with a number of rows
#' of \code{k}, the 2.5\% and 97.5\% quantiles of each \code{mu}.
#' \item \code{qsigmaf} a \code{matrix} of \code{numerical} with a number of
#' rows of \code{k}, the 2.5\% and 97.5\% quantiles of the variance of the
#' forward reads for each nucleosome.
#' \item \code{qsigmar} a \code{matrix} of \code{numerical} with a number of
#' rows of \code{k}, the 2.5\% and 97.5\% quantiles of the variance the
#' reverse reads for each nucleosome.
#' \item \code{qdelta} a \code{matrix} of \code{numerical} with a number of
#' rows of \code{k}, the 2.5\% and 97.5\% quantiles of the distance between
#' the maxima of the forward and reverse reads
#' position densities for each nucleosome.
#' \item \code{qdf} a \code{matrix} of \code{numerical} with a number of
#' rows of \code{k}, the 2.5\% and 97.5\% quantiles of the degrees of freedom
#' for each nucleosome.
#' \item \code{qw} a \code{matrix} of \code{numerical} with a number of rows
#' of \code{k}, the 2.5\% and 97.5\% quantiles of the weight for each
#' nucleosome.
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{rjmcmc}} {for profiling of nucleosome positions}
#' }
#'
#' @usage data(RJMCMC_result)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading dataset
#' data(RJMCMC_result)
#' data(reads_demo)
#'
#' ## Results before post-treatment
#' RJMCMC_result$mu
#'
#' ## Post-treatment function which merged closely positioned nucleosomes
#' postResult <- postTreatment(startPosForwardReads = reads_demo$readsForward,
#' startPosReverseReads = reads_demo$readsReverse, RJMCMC_result, 74, 73500)
#'
#' ## Results after post-treatment
#' postResult
#'
NULL


#' Simulated dataset of reads generated by \code{nucleoSim} package
#' (for demo purpose).
#'
#' A \code{list} of \code{class}
#' "syntheticNucReads" which contains the information about synthetic reads
#' related to nucleosomes. The datset has been created using a total of 300
#' well-positioned nucleosomes, 30 fuzzy nucleosomes with variance of reads
#' following a Normal distribution.
#'
#' @name syntheticNucleosomeReads
#'
#' @docType data
#'
#' @aliases syntheticNucleosomeReads
#'
#' @format A \code{list} containing:
#' \itemize{
#'     \item \code{call} the called that generated the dataset.
#'     \item \code{dataIP} a \code{data.frame} with the chromosome name, the
#'     starting and ending positions and the direction of all forward and
#'     reverse reads for all well-positioned and fuzzy nucleosomes. Paired-end
#'     reads are identified with an unique id.
#'     \item \code{wp} a \code{data.frame} with the positions of all the
#'     well-positioned nucleosomes, as well as the number of paired-reads
#'     associated to each one.
#'     \item \code{fuz} a \code{data.frame} with the positions of all the
#'     fuzzy nucleosomes, as well as the number of paired-reads associated
#'     to each one.
#'     \item \code{paired} a \code{data.frame} with the starting and ending
#'     positions of the reads used to generate the paired-end reads.
#'     Paired-end reads are identified with an unique id.
#' }
#'
NULL
