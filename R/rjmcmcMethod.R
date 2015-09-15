#' @title Nucleosome positioning mapping
#'
#' @description Use of a fully Bayesian hierarchical model for genome-wide
#' profiling of nucleosome positions based on high-throughput short-read
#' data (MNase-Seq data).
#'
#' @param startPosForwardReads a \code{vector} of \code{numeric}, the
#' start position of all the forward reads.
#'
#' @param startPosReverseReads a \code{vector} of \code{numeric}, the
#' start position of all the reverse reads. Beware that the start position of
#' a reverse read is always higher that the end positition.
#'
#' @param nbrIterations a positive \code{integer} or \code{numeric}, the
#' number of iterations. Non-integer values of
#' \code{nbrIterations} will be casted to \code{integer} and truncated towards
#' zero.
#'
#' @param kMax a positive \code{integer} or \code{numeric}, the maximum number
#' of nucleosomes per region. Non-integer values
#' of \code{kMax} will be casted to \code{integer} and truncated towards zero.
#'
#' @param lambda a positive \code{numeric}, the theorical mean
#' of the Poisson distribution.
#'
#' @param minInterval a \code{numeric}, the minimum distance between two
#' nucleosomes.
#'
#' @param maxInterval a \code{numeric}, the maximum distance between two
#' nucleosomes.
#'
#' @param minReads a positive \code{integer} or \code{numeric}, the minimum
#' number of reads in a potential canditate region. Non-integer values
#' of \code{minReads} will be casted to \code{integer} and truncated towards
#' zero.
#'
#' @param verbose a \code{logical} indicating if extra information must be
#' printed or not. Default: \code{FALSE}.
#'
#' @return a \code{list} of \code{class} "rjmcmcNucleosomes" containing :
#' \itemize{
#' \item \code{call} the matched call.
#' \item \code{K} a \code{vector} of \code{integer}, the number of
#' the nucleosomes for each iteration.
#' \item \code{k} a \code{integer}, the number of nucleosomes.
#' \item \code{mu} a \code{vector} of \code{numeric} of length
#' \code{k}, the positions of the nucleosomes.
#' \item \code{sigmaf} a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the forward reads for each nucleosome.
#' \item \code{sigmar} a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the reverse reads for each nucleosome.
#' \item \code{delta} a \code{vector} of \code{numeric} of length
#' \code{k}, the distance between the maxima of the forward and reverse reads
#' position densities for each nucleosome.
#' \item \code{dl} TODO
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
#' \item \code{qdl} TODO
#' \item \code{qw} a \code{matrix} of \code{numerical} with a number of rows
#' of \code{k}, the 2.5\% and 97.5\% quantiles of each \code{w}.
#' }
#'
#' @examples
#'
#' ## Loading dataset
#' data(reads_demo)
#'
#' ## Nucleosome positioning
#' result <- RJMCMC(startPosForwardReads = reads_demo$readsForward,
#'          startPosReverseReads = reads_demo$readsReverse,
#'          nbrIterations = 1000, lambda = 2, kMax = 30,
#'          minInterval = 146, maxInterval = 292, minReads = 5)
#'
#' ## Print the number of nucleosomes
#' result$k
#'
#' ## Print the position of nucleosomes
#' result$mu
#'
#' @importFrom MCMCpack ddirichlet rdirichlet
#' @importFrom stats dmultinom dpois var rmultinom dt quantile
#' @importFrom S4Vectors Rle
#' @import BiocGenerics
#' @author Rawane Samb
#' @export
RJMCMC <- function(startPosForwardReads, startPosReverseReads,
                    nbrIterations, kMax, lambda,
                    minInterval, maxInterval, minReads,
                    verbose = FALSE)
{
    # Get call information
    cl <- match.call()

    # Parameters validation
    validateParameters(startPosForwardReads = startPosForwardReads,
                                startPosReverseReads = startPosReverseReads,
                                nbrIterations = nbrIterations,
                                kMax = kMax,
                                lambda = lambda,
                                minInterval = minInterval,
                                maxInterval = maxInterval,
                                minReads = minReads,
                                verbose = verbose)

    # Casting specific inputs as integer
    minReads        <- as.integer(minReads)
    nbrIterations   <- as.integer(nbrIterations)
    kMax            <- as.integer(kMax)

    ##############################################################
    #### Parameter Initialization                             ####
    ##############################################################

    y               <- c(startPosForwardReads, startPosReverseReads)
    nf              <- length(startPosForwardReads)
    nr              <- length(startPosReverseReads)
    nbrReads        <- nf + nr

    # Order reads an mark reverse reads as -1 in a new vector
    d       <- c(rep(1, nf), rep(-1, nr))
    yOrder  <- order(y)
    y       <- y[yOrder]
    d       <- d[yOrder]
    rm(yOrder)

    ## ASTRID : voir si zeta, detamin, deltamax devraient etre des integer
    zeta            <- 147
    deltamin        <- 142
    deltamax        <- 152

    # Max and min read positions
    minReadPos <- min(y)
    maxReadPos <- max(y)

    paramValues = list(startPSF =  startPosForwardReads,
                        startPSR = startPosReverseReads,
                        kmax = kMax,
                        lambda = lambda,
                        minReads = minReads,
                        y = y, nr = nr, nf = nf, nbrReads = nbrReads,
                        zeta = zeta, deltamin = deltamin, deltamax = deltamax,
                        minReadPos = minReadPos, maxReadPos = maxReadPos,
                        nbrIterations = nbrIterations)

    # Vector of the number of nucleosomes (integer values)
    k               <- Rle(rep(0L, nbrIterations))

    # Vector of the position of the nucleosomes
    mu              <- matrix(0, nrow = nbrIterations, ncol = kMax)

    sigmaf          <- matrix(0, nrow = nbrIterations, ncol = kMax)

    sigmar          <- matrix(0, nrow = nbrIterations, ncol = kMax)

    delta           <- matrix(0, nrow = nbrIterations, ncol = kMax)

    w               <- matrix(0, nrow = nbrIterations, ncol = kMax)

    dl              <- matrix(0L, nrow = nbrIterations, ncol = kMax)

    k[1]            <- 1L

    mu[1, 1]        <- runif(1, minReadPos, maxReadPos)

    sigmaf[1, 1]    <- 1
    sigmar[1, 1]    <- 1
    delta[1, 1]     <- runif(1, 0, 2*(mu[1, 1] - minReadPos))
    w[1, 1]         <- 1
    dl[1, 1]        <- 3

    nbrIterations   <- ifelse(nbrReads <= 10, 1000, nbrIterations)

    kValue          <- as.integer(k[1])


    ## Initialization of current values used in the loop
    muValue       <- mu[1,]
    sigmafValue   <- sigmaf[1,]
    sigmarValue   <- sigmar[1,]
    deltaValue    <- delta[1,]
    wValue        <- w[1,]
    dlValue       <- dl[1,]
    aValue        <- rep(0, kMax + 1L)
    aValue[1]     <- minReadPos
    aValue[as.integer(k[1]) + 1L]  <- maxReadPos
    dimValue       <- rep(0, kMax)
    dimValue[1]    <- nbrReads
    ktildeValue    <- as.integer(k[1])

    for (i in 2:nbrIterations) {

        ## List of current values that is going to be passed to sub-functions
        varTilde <- list()

        if (kValue == 1L) {
            ## CASE : Number of nucleosomes equal to 1
            u <- runif(1)

            if (u <= 0.5) {

                # Birth move in case k=1
                varTilde <- birthMoveK1(paramValues, kValue, muValue,
                                        sigmafValue, sigmarValue, deltaValue,
                                        wValue, dlValue, aValue, dimValue)
            } ### end of B move in case k=1
            else {
                ###Metropolis-Hastings move
                varTilde <- mhMoveK1(paramValues, kValue, muValue,
                                        sigmafValue, sigmarValue, deltaValue,
                                        wValue, dlValue, aValue, dimValue)

            }    ### end of M-H move in case k=1
        }  ### end of test of k=1
        else {
            ## CASE : Number of nucleosomes larger than 1
            u<-runif(1)

            if (u <= Dk(kValue, lambda, kMax)) {
                ### Death move
                varTilde <- deathMove(paramValues, kValue, muValue,
                                        sigmafValue, sigmarValue, deltaValue,
                                        wValue, dlValue, aValue, dimValue )

            } ### end of Death move
            else {

                if (u <= (Dk(kValue, lambda, kMax) +
                                Bk(kValue, lambda, kMax))) {
                    ### Birth move
                    varTilde <- birthMove(paramValues, kValue, muValue,
                                            sigmafValue, sigmarValue,
                                            deltaValue, wValue, dlValue,
                                            aValue, dimValue)

                } ### end of Birth move
                else {
                    ### Metropolis-Hastings move
                    varTilde <- mhMove(paramValues, kValue, muValue,
                                            sigmafValue, sigmarValue,
                                            deltaValue, wValue, dlValue,
                                            aValue, dimValue)

                } ### end of Metropolis-Hastings move
            } #end of else of Birth move and Metropolis-Hastings move
        } ###end of moves in case larger than 1

        v <- runif(1)      #Acceptation/rejet

        if (varTilde$rho >= v && varTilde$k <= kMax) {
            # Acceptation, so values are updated
            kValue          <- varTilde$k
            maxValue        <- as.integer(kValue)
            zeroVector      <- rep(0, kMax - maxValue)
            muValue         <- c(varTilde$mu[1:maxValue], zeroVector)
            sigmafValue     <- c(varTilde$sigmaf[1:maxValue], zeroVector)
            sigmarValue     <- c(varTilde$sigmar[1:maxValue], zeroVector)
            deltaValue      <- c(varTilde$delta[1:maxValue], zeroVector)
            dlValue         <- c(varTilde$dl[1:maxValue], zeroVector)
            wValue          <- c(varTilde$w[1:maxValue], zeroVector)
            dimValue        <- c(varTilde$dim[1:maxValue], zeroVector)
            aValue          <- c(varTilde$a[1:(maxValue + 1)], zeroVector)
        }

        # Si on veut faire le merge on peut le faire ici
        # et le conserver sans changer les varValue
        k[i]                    <- kValue
        maxValue                <- as.integer(k[i])
        mu[i, 1:maxValue]       <- muValue[1:maxValue]
        sigmaf[i, 1:maxValue]   <- sigmafValue[1:maxValue]
        sigmar[i, 1:maxValue]   <- sigmarValue[1:maxValue]
        delta[i, 1:maxValue]    <- deltaValue[1:maxValue]
        dl[i, 1:maxValue]       <- dlValue[1:maxValue]
        w[i, 1:maxValue]        <- wValue[1:maxValue]

        ## Set the new value of kValue for the next iteration

    } ###end of boucle RJMCMC

    liste <- rep(list(NULL), nbrIterations)

    for (i in 1:nbrIterations)
    {
        kVal <- as.integer(k[i])
        new.list <- list(
                        k      = kVal,
                        mu     = mu[i, 1:kVal],
                        sigmaf = sigmaf[i, 1:kVal],
                        sigmar = sigmar[i, 1:kVal],
                        delta  = delta[i, 1:kVal],
                        dl     = dl[i, 1:kVal],
                        w      = w[i, 1:kVal]
        )

        if (verbose) {
            print("TODO")
        }

        listeUpdate <- mergeNucleosomes(startPosForwardReads,
                                        startPosReverseReads, y, new.list,
                                        minInterval,
                                        maxInterval, minReads)
        kVal                <- listeUpdate$k
        k[i]                <- kVal
        mu[i, ]             <- c(listeUpdate$mu, rep(0, kMax - kVal))
        sigmaf[i, ]   <- c(listeUpdate$sigmaf, rep(0, kMax - kVal))
        sigmar[i, ]   <- c(listeUpdate$sigmar, rep(0, kMax - kVal))
        delta[i, ]    <- c(listeUpdate$delta, rep(0, kMax - kVal))
        w[i, ]        <- c(listeUpdate$w, rep(0, kMax - kVal))
        dl[i, ]       <- c(listeUpdate$dl, rep(0, kMax - kVal))
    }

    ## Astrid : Si la fonction mode() retourne NA, le cas n'est pas gere
    ## Getting the number of nucleosomes with the highest frequency
    km          <- elementWithHighestMode(as.integer(k))
    kPositions  <- which(as.integer(k) == km)

    mu_hat     <- colMeans(mu[kPositions, 1:km, drop=FALSE])
    sigmaf_hat <- colMeans(sigmaf[kPositions, 1:km, drop=FALSE])
    sigmar_hat <- colMeans(sigmar[kPositions, 1:km, drop=FALSE])
    w_hat      <- colMeans(w[kPositions, 1:km, drop=FALSE])
    delta_hat  <- colMeans(delta[kPositions, 1:km, drop=FALSE])
    dl_hat     <- round(colMeans(dl[kPositions, 1:km, drop=FALSE]))

    # Getting 2.5% and 97.5% quantiles for each important data type
    qmu     <- t(apply(mu[, 1:km, drop=FALSE], MARGIN=2, FUN=quantile,
                        probs=c(0.025, 0.975)))
    qsigmaf <- t(apply(sigmaf[, 1:km, drop=FALSE], MARGIN=2, FUN=quantile,
                        probs=c(0.025, 0.975)))
    qsigmar <- t(apply(sigmar[, 1:km, drop=FALSE], MARGIN=2, FUN=quantile,
                        probs=c(0.025, 0.975)))
    qdelta  <- t(apply(delta[, 1:km, drop=FALSE], MARGIN=2, FUN=quantile,
                        probs=c(0.025, 0.975)))
    qdl     <- t(apply(dl[, 1:km, drop=FALSE], MARGIN=2, FUN=quantile,
                        probs=c(0.025, 0.975)))
    qw      <- t(apply(w[, 1:km, drop=FALSE], MARGIN=2, FUN=quantile,
                        probs=c(0.025, 0.975)))

    # Create the final list
    result <- list(
        call    = cl,
        K       = as.integer(k),
        k       = km,
        mu      = mu_hat,
        sigmaf  = sigmaf_hat,
        sigmar  = sigmar_hat,
        delta   = delta_hat,
        dl      = dl_hat,
        w       = w_hat,
        qmu     = qmu,
        qsigmaf = qsigmaf,
        qsigmar = qsigmar,
        qdelta  = qdelta,
        qdl     = qdl,
        qw      = qw
    )

    class(result)<-"rjmcmcNucleosomes"

    return(result)

}
