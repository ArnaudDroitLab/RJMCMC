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
#' @param runSplitFunction a \code{logical} indicating if the split function
#' is run or skipt. Default: \code{TRUE}.
#'
#' @param runMergeFunction a \code{logical} indicating if the merge function
#' is run or skipt. Default: \code{TRUE}.
#'
#' @param adaptIterationsToReads a \code{logical} indicating if the number
#' of iterations must be modified in function of the number of reads.
#' Default: \code{TRUE}.
#'
#' @return a \code{list} of \code{class} "rjmcmcNucleosomes" containing :
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
#' @examples
#'
#' ## Loading dataset
#' data(reads_demo)
#'
#' ## Nucleosome positioning, running both merge and split functions
#' result <- rjmcmc(startPosForwardReads = reads_demo$readsForward,
#'          startPosReverseReads = reads_demo$readsReverse,
#'          nbrIterations = 1000, lambda = 2, kMax = 30,
#'          minInterval = 146, maxInterval = 292, minReads = 5)
#'
#' ## Print the final estimation of the number of nucleosomes
#' result$k
#'
#' ## Print the position of nucleosomes
#' result$mu
#'
#' @importFrom MCMCpack ddirichlet rdirichlet
#' @importFrom stats dmultinom dpois var rmultinom dt quantile
#' @importFrom IRanges IRanges
#' @import BiocGenerics
#' @author Rawane Samb, Pascal Belleau, Astrid Deschenes
#' @export
rjmcmc <- function(startPosForwardReads, startPosReverseReads,
                    nbrIterations, kMax, lambda,
                    minInterval, maxInterval, minReads,
                    runSplitFunction = TRUE,
                    runMergeFunction = TRUE,
                    adaptIterationsToReads = TRUE) {

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
                            runSplitFunction = runSplitFunction,
                            runMergeFunction = runMergeFunction,
                            adaptIterationsToReads = adaptIterationsToReads)

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

    ## Fixed parameters
    zeta            <- 147
    deltamin        <- 142
    deltamax        <- 152

    # Max and min read positions
    minReadPos <- min(y)
    maxReadPos <- max(y)

    # Adapt the number of iterations
    if (adaptIterationsToReads) {
        nbrIterations <- ifelse(nbrReads <= 12, 1000, nbrIterations)
    }

    # List of fixed parameters
    paramValues <- list(startPSF = startPosForwardReads,
                        startPSR = startPosReverseReads,
                        kmax = kMax,
                        lambda = lambda,
                        minReads = minReads,
                        y = y, nr = nr, nf = nf, nbrReads = nbrReads,
                        zeta = zeta, deltamin = deltamin, deltamax = deltamax,
                        minReadPos = minReadPos, maxReadPos = maxReadPos)

    # Vector of the number of nucleosomes (integer values)
    k               <- rep(0L, nbrIterations)

    # Vector of the position of the nucleosomes
    mu              <- matrix(0, nrow = nbrIterations, ncol = kMax)

    sigmaf          <- matrix(0, nrow = nbrIterations, ncol = kMax)

    sigmar          <- matrix(0, nrow = nbrIterations, ncol = kMax)

    delta           <- matrix(0, nrow = nbrIterations, ncol = kMax)

    w               <- matrix(0, nrow = nbrIterations, ncol = kMax)
    # Vector of the degrees of freedom of the nucleosomes
    df              <- matrix(0L, nrow = nbrIterations, ncol = kMax)

    k[1]            <- 1L

    mu[1, 1]        <- runif(1, minReadPos, maxReadPos)

    sigmaf[1, 1]    <- 1
    sigmar[1, 1]    <- 1
    delta[1, 1]     <- runif(1, 0, 2*(mu[1, 1] - minReadPos))
    w[1, 1]         <- 1
    df[1, 1]        <- 3

    kValue          <- as.integer(k[1])

    ## Initialization of current values used in the loop
    muValue                         <- mu[1,]
    sigmafValue                     <- sigmaf[1,]
    sigmarValue                     <- sigmar[1,]
    deltaValue                      <- delta[1,]
    wValue                          <- w[1,]
    dfValue                         <- df[1,]
    aValue                          <- rep(0, kMax + 1L)
    aValue[1]                       <- minReadPos
    aValue[as.integer(k[1]) + 1L]   <- maxReadPos
    dimValue                        <- rep(0, kMax)
    dimValue[1]                     <- nbrReads
    ktildeValue                     <- as.integer(k[1])

    for (i in 2:nbrIterations) {
        ## List of current values that is going to be passed to sub-functions
        varTilde <- list()

        if (kValue == 1L) {
            ## CASE : Number of nucleosomes equal to 1
            u <- runif(1)

            if (u <= 0.5) {
                ## Birth move in case k=1
                varTilde <- birthMoveK1(paramValues, kValue, muValue,
                                        sigmafValue, sigmarValue, deltaValue,
                                        wValue, dfValue, aValue, dimValue)
            } ## end of Birth move in case k=1
            else {
                ## Metropolis-Hastings move
                varTilde <- mhMoveK1(paramValues, kValue, muValue,
                                        sigmafValue, sigmarValue, deltaValue,
                                        wValue, dfValue, aValue, dimValue)

            } ## end Metropolis-Hastings move
        }  ## end CASE : Number of nucleosomes equal to 1
        else {
            ## CASE : Number of nucleosomes larger than 1
            u<-runif(1)

            if (u <= Dk(kValue, lambda, kMax)) {
                ## Death move
                varTilde <- deathMove(paramValues, kValue, muValue,
                                        sigmafValue, sigmarValue, deltaValue,
                                        wValue, dfValue, aValue, dimValue)

            } ## end of Death move
            else {

                if (u <= (Dk(kValue, lambda, kMax) +
                                Bk(kValue, lambda, kMax))) {
                    ### Birth move
                    varTilde <- birthMove(paramValues, kValue, muValue,
                                            sigmafValue, sigmarValue,
                                            deltaValue, wValue, dfValue,
                                            aValue, dimValue)

                } ## end of Birth move
                else {
                    ## Metropolis-Hastings move
                    varTilde <- mhMove(paramValues, kValue, muValue,
                                            sigmafValue, sigmarValue,
                                            deltaValue, wValue, dfValue,
                                            aValue, dimValue)

                } ## end of Metropolis-Hastings move
            } ## end of else
        } ## end of moves in case larger than 1

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
            dfValue         <- c(varTilde$df[1:maxValue], zeroVector)
            wValue          <- c(varTilde$w[1:maxValue], zeroVector)
            dimValue        <- c(varTilde$dim[1:maxValue], zeroVector)
            aValue          <- c(varTilde$a[1:(maxValue + 1)], zeroVector)
        }

        ## Prepared list used by merge and split function
        kVal <- kValue
        listUpdate <- list(k       = kVal,
                            mu     = muValue[1:kVal],
                            sigmaf = sigmafValue[1:kVal],
                            sigmar = sigmarValue[1:kVal],
                            delta  = deltaValue[1:kVal],
                            df     = dfValue[1:kVal],
                            w      = wValue[1:kVal]
        )

        ## Run merging function if activated
        if (runMergeFunction) {
            listUpdate <- mergeNucleosomes(startPosForwardReads,
                                    startPosReverseReads, y, listUpdate,
                                    minInterval, maxInterval, minReads)
        }

        ## Run split function if activated
        if (runSplitFunction) {
            listUpdate <- splitNucleosome(startPosForwardReads,
                                    startPosReverseReads, y, listUpdate,
                                    minInterval, maxInterval, minReads)
        }

        ## Assign resulting values for this iteration
        kVal          <- listUpdate$k
        k[i]          <- kVal
        naVector      <- rep(NA, kMax - kVal)
        mu[i, ]       <- c(listUpdate$mu, naVector)
        sigmaf[i, ]   <- c(listUpdate$sigmaf, naVector)
        sigmar[i, ]   <- c(listUpdate$sigmar, naVector)
        delta[i, ]    <- c(listUpdate$delta, naVector)
        w[i, ]        <- c(listUpdate$w, naVector)
        df[i, ]       <- c(listUpdate$df, naVector)

    } ###end of boucle RJMCMC

    ## ATTENTION: Beware that the potential return of NA is not handled
    ## Getting the number of nucleosomes with the highest frequency
    km          <- elementWithHighestMode(as.integer(k))
    kPositions  <- which(as.integer(k) == km)

    mu_hat     <- colMeans(mu[kPositions, 1:km, drop = FALSE])
    sigmaf_hat <- colMeans(sigmaf[kPositions, 1:km, drop = FALSE])
    sigmar_hat <- colMeans(sigmar[kPositions, 1:km, drop = FALSE])
    w_hat      <- colMeans(w[kPositions, 1:km, drop = FALSE])
    delta_hat  <- colMeans(delta[kPositions, 1:km, drop = FALSE])
    df_hat     <- round(colMeans(df[kPositions, 1:km, drop = FALSE]))

    # Getting 2.5% and 97.5% quantiles for each important data type
    qmu     <- t(apply(mu[kPositions, 1:km, drop = FALSE], MARGIN = 2,
                        FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE))
    qsigmaf <- t(apply(sigmaf[kPositions, 1:km, drop = FALSE], MARGIN = 2,
                        FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE))
    qsigmar <- t(apply(sigmar[kPositions, 1:km, drop = FALSE], MARGIN = 2,
                        FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE))
    qdelta  <- t(apply(delta[kPositions, 1:km, drop = FALSE], MARGIN = 2,
                        FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE))
    qdf     <- t(apply(df[kPositions, 1:km, drop = FALSE], MARGIN = 2,
                        FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE))
    qw      <- t(apply(w[kPositions, 1:km, drop = FALSE], MARGIN = 2,
                        FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE))

    # Create the final list
    result <- list(
        call    = cl,
        K       = k,
        k       = km,
        mu      = mu_hat,
        sigmaf  = sigmaf_hat,
        sigmar  = sigmar_hat,
        delta   = delta_hat,
        df      = df_hat,
        w       = w_hat,
        qmu     = qmu,
        qsigmaf = qsigmaf,
        qsigmar = qsigmar,
        qdelta  = qdelta,
        qdf     = qdf,
        qw      = qw
    )

    class(result)<-"rjmcmcNucleosomes"

    return(result)
}
