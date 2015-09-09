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
#' @param kmax a positive \code{integer} or \code{numeric}, the maximum number
#' of nucleosomes per region. Non-integer values
#' of \code{kmax} will be casted to \code{integer} and truncated towards zero.
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
#'          nbrIterations = 1000, lambda = 2, kmax = 30,
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
                    nbrIterations, kmax, lambda,
                    minInterval, maxInterval, minReads)
{
    # Get call information
    cl <- match.call()

    # Parameters validation
    validateParameters(startPosForwardReads, startPosReverseReads,
                            nbrIterations, kmax, lambda, minInterval,
                            maxInterval, minReads)

    # Casting specific inputs as integer
    minReads        <- as.integer(minReads)
    nbrIterations   <- as.integer(nbrIterations)
    kmax            <- as.integer(kmax)

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

    paramValues = list(startPSF =  startPosForwardReads, startPSR = startPosReverseReads
                       , kmax = kmax, lambda = lambda, minReads = minReads
                       , y = y, nr = nr, nf =nf, nbrReads = nbrReads , zeta = zeta
                       , deltamin = deltamin, deltamax = deltamax
                       , minReadPos= minReadPos, maxReadPos = maxReadPos )

    # Vector of the number of nucleosomes (integer values)
    k               <- Rle(rep(0L, nbrIterations))

    # Vector of the position of the nucleosomes
    mu              <- matrix(0, nrow = nbrIterations, ncol = kmax)



    sigmaf          <- matrix(0, nrow = nbrIterations, ncol = kmax)

    sigmar          <- matrix(0, nrow = nbrIterations, ncol = kmax)

    delta           <- matrix(0, nrow = nbrIterations, ncol = kmax)

    w               <- matrix(0, nrow = nbrIterations, ncol = kmax)

    dl              <- matrix(0L, nrow = nbrIterations, ncol = kmax)

    # Note use at the end

#     a               <- matrix(0, nrow = nbrIterations, ncol = kmax + 1L)
#     dim             <- matrix(0L, nrow = nbrIterations, ncol = kmax)
#     rhob            <- rep(0, nbrIterations)
#     rhod            <- rep(0, nbrIterations)
#     rhomh           <- rep(0, nbrIterations)
#     Kn              <- rep(0, nbrIterations)



    # Tilde tmp
#     ktilde          <- rep(0L, nbrIterations)
#     mutilde         <- matrix(0, nrow = nbrIterations, ncol = kmax)
#     sigmaftilde     <- matrix(0, nrow = nbrIterations, ncol = kmax)
#     sigmartilde     <- matrix(0, nrow = nbrIterations, ncol = kmax)
#     deltatilde      <- matrix(0, nrow = nbrIterations, ncol = kmax)
#     wtilde          <- matrix(0, nrow = nbrIterations, ncol = kmax)
#     atilde          <- matrix(0, nrow = nbrIterations, ncol = kmax + 1L)
#     dimtilde        <- matrix(0, nrow = nbrIterations, ncol = kmax)
#     dltilde         <- matrix(3, nrow = nbrIterations, ncol = kmax)


    k[1]            <- 1L

    mu[1, 1]        <- runif(1, minReadPos, maxReadPos)

    sigmaf[1, 1]    <- 1
    sigmar[1, 1]    <- 1
    delta[1, 1]     <- runif(1, 0, 2*(mu[1,1] - minReadPos))
    w[1, 1]         <- 1
    dl[1, 1]        <- 3

#     a[1, 1]         <- minReadPos
#     a[1, as.integer(k[1]) + 1L]  <- maxReadPos
#
#     dim[1,1]        <- nbrReads


#     Kaf             <- matrix(0, nrow = nf, ncol = kmax)
#     Kbf             <- matrix(0, nrow = nf, ncol = kmax)
#     Kar             <- matrix(0, nrow = nr, ncol = kmax)
#     Kbr             <- matrix(0, nrow = nr, ncol = kmax)
#
#     Y1f             <- rep(0, nf)
#     Y2f             <- rep(0, nf)
#     Y1r             <- rep(0, nr)
#     Y2r             <- rep(0, nr)

    nbrIterations   <- ifelse(nbrReads <= 10, 1000, nbrIterations)

    kValue <- as.integer(k[1])

    muValue       <- mu[1,]
    sigmafValue   <- sigmaf[1,]
    sigmarValue   <- sigmar[1,]
    deltaValue    <- delta[1,]
    wValue        <-  w[1,]
    dlValue       <- dl[1,]
    aValue        <- rep(0, kmax)
    aValue[1, 1]  <- minReadPos
    aValue[1, as.integer(k[1]) + 1L]  <- maxReadPos
    dimValue       <- rep(0, kmax)
    dimValue[1, 1] <- nbrReads

    for (i in 2:nbrIterations) {

        ## Current number of nucleosomes
        ## kValue <- k[i-1]
        varTilde <- list()
        if (kValue == 1L) {
            ## CASE : Number of nucleosomes equal to 1
            u <- runif(1)

            if (u <= 0.5) {

                # Birth move in case k=1

                varTilde <- birthMoveK1(paramValues, kValue, muValue
                                        , sigmafValue, sigmarValue, deltaValue
                                        , wValue, dlValue, aValue, dimValue )
            } ### end of B move in case k=1
            else {

                ###Metropolis-Hastings move
                varTilde <- m-hMoveK1(paramValues, kValue, muValue
                                      , sigmafValue, sigmarValue, deltaValue
                                      , wValue, dlValue, aValue, dimValue )

            }    ### end of M-H move in case k=1
        }  ### end of test of k=1

        else {

            u<-runif(1)

            if (u <= Dk(kValue, lambda, kmax)) {

                ### Death move
                varTilde <- deathMove(paramValues, kValue, muValue
                                      , sigmafValue, sigmarValue, deltaValue
                                      , wValue, dlValue, aValue, dimValue )

            } else {

                if (u <= (Dk(kValue, lambda, kmax) + Bk(kValue, lambda, kmax))) {

                    #Birth move
                    varTilde <- birthMove(paramValues, kValue, muValue
                                          , sigmafValue, sigmarValue, deltaValue
                                          , wValue, dlValue, aValue, dimValue )

                } else {
                    ### Metropolis-Hastings move
                    varTilde <- m-hMove(paramValues, kValue, muValue
                                          , sigmafValue, sigmarValue, deltaValue
                                          , wValue, dlValue, aValue, dimValue )

                }
            } #end of else of B move and M-H move
        } ###end of moves in case k>=2

        v <- runif(1)      #Acceptation/rejet

        if (varTilde$rho >= v && varTilde$ktilde <= kmax) {
            kValue                    <- varTilde$k
            maxValue                <- as.integer(kValue)
            muValue[ 1:maxValue]       <- varTilde$mu[ 1:maxValue]
            sigmafValue[ 1:maxValue]   <- varTilde$sigmaf[ 1:maxValue]
            sigmarValue[ 1:maxValue]   <- varTilde$sigmar[ 1:maxValue]
            deltaValue[ 1:maxValue]    <- varTilde$delta[ 1:maxValue]
            dlValue[ 1:maxValue]       <- varTilde$dl[ 1:maxValue]
            wValue[ 1:maxValue]        <- varTilde$w[ 1:maxValue]
            dimValue[ 1:maxValue]       <- varTilde$dim[ 1:maxValue]
            aValue[ 1:(maxValue+1)]     <- varTilde$a[ 1:(maxValue + 1)]
        }

        # Si on veut faire le merge on peut le faire ici
        # et le conserver sans changer les varValue
        k[i]                    <- kValue
        maxValue                <- as.integer(k[i])
        mu[i, 1:maxValue]       <- muValue[ 1:maxValue]
        sigmaf[i, 1:maxValue]   <- sigmafValue[ 1:maxValue]
        sigmar[i, 1:maxValue]   <- sigmarValue[ 1:maxValue]
        delta[i, 1:maxValue]    <- deltaValue[ 1:maxValue]
        dl[i,1:maxValue]        <- dlValue[ 1:maxValue]
        w[i,1:maxValue]         <- wValue[ 1:maxValue]

        ## Set the new value of kValue for the next iteration

        #kValue <- as.integer(k[i])


    } ###end of boucle RJMCMC

    liste <- rep(list(NULL),nbrIterations)

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

        listeUpdate <- mergeNucleosomes(startPosForwardReads,
                                        startPosReverseReads,y, new.list,
                                        minInterval,
                                        maxInterval, minReads)
        kVal                <- listeUpdate$k
        k[i]                <- kVal
        mu[i, ]       <- c(listeUpdate$mu, rep(0, kmax - kVal))
        sigmaf[i, 1:kVal]   <- c(listeUpdate$sigmaf, rep(0, kmax - kVal))
        sigmar[i, 1:kVal]   <- c(listeUpdate$sigmar, rep(0, kmax - kVal))
        delta[i, 1:kVal]    <- c(listeUpdate$delta, rep(0, kmax - kVal))
        w[i, 1:kVal]        <- c(listeUpdate$w, rep(0, kmax - kVal))
        dl[i, 1:kVal]       <- c(listeUpdate$dl, rep(0, kmax - kVal))
    }

#     kmax        <- max(kmax, sapply(1:nbrIterations,
#                                         function(i){liste[[i]]$k}))
#     mu          <- matrix(0, nrow=nbrIterations, ncol=kmax)
#     sigmaf      <- matrix(0, nrow=nbrIterations, ncol=kmax)
#     sigmar      <- matrix(0, nrow=nbrIterations, ncol=kmax)
#     delta       <- matrix(0, nrow=nbrIterations, ncol=kmax)
#     w           <- matrix(0, nrow=nbrIterations, ncol=kmax)
#     dl          <- matrix(0, nrow=nbrIterations, ncol=kmax)
#
#     for (i in 1:nbrIterations)
#     {
#         kVal                <- liste[[i]]$k
#         k[i]                <- kVal
#         mu[i, 1:kVal]       <- liste[[i]]$mu
#         sigmaf[i, 1:kVal]   <- liste[[i]]$sigmaf
#         sigmar[i, 1:kVal]   <- liste[[i]]$sigmar
#         delta[i, 1:kVal]    <- liste[[i]]$delta
#         w[i, 1:kVal]        <- liste[[i]]$w
#         dl[i, 1:kVal]       <- liste[[i]]$dl
#     }

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
