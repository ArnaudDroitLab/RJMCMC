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
#' \item \code{qsigmaf} TODO
#' \item \code{qsigmar} TODO
#' \item \code{qdelta} TODO
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
#' @import BiocGenerics
#' @author Rawane Samb
#' @export
RJMCMC <- function(startPosForwardReads, startPosReverseReads,
                    nbrIterations, kmax, lambda,
                    minInterval, maxInterval, minReads)
{
    ## ASTRID : voir si minInterval, maxInterval
    ## ne pourraient pas etre des integers
    ## ASTRID : il faudrait aussi penser au nom des variables

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
    d <- c(rep(1, nf), rep(-1, nr))
    yOrder <- order(y)
    y <- y[yOrder]
    d <- d[yOrder]
    rm(yOrder)

    ## ASTRID : voir si zeta, detamin, deltamax devraient etre des integer
    zeta            <- 147
    deltamin        <- 142
    deltamax        <- 152

    # Vector of the number of nucleosomes (integer values)
    k               <- rep(0L, nbrIterations)
    ktilde          <- rep(0L, nbrIterations)

    # Max and min read positions
    minReadPos <- min(y)
    maxReadPos <- max(y)

    # Vector of the position of the nucleosomes
    mu              <- matrix(0, nrow = nbrIterations, ncol = kmax)
    mutilde         <- matrix(0, nrow = nbrIterations, ncol = kmax)

    sigmaftilde     <- matrix(0, nrow = nbrIterations, ncol = kmax)
    sigmaf          <- matrix(0, nrow = nbrIterations, ncol = kmax)
    sigmartilde     <- matrix(0, nrow = nbrIterations, ncol = kmax)
    sigmar          <- matrix(0, nrow = nbrIterations, ncol = kmax)
    deltatilde      <- matrix(0, nrow = nbrIterations, ncol = kmax)
    delta           <- matrix(0, nrow = nbrIterations, ncol = kmax)
    wtilde          <- matrix(0, nrow = nbrIterations, ncol = kmax)
    w               <- matrix(0, nrow = nbrIterations, ncol = kmax)
    a               <- matrix(0, nrow = nbrIterations, ncol = kmax + 1L)
    atilde          <- matrix(0, nrow = nbrIterations, ncol = kmax + 1L)
    dimtilde        <- matrix(0, nrow = nbrIterations, ncol = kmax)
    dim             <- matrix(0, nrow = nbrIterations, ncol = kmax)
    dl              <- matrix(0, nrow = nbrIterations, ncol = kmax)
    dltilde         <- matrix(3, nrow = nbrIterations, ncol = kmax)

    k[1]            <- 1L

    mu[1, 1]        <- runif(1, minReadPos, maxReadPos)

    sigmaf[1, 1]    <- 1
    sigmar[1, 1]    <- 1
    delta[1, 1]     <- runif(1, 0, 2*(mu[1,1] - minReadPos))
    w[1, 1]         <- 1
    dl[1, 1]        <- 3

    a[1, 1]         <- minReadPos
    a[1, k[1] + 1]  <- maxReadPos

    dim[1,1]        <- nbrReads

    rhob            <- rep(0, nbrIterations)
    rhod            <- rep(0, nbrIterations)
    rhomh           <- rep(0, nbrIterations)
    # Kn1             <- rep(0, nbrIterations)
    # Kn2             <- rep(0, nbrIterations)
    Kn              <- rep(0, nbrIterations)
    # Ln1             <- rep(0, nbrIterations)
    # Ln2             <- rep(0, nbrIterations)
    # Ln              <- rep(0, nbrIterations)

    Kaf             <- matrix(0, nrow = nf, ncol = kmax)
    Kbf             <- matrix(0, nrow = nf, ncol = kmax)
    Kar             <- matrix(0, nrow = nr, ncol = kmax)
    Kbr             <- matrix(0, nrow = nr, ncol = kmax)

    Y1f             <- rep(0, nf)
    Y2f             <- rep(0, nf)
    Y1r             <- rep(0, nr)
    Y2r             <- rep(0, nr)

    nbrIterations   <- ifelse(nbrReads <= 10, 1000, nbrIterations)

    for (i in 2:nbrIterations) {

        ## Current number of nucleosomes
        kValue <- k[i-1]

        if (kValue == 1L) {
            ## CASE : Number of nucleosomes equal to 1
            u <- runif(1)

            if (u <= 0.5) {
                ktilde[i] <- kValue + 1L
                count  <- 1L
                repeat {
                    j <- sample(1:kValue, 1)
                    mutilde[i, j] <- runif(1, minReadPos, mu[i-1, j])
                    mutilde[i, 1:ktilde[i]] <- sort(c(mu[i-1, 1:kValue],
                                                        mutilde[i, j]))

                    atilde[i, j+1] <- runif(1, mutilde[i, j], mutilde[i, j+1])
                    atilde[i, 1:(ktilde[i]+1)] <- sort(c(a[i-1, 1:ktilde[i]],
                                                            atilde[i, j+1]))
                    atilde[i, 1]                <- minReadPos
                    atilde[i, (ktilde[i]+1)]    <- maxReadPos

                    dimtilde[i, 1] <- length(y[atilde[i, 1] <= y &
                                                    y < atilde[i, 2]])
                    dimtilde[i, ktilde[i]] <- length(y[atilde[i, ktilde[i]] <=
                                                        y & y <= maxReadPos])
                    if (ktilde[i] > 2) {
                        for (m in 2:(ktilde[i]-1)) {
                            dimtilde[i,m] <- length(y[(atilde[i, m] <= y &
                                                    y < atilde[i, m+1])]) }
                    }
                    Pr <- min(dimtilde[i, 1:ktilde[i]])

                    ybar <- mean(y[atilde[i, j] <= y & y <= atilde[i, j+1]])
                    classesf <- y[atilde[i, j] <= y & y <= ybar]
                    classesr <- y[ybar <= y & y <= atilde[i, j+1]]

                    Lf <- length(classesf[!duplicated(classesf)])
                    Lr <- length(classesr[!duplicated(classesr)])


                    count <- count + 1L

                    if ((Pr > 1 & Lf > 1 & Lr > 1)  ||
                                            count == 1000L) break()
                }

                if (count == 1000L) {
                    rhob[i] <- 0
                } else  {

                    dltilde[i, j] <- sample(3:30, 1)

                    sigmaftilde[i, j] <- ifelse(Lf > 1,
                            var(classesf) * (dltilde[i, j] - 2)/dltilde[i, j],
                            sigmaf[i-1, j])
                    sigmartilde[i, j] <- ifelse(Lr > 1,
                            var(classesr) * (dltilde[i, j] - 2)/dltilde[i, j],
                            sigmar[i-1, j])

                    sigmaftilde[i, 1:ktilde[i]] <- c(sigmaf[i-1, 1:k[i-1]],
                                                        sigmaftilde[i, j])
                    sigmartilde[i, 1:ktilde[i]] <- c(sigmar[i-1, 1:k[i-1]],
                                                        sigmartilde[i, j])

                    deltatilde[i, j] <- tnormale(zeta,
                            1/(sigmaftilde[i, j]^{-1} + sigmartilde[i,j]^{-1}),
                            deltamin, deltamax)
                    deltatilde[i, 1:ktilde[i]] <- c(deltatilde[i, j],
                                                        delta[i-1, 1:kValue])

                    alpha                   <- rep(1, kValue)
                    alphatilde              <- rep(1, ktilde[i])
                    alphaproptilde          <- rep(1, ktilde[i])
                    alphaprop               <- rep(1, kValue)
                    wtilde[i, 1:ktilde[i]]  <- rdirichlet(1, alphaproptilde)
                    ennetilde               <- dimtilde[i, 1:ktilde[i]]
                    enne <- rmultinom(1, nbrReads, w[i-1, 1:kValue])

                    #Rapport de vraisemblance

                    for (m in 1:k[i-1]) {
                        Kaf[,m] <- (wtilde[i, m]*(1/sqrt(sigmaftilde[i,m]))*dt((startPosForwardReads-mutilde[i,m]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,m]),dltilde[i,m]))
                        Kbf[,m] <- (w[i-1, m]*(1/sqrt(sigmaf[i-1,m]))*dt((startPosForwardReads-mu[i-1,m]+delta[i-1,m]/2)/sqrt(sigmaf[i-1,m]),dl[i-1,m]))
                    }
                    Kaf[,ktilde[i]] <- (wtilde[i, ktilde[i]]*(1/sqrt(sigmaftilde[i,ktilde[i]]))*dt(( startPosForwardReads -mutilde[i,ktilde[i]]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,ktilde[i]]),dltilde[i,m]))
                    for (s in 1:nf) {
                        Y1f[s] <- log(sum(Kaf[s, 1:ktilde[i]]))
                        Y2f[s] <- log(sum(Kbf[s, 1:k[i - 1]]))
                    }

                    for (m in 1:k[i-1]) {
                        Kar[,m] <- (wtilde[i, m]*(1/sqrt(sigmartilde[i,m]))*dt((startPosReverseReads-mutilde[i,m]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,m]),dltilde[i,m]))
                        Kbr[,m] <- (w[i-1, m]*(1/sqrt(sigmar[i-1,m]))*dt((startPosReverseReads-mu[i-1,m]-delta[i-1,m]/2)/sqrt(sigmar[i-1,m]),dl[i-1,m]))
                    }
                    Kar[,ktilde[i]] <- (wtilde[i, ktilde[i]]*(1/sqrt(sigmartilde[i,ktilde[i]]))*dt(( startPosReverseReads -mutilde[i,ktilde[i]]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,ktilde[i]]),dltilde[i,m]))
                    for (s in 1:nr) {
                        Y1r[s] <- log(sum(Kar[s, 1:ktilde[i]]))
                        Y2r[s] <- log(sum(Kbr[s, 1:k[i - 1]]))
                    }

                    Kn <- sum(Y1f) + sum(Y1r)
                    Kd <- sum(Y2f) + sum(Y2r)

                    q <- Kn - Kd
                    rap.q <- exp(q)

                    rap.vrais <- rap.q

                    if (j == 1) {
                        qalloc <- 1/(mu[i-1, j] - minReadPos) # Density of mutilde[i,j]
                    } else {
                        qalloc <- 1/(mu[i-1, j] - mu[i-1, j-1])
                    }

                    rap.priormu <- (priorMuDensity(mutilde[i, 1:ktilde[i]],y)/priorMuDensity(mu[i-1,1:k[i-1]],y))
                    rap.priorw <- (ddirichlet(wtilde[i, 1:ktilde[i]],alphatilde)/ddirichlet(w[i-1,1:k[i-1]],alpha) )
                    rap.priorenne <- dmultinom(ennetilde, nbrReads,wtilde[i,1:ktilde[i]])/dmultinom(dim[i-1,1:k[i-1]], nbrReads,w[i-1,1:k[i-1]])
                    rap.priork <- (dpois(ktilde[i], lambda)/dpois(k[i - 1], lambda))
                    rap.propmu <- (1/(qalloc))
                    rap.propw <- (ddirichlet(w[i-1, 1:k[i-1]], alphaprop)/ddirichlet(wtilde[i, 1:ktilde[i]], alphaproptilde))

                    rap.prior <- rap.priormu * rap.priorw * rap.priorenne * rap.priork
                    rap.prop <-  rap.propmu  * rap.propw

                    rhob[i] <- min(1,(rap.vrais) * (rap.prior) * (rap.prop ) * (Dk(ktilde[i],lambda,kmax)/Bk(k[i - 1], lambda, kmax)))

                }

                rhob[i] <- ifelse (is.na(rhob[i]), 0, rhob[i])

                v <- runif(1)      #Acceptation/rejet du Birth move

                if (rhob[i] >= v && ktilde[i] <= kmax) {
                    k[i]                    <- ktilde[i]
                    maxValue                <- k[i]
                    mu[i, 1:maxValue]       <- mutilde[i, 1:maxValue]
                    sigmaf[i, 1:maxValue]   <- sigmaftilde[i, 1:maxValue]
                    sigmar[i, 1:maxValue]   <- sigmartilde[i, 1:maxValue]
                    delta[i, 1:maxValue]    <- deltatilde[i, 1:maxValue]
                    dl[i, 1:maxValue]       <- dltilde[i, 1:maxValue]
                    w[i, 1:maxValue]        <- wtilde[i, 1:maxValue]
                    dim[i,1:maxValue]       <- dimtilde[i, 1:maxValue]
                    a[i,1:(maxValue+1)]     <- atilde[i, 1:(maxValue + 1)]
                } else {
                    k[i]                    <- k[i-1]
                    maxValue                <- k[i]
                    mu[i, 1:maxValue]       <- mu[i-1, 1:maxValue]
                    sigmaf[i, 1:maxValue]   <- sigmaf[i-1, 1:maxValue]
                    sigmar[i, 1:maxValue]   <- sigmar[i-1, 1:maxValue]
                    delta[i, 1:maxValue]    <- delta[i-1, 1:maxValue]
                    dl[i,1:maxValue]        <- dl[i-1, 1:maxValue]
                    w[i,1:maxValue]         <- w[i-1, 1:maxValue]
                    dim[i,1:maxValue]       <- dim[i-1, 1:maxValue]
                    a[i,1:(maxValue + 1)]   <- a[i-1, 1:(maxValue+1)]
                }
            } ### end of B move in case k=1

            else {

                ###Metropolis-Hastings move

                ktilde[i] <- k[i-1]
                count     <- 1L
                repeat {
                    j <- sample(1:k[i-1], 1)
                    mutilde[i,j] <- runif(1, mu[i-1,j], maxReadPos)
                    mutilde[i,1:ktilde[i]] <- sort(c(mutilde[i,1:ktilde[i]]))

                    atilde[i,j]   <- minReadPos
                    atilde[i,j+1] <- maxReadPos

                    dimtilde[i,1] <- length(y[atilde[i,1]<=y & y< atilde[i,2]])
                    dimtilde[i,ktilde[i]] <- length(y[atilde[i,ktilde[i]] <= y & y <= maxReadPos])
                    if (ktilde[i]>2) {
                        for (m in 2: (ktilde[i]-1)) {
                            dimtilde[i,m] <-length(y[(atilde[i, m] <=y & y < atilde[i, m + 1])])}}
                    Pr <- min(dimtilde[i,1:ktilde[i]])

                    ybar <- mean(y[atilde[i,j]<=y & y<=atilde[i,j+1]])
                    classesf <- y[atilde[i,j]<=y & y<=ybar]
                    classesr <- y[ybar <= y & y <= atilde[i, j + 1]]

                    Lf <- length(classesf[!duplicated(classesf)])
                    Lr <- length(classesr[!duplicated(classesr)])
                    count <- count + 1L

                    if ( (Pr > 1 & Lf > 1 & Lr > 1)  ||
                                    count == 1000L) break()
                }

                if (count == 1000L) {
                    rhomh[i] <- 0
                } else {

                    sigmaftilde[i, 1:ktilde[i]] <- sigmaf[i-1, 1:k[i - 1]]
                    sigmartilde[i, 1:ktilde[i]] <- sigmar[i-1, 1:k[i - 1]]

                    dltilde[i, j] <- sample(3:30,1)

                    sigmaftilde[i, j] <- ifelse(Lf>1, var(classesf)*(dltilde[i,j]-2)/dltilde[i,j], sigmaf[i-1,j])
                    sigmartilde[i, j] <- ifelse(Lr>1, var(classesr)*(dltilde[i,j]-2)/dltilde[i,j], sigmar[i-1,j])

                    deltatilde[i, 1:ktilde[i]] <- delta[i-1,1:k[i-1]]

                    deltatilde[i, j] <- tnormale(zeta, 1/(sigmaftilde[i,j]^{-1}+sigmartilde[i,j]^{-1}), deltamin, deltamax)

                    alpha <- rep(1,k[i-1])
                    alphatilde <- rep(1,ktilde[i])
                    ennetilde <- dimtilde[i,1:ktilde[i]]
                    alphaproptilde <- rep(1,ktilde[i])
                    alphaprop <- rep(1,k[i-1])
                    wtilde[i, 1:ktilde[i]] <- rdirichlet(1,alphaproptilde)

                    ### calcul du rapport de vraisemblance de M-H move

                    for (m in 1:ktilde[i]) {
                        Kaf[,m] <- (wtilde[i, m]*(1/sqrt(sigmaftilde[i,m]))*dt(( startPosForwardReads -mutilde[i,m]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,m]),dltilde[i,m]))
                        Kbf[,m] <- (w[i-1, m]*(1/sqrt(sigmaf[i-1,m]))*dt(( startPosForwardReads -mu[i-1,m]+delta[i-1,m]/2)/sqrt(sigmaf[i-1,m]),dl[i-1,m]))
                    }
                    for (s in 1:nf) {
                        Y1f[s] <-log(sum(Kaf[s, 1:ktilde[i]]))
                        Y2f[s] <-log(sum(Kbf[s, 1:k[i - 1]]))
                    }

                    for (m in 1:ktilde[i]) {
                        Kar[,m] <- (wtilde[i, m]*(1/sqrt(sigmartilde[i,m]))*dt((startPosReverseReads -mutilde[i,m]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,m]),dltilde[i,m]))
                        Kbr[,m] <- (w[i-1, m]*(1/sqrt(sigmar[i-1,m]))*dt((startPosReverseReads -mu[i-1,m]-delta[i-1,m]/2)/sqrt(sigmar[i-1,m]),dl[i-1,m]))
                    }
                    for (s in 1:nr) {
                        Y1r[s] <- log(sum(Kar[s, 1:ktilde[i]]))
                        Y2r[s] <- log(sum(Kbr[s, 1:k[i-1]]))
                    }

                    Kn <- sum(Y1f) + sum(Y1r)
                    Kd <- sum(Y2f) + sum(Y2r)

                    q <- Kn - Kd
                    rap.q <- exp(q)

                    rap.vrais <- rap.q

                    rap.priormu <- (priorMuDensity(mutilde[i, 1:ktilde[i]], y)/priorMuDensity(mu[i-1, 1:k[i-1]], y))
                    rap.priorw <- (ddirichlet(wtilde[i,1:ktilde[i]],alphatilde)/ddirichlet(w[i-1,1:k[i-1]],alpha) )
                    rap.priorenne <- dmultinom(ennetilde, nbrReads,wtilde[i,1:ktilde[i]])/dmultinom(dim[i-1,1:k[i-1]], nbrReads,w[i-1,1:k[i-1]])
                    rap.priork <- 1
                    rap.propmu <- 1
                    rap.propw <- (ddirichlet(w[i-1, 1:k[i - 1]],alphaprop)/ddirichlet(wtilde[i, 1:ktilde[i]], alphaproptilde))

                    rap.prior <- rap.priormu * rap.priorw * rap.priorenne * rap.priork
                    rap.prop <-  rap.propmu  * rap.propw
                    rhomh[i] <- min(1, rap.vrais * (rap.prior)  *  (rap.prop))

                }

#                 rhomh[i] <- ifelse( is.na(rhomh[i]) == FALSE, rhomh[i], 0)
                rhomh[i] <- ifelse(is.na(rhomh[i]), 0, rhomh[i])

                v = runif(1)

                if (rhomh[i] >= v && ktilde[i] <= kmax) {
                    k[i]                    <- ktilde[i]
                    maxValue                <- as.integer(k[i])
                    mu[i, 1:maxValue]       <- mutilde[i, 1:maxValue]
                    sigmaf[i, 1:maxValue]   <- sigmaftilde[i, 1:maxValue]
                    sigmar[i, 1:maxValue]   <- sigmartilde[i, 1:maxValue]
                    delta[i, 1:maxValue]    <- deltatilde[i, 1:maxValue]
                    dl[i, 1:maxValue]       <- dltilde[i, 1:maxValue]
                    w[i, 1:maxValue]        <- wtilde[i, 1:maxValue]
                    dim[i, 1:maxValue]      <- dimtilde[i, 1:maxValue]
                    a[i, 1:(maxValue + 1)]  <- atilde[i, 1:(maxValue + 1)]
                } else {
                    k[i] <- k[i-1]
                    maxValue                <- as.integer(k[i])
                    mu[i, 1:maxValue]       <- mu[i-1, 1:maxValue]
                    sigmaf[i, 1:maxValue]   <- sigmaf[i-1, 1:maxValue]
                    sigmar[i, 1:maxValue]   <- sigmar[i-1, 1:maxValue]
                    delta[i, 1:maxValue]    <- delta[i-1, 1:maxValue]
                    dl[i, 1:maxValue]       <- dl[i-1, 1:maxValue]
                    w[i, 1:maxValue]        <- w[i-1, 1:maxValue]
                    dim[i, 1:maxValue]      <- dim[i-1, 1:maxValue]
                    a[i, 1:(maxValue + 1)]  <- a[i-1, 1:(maxValue + 1)]
                }
            }    ### end of M-H move in case k=1
        }  ### end of test of k=1

        else {

            u<-runif(1)

            if (u <= Dk(kValue, lambda, kmax)) {

                ### Death move

                ktilde[i] <- kValue - 1L
                count  <- 1L
                repeat {

                    j <- sample(1:k[i-1],1)
                    X <- mu[i-1,1:k[i-1]]
                    mutilde[i,1:ktilde[i]] <- sort(c(X[-j]))

                    Yf <- sigmaf[i-1,1:k[i-1]]
                    sigmaftilde[i,1:ktilde[i]] <- Yf[-j]

                    Yr <- sigmar[i-1,1:k[i-1]]
                    sigmartilde[i,1:ktilde[i]] <- Yr[-j]

                    Delta <- delta[i-1,1:k[i-1]]
                    deltatilde[i,1:ktilde[i]] <- Delta[-j]

                    Dl <- dl[i-1,1:k[i-1]]
                    dltilde[i,1:ktilde[i]] <- Dl[-j]

                    Z <- a[i-1,1:(k[i-1]+1)]
                    if (j==k[i-1]) {atilde[i,1:(ktilde[i]+1)] <- Z[-j]}
                    else { atilde[i,1:(ktilde[i]+1)] <- Z[-(j+1)] }
                    atilde[i,1]             <- minReadPos
                    atilde[i,ktilde[i]+1]   <- maxReadPos

                    dimtilde[i,1] <- length(y[atilde[i,1]<=y & y<atilde[i,2]])
                    dimtilde[i,ktilde[i]] <- length(y[atilde[i,ktilde[i]]<=y & y<= maxReadPos])
                    if (ktilde[i]>2) {
                        for (m in 2: (ktilde[i]-1)) {
                            dimtilde[i,m] <- length(y[atilde[i,m]<=y & y<atilde[i,m+1]]) }}
                    Pr <- min(dimtilde[i,1:ktilde[i]])

                    ybar <- mean(y[atilde[i,j]<=y & y<=atilde[i,j+1]])
                    classesf <- y[atilde[i,j]<=y & y<=ybar]
                    classesr <- y[ybar<=y & y<=atilde[i,j+1]]

                    Lf <- length(classesf[!duplicated(classesf)])
                    Lr <- length(classesr[!duplicated(classesr)])
                    count <- count + 1L

                    if ( (Pr>1 & Lf>1 & Lr>1)  ||
                                count == 1000L) break()

                }

                if (count == 1000L) {
                    rhod[i] <- 0
                } else {
                    alpha <- rep(1, k[i-1])
                    alphatilde <- rep(1, ktilde[i])
                    alphaproptilde <- rep(1, ktilde[i])
                    alphaprop <- rep(1, k[i-1])
                    ennetilde <- dimtilde[i,1:ktilde[i]]
                    enne <- rmultinom(1, nbrReads, w[i-1,1:k[i-1]])
                    wtilde[i,1:ktilde[i]] <- rdirichlet(1, alphaproptilde)

                    ### Rapport de vraisemblance ###

                    for (m in 1:ktilde[i]) {
                        Kaf[,m] <- (wtilde[i,m]*(1/sqrt(sigmaftilde[i,m]))*dt(( startPosForwardReads -mutilde[i,m]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,m]),dltilde[i,m]))
                        Kbf[,m] <- (w[i-1,m]*(1/sqrt(sigmaf[i-1,m]))*dt(( startPosForwardReads -mu[i-1,m]+delta[i-1,m]/2)/sqrt(sigmaf[i-1,m]),dl[i-1,m]))
                    }
                    Kbf[,k[i-1]] <- (w[i-1,k[i-1]]*(1/sqrt(sigmaf[i-1,k[i-1]]))*dt(( startPosForwardReads -mu[i-1,k[i-1]]+delta[i-1,m]/2)/sqrt(sigmaf[i-1,k[i-1]]),dltilde[i,m]))

                    for (s in 1:nf) {
                        Y1f[s] <- log(sum(Kaf[s,1:ktilde[i]]))
                        Y2f[s] <- log(sum(Kbf[s,1:k[i-1]]))
                    }

                    for (m in 1:ktilde[i]) {
                        Kar[,m] <- (wtilde[i,m]*(1/sqrt(sigmartilde[i,m]))*dt(( startPosReverseReads -mutilde[i,m]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,m]),dltilde[i,m]))
                        Kbr[,m] <- (w[i-1,m]*(1/sqrt(sigmar[i-1,m]))*dt(( startPosReverseReads -mu[i-1,m]-delta[i-1,m]/2)/sqrt(sigmar[i-1,m]),dl[i-1,m]))
                    }
                    Kbr[,k[i-1]] <- (w[i-1,k[i-1]]*(1/sqrt(sigmar[i-1,k[i-1]]))*dt(( startPosReverseReads -mu[i-1,k[i-1]]-delta[i-1,m]/2)/sqrt(sigmar[i-1,k[i-1]]),dl[i-1,m]))

                    for (s in 1:nr) {
                        Y1r[s] <- log(sum(Kar[s,1:ktilde[i]]))
                        Y2r[s] <- log(sum(Kbr[s,1:k[i-1]]))
                    }

                    Kn <- sum(Y1f) + sum(Y1r)
                    Kd <- sum(Y2f) + sum(Y2r)

                    q <- Kn - Kd
                    rap.q <- exp(q)

                    rap.vrais <- rap.q

                    # Density of mutilde[i,j]
                    if (j == 1) {
                        qalloc <- 1/(mu[i-1,j+1] - minReadPos)
                    } else {
                        if (j == k[i-1]) {
                            qalloc <- 1/(maxReadPos - mu[i-1, j-1])
                        } else {
                            qalloc <- 1/(mu[i-1, j+1] - mu[i-1, j-1])
                        }
                    }

                    rap.priormu <- (priorMuDensity(mutilde[i,1:ktilde[i]],y)/priorMuDensity(mu[i-1,1:k[i-1]],y))
                    rap.priorw <- (ddirichlet(wtilde[i,1:ktilde[i]],alphatilde)/ddirichlet(w[i-1,1:k[i-1]],alpha) )
                    rap.priorenne <- dmultinom(ennetilde, nbrReads,wtilde[i,1:ktilde[i]])/dmultinom(dim[i-1,1:k[i-1]], nbrReads,w[i-1,1:k[i-1]])
                    rap.priork <- (dpois(ktilde[i],lambda)/dpois(k[i-1],lambda))
                    rap.propmu <- (qalloc)
                    rap.propw <- (ddirichlet(wtilde[i,1:ktilde[i]],alphaproptilde)/ddirichlet(w[i-1,1:k[i-1]],alphaprop))

                    rap.prior <- rap.priormu * rap.priorw * rap.priorenne * rap.priork
                    rap.prop <- rap.propmu  * rap.propw

                    rhod[i] <- min(1,(rap.vrais) * (rap.prior)  *  (rap.prop) * (Bk(ktilde[i],lambda,kmax)/Dk(k[i-1],lambda,kmax)))

                }

                rhod[i] <-  ifelse( is.na(rhod[i]) == FALSE, rhod[i], 0)

                v <- runif(1)      #Acceptation/rejet du Death move

                if (rhod[i] >= v && ktilde[i] <= kmax) {
                    k[i]                    <- ktilde[i]
                    maxValue                <- as.integer(k[i])
                    mu[i, 1:maxValue]       <- mutilde[i, 1:maxValue]
                    sigmaf[i, 1:maxValue]   <- sigmaftilde[i, 1:maxValue]
                    sigmar[i, 1:maxValue]   <- sigmartilde[i, 1:maxValue]
                    delta[i, 1:maxValue]    <- deltatilde[i, 1:maxValue]
                    dl[i, 1:maxValue]       <- dltilde[i, 1:maxValue]
                    w[i, 1:maxValue]        <- wtilde[i, 1:maxValue]
                    dim[i, 1:maxValue]      <- dimtilde[i, 1:maxValue]
                    a[i, 1:(maxValue + 1)]  <- atilde[i, 1:(maxValue + 1)]
                } else {
                    k[i]                    <- k[i-1]
                    maxValue                <- as.integer(k[i])
                    mu[i, 1:maxValue]       <- mu[i-1, 1:maxValue]
                    sigmaf[i, 1:maxValue]   <- sigmaf[i-1, 1:maxValue]
                    sigmar[i, 1:maxValue]   <- sigmar[i-1, 1:maxValue]
                    delta[i, 1:maxValue]    <- delta[i-1, 1:maxValue]
                    dl[i, 1:maxValue]       <- dl[i-1, 1:maxValue]
                    w[i, 1:maxValue]        <- w[i-1, 1:maxValue]
                    dim[i, 1:maxValue]      <- dim[i-1, 1:maxValue]
                    a[i, 1:(maxValue + 1)]  <- a[i-1, 1:(maxValue + 1)]
                }
            }

            else {  if (u <= (Dk(k[i-1], lambda, kmax) + Bk(k[i-1], lambda, kmax))) {

                #Birth move

                ktilde[i] <- k[i-1] + 1L
                count <- 1L
                repeat {
                    j <- sample(1:k[i-1],1)
                    if (j == 1) {mutilde[i,j] <- runif(1, minReadPos, mu[i-1,j])}
                    else {
                        mutilde[i,j] <- runif(1,mu[i-1,j-1],mu[i-1,j])}
                    mutilde[i,1:ktilde[i]] <- sort(c(mu[i-1,1:k[i-1]],mutilde[i,j]))

                    atilde[i,j+1] <- ifelse(j<ktilde[i-1], runif(1,mutilde[i,j],mutilde[i,j+1]), runif(1,mutilde[i,j],maxReadPos))
                    atilde[i,1:(ktilde[i]+1)] <- sort(c(a[i-1,1:ktilde[i]],atilde[i,j+1]))
                    if (j == 1){
                        atilde[i,j] <- minReadPos}
                    else {
                        atilde[i,j] <- runif(1,mutilde[i,j-1],mutilde[i,j])}
                    atilde[i,1]             <- minReadPos
                    atilde[i,ktilde[i]+1]   <- maxReadPos

                    dimtilde[i,1] <- length(y[atilde[i,1]<=y & y<atilde[i,2]])
                    dimtilde[i,ktilde[i]] <- length(y[atilde[i,ktilde[i]]<=y & y<=maxReadPos])
                    if (ktilde[i] > 2) {
                        for (m in 2: (ktilde[i]-1)) {
                            dimtilde[i,m] <- length(y[atilde[i,m]<=y & y<atilde[i,m+1]])
                        }
                    }
                    Pr <- min(dimtilde[i,1:ktilde[i]])

                    ybar <- mean(y[atilde[i,j] <= y & y <= atilde[i,j+1]])
                    classesf <- y[atilde[i,j] <= y & y <= ybar]
                    classesr <- y[ybar <= y & y <= atilde[i,j+1]]

                    Lf <- length(classesf[!duplicated(classesf)])
                    Lr <- length(classesr[!duplicated(classesr)])
                    count <- count + 1L
                    if ( (Pr>1 & Lf>1 & Lr>1)  ||
                                    count == 1000L) break()

                }

                if (count == 1000L) {
                    rhob[i] <- 0
                } else {
                    dltilde[i,j] <- sample(3:30,1)
                    sigmaftilde[i,j] <- ifelse(Lf > 1, var(classesf)*(dltilde[i,j]-2)/dltilde[i,j], sigmaf[i-1,j])
                    sigmartilde[i,j] <- ifelse(Lr > 1, var(classesr)*(dltilde[i,j]-2)/dltilde[i,j], sigmar[i-1,j])

                    if (j == 1) {
                        sigmaftilde[i,1:ktilde[i]]<- c(sigmaf[i-1,1:k[i-1]],sigmaftilde[i,j])
                        sigmartilde[i,1:ktilde[i]]<- c(sigmar[i-1,1:k[i-1]],sigmartilde[i,j] ) }
                    else {
                        sigmaftilde[i,1:ktilde[i]] <- c(sigmaf[i-1,1:(j-1)],sigmaftilde[i,j],sigmaf[i-1,j:k[i-1]])
                        sigmartilde[i,1:ktilde[i]] <- c(sigmar[i-1,1:(j-1)],sigmartilde[i,j],sigmar[i-1,j:k[i-1]]) }

                    deltatilde[i,j] <- tnormale(zeta, 1/(sigmaftilde[i,j]^{-1}+sigmartilde[i,j]^{-1}), deltamin, deltamax)
                    if (j == 1) {
                        deltatilde[i, 1:ktilde[i]] <- c(deltatilde[i,j], delta[i-1,1:k[i-1]]) }
                    else if (j == k[i-1]) {
                        deltatilde[i, 1:ktilde[i]] <- c(delta[i-1,1:k[i-1]], deltatilde[i,j]) }
                    else {
                        deltatilde[i, 1:ktilde[i]] <- c(delta[i-1,1:(j-1)], deltatilde[i,j], delta[i-1,j:k[i-1]]) }

                    alpha <- rep(1, k[i-1])
                    alphatilde <- rep(1, ktilde[i])
                    alphaproptilde <- rep(1, ktilde[i])
                    alphaprop <- rep(1, k[i-1])
                    ennetilde <- dimtilde[i, 1:ktilde[i]]
                    enne <- rmultinom(1, nbrReads, w[i-1, 1:k[i-1]])
                    wtilde[i,1:ktilde[i]] <- rdirichlet(1, alphaproptilde)

                    ### Rapport de vraisemblance ###

                    for (m in 1:k[i-1]) {
                        Kaf[,m]<- (wtilde[i,m]*(1/sqrt(sigmaftilde[i,m]))*dt(( startPosForwardReads -mutilde[i,m]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,m]),dltilde[i,m]))
                        Kbf[,m]<- (w[i-1,m]*(1/sqrt(sigmaf[i-1,m]))*dt(( startPosForwardReads -mu[i-1,m]+delta[i-1,m]/2)/sqrt(sigmaf[i-1,m]),dl[i-1,m]))
                    }
                    Kaf[,ktilde[i]]<- (wtilde[i,ktilde[i]]*(1/sqrt(sigmaftilde[i,ktilde[i]]))*dt(( startPosForwardReads -mutilde[i,ktilde[i]]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,ktilde[i]]),dltilde[i,m]))
                    for (s in 1:nf) {
                        Y1f[s] <-log(sum(Kaf[s,1:ktilde[i]]))
                        Y2f[s] <-log(sum(Kbf[s,1:k[i-1]]))
                    }

                    for (m in 1:k[i-1]) {
                        Kar[,m]<- (wtilde[i,m]*(1/sqrt(sigmartilde[i,m]))*dt(( startPosReverseReads -mutilde[i,m]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,m]),dltilde[i,m]))
                        Kbr[,m]<- (w[i-1,m]*(1/sqrt(sigmar[i-1,m]))*dt(( startPosReverseReads -mu[i-1,m]-delta[i-1,m]/2)/sqrt(sigmar[i-1,m]),dl[i-1,m]))
                    }
                    Kar[,ktilde[i]]<- (wtilde[i,ktilde[i]]*(1/sqrt(sigmartilde[i,ktilde[i]]))*dt(( startPosReverseReads -mutilde[i,ktilde[i]]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,ktilde[i]]),dltilde[i,m]))
                    for (s in 1:nr) {
                        Y1r[s] <-log(sum(Kar[s,1:ktilde[i]]))
                        Y2r[s] <-log(sum(Kbr[s,1:k[i-1]]))
                    }

                    Kn <- sum(Y1f)+sum(Y1r)
                    Kd <- sum(Y2f)+sum(Y2r)

                    q <- Kn - Kd
                    rap.q<-exp(q)

                    rap.vrais <- rap.q

                    #Density of
                    mutilde[i,j]
                    if (j==1) {
                        qalloc <- 1/(mu[i-1, j] - minReadPos)
                    } else {
                        qalloc <- 1/(mu[i-1, j] - mu[i-1, j-1])
                    }

                    rap.priormu <- (priorMuDensity(mutilde[i,1:ktilde[i]],y)/priorMuDensity(mu[i-1,1:k[i-1]],y))
                    rap.priorw <- (ddirichlet(wtilde[i,1:ktilde[i]],alphatilde)/ddirichlet(w[i-1,1:k[i-1]],alpha) )
                    rap.priorenne <- dmultinom(ennetilde, nbrReads,wtilde[i,1:ktilde[i]])/dmultinom(dim[i-1,1:k[i-1]], nbrReads,w[i-1,1:k[i-1]])
                    rap.priork <- (dpois(ktilde[i],lambda)/dpois(k[i-1],lambda))
                    rap.propmu <- (1/(qalloc))
                    rap.propw <- (ddirichlet(w[i-1,1:k[i-1]],alphaprop)/ddirichlet(wtilde[i,1:ktilde[i]],alphaproptilde))

                    rap.prior <- rap.priormu * rap.priorw * rap.priorenne * rap.priork
                    rap.prop <-  rap.propmu  * rap.propw
                    rhob[i] <- min(1, (rap.vrais) * (rap.prior)  *  (rap.prop ) * (Dk(ktilde[i],lambda,kmax)/Bk(k[i-1],lambda,kmax)))
                }

                v <- runif(1)      # Acceptation/rejet du Birth move

                rhob[i] <- ifelse(is.na(rhob[i]), 0, rhob[i])

                if (rhob[i] >= v && ktilde[i] <= kmax) {
                    k[i]                    <- ktilde[i]
                    maxValue                <- k[i]
                    mu[i, 1:maxValue]       <- mutilde[i, 1:maxValue]
                    sigmaf[i, 1:maxValue]   <- sigmaftilde[i, 1:maxValue]
                    sigmar[i, 1:maxValue]   <- sigmartilde[i, 1:maxValue]
                    delta[i, 1:maxValue]    <- deltatilde[i, 1:maxValue]
                    dl[i, 1:maxValue]       <- dltilde[i, 1:maxValue]
                    w[i, 1:maxValue]        <- wtilde[i, 1:maxValue]
                    dim[i, 1:maxValue]      <- dimtilde[i, 1:maxValue]
                    a[i, 1:(maxValue + 1)]  <- atilde[i, 1:(maxValue + 1)]
                } else {
                    k[i]                    <- k[i-1]
                    maxValue                <- k[i]
                    mu[i, 1:maxValue]       <- mu[i-1, 1:maxValue]
                    sigmaf[i, 1:maxValue]   <- sigmaf[i-1, 1:maxValue]
                    sigmar[i, 1:maxValue]   <- sigmar[i-1, 1:maxValue]
                    delta[i, 1:maxValue]    <- delta[i-1, 1:maxValue]
                    dl[i, 1:maxValue]       <- dl[i-1, 1:maxValue]
                    w[i, 1:maxValue]        <- w[i-1, 1:maxValue]
                    dim[i,1:maxValue]       <- dim[i-1,1:maxValue]
                    a[i,1:(maxValue + 1)]   <- a[i-1, 1:(maxValue + 1)]
                }
            }
            else {

                ### Metropolis-Hastings move

                ktilde[i] <- k[i-1]
                count  <- 1L
                repeat {
                    j <- sample(2:k[i-1],1)
                    mutilde[i,1:ktilde[i]] <- mu[i-1,1:k[i-1]]
                    if (j==1) {
                        mutilde[i,j] <- runif(1, minReadPos, mu[i-1,j+1])}
                    else {if (j==ktilde[i]) {
                        mutilde[i,j] <- runif(1,mu[i-1,j-1],maxReadPos)}
                        else { mutilde[i,j] <- runif(1,mu[i-1,j-1],mu[i-1,j+1]) }}
                    mutilde[i,1:ktilde[i]] <- sort(c(mutilde[i,1:ktilde[i]]))

                    atilde[i,1:(ktilde[i]+1)] <- sort(c(a[i-1,1:(k[i-1]+1)]))
                    if (j==ktilde[i]) {
                        atilde[i,j] <- runif(1,mutilde[i,j],maxReadPos)
                        atilde[i,j+1] <- maxReadPos }
                    else { if (j==1) {
                        atilde[i,j] <- minReadPos
                        atilde[i,j+1] <- runif(1,mutilde[i,j],mutilde[i,j+1])}
                        else {
                            atilde[i,j] <- runif(1,mutilde[i,j-1],mutilde[i,j])
                            atilde[i,j+1] <- runif(1,mutilde[i,j],mutilde[i,j+1])  }}
                    atilde[i,1] <- minReadPos
                    atilde[i,ktilde[i]+1] <- maxReadPos

                    dimtilde[i,1] <- length(y[atilde[i,1]<=y & y<atilde[i,2]])
                    dimtilde[i,ktilde[i]] <- length(y[atilde[i,ktilde[i]]<=y & y<=maxReadPos])
                    if (ktilde[i]>2) {
                        for (m in 2: (ktilde[i]-1)) {
                            dimtilde[i,m]<-length(y[atilde[i,m]<=y & y<atilde[i,m+1]]) }}
                    Pr <- min(dimtilde[i,1:ktilde[i]])

                    ybar <- mean(y[atilde[i,j]<=y & y<=atilde[i,j+1]])
                    classesf <- y[atilde[i,j]<=y & y<=ybar]
                    classesr <- y[ybar<=y & y<=atilde[i,j+1]]

                    Lf <- length(classesf[!duplicated(classesf)])
                    Lr <- length(classesr[!duplicated(classesr)])
                    count <- count + 1L
                    if ( (Pr > 1 & Lf > 1 & Lr > 1)
                                    || count == 1000L) break()
                }

                if (count == 1000L) {
                    rhomh[i] <- 0
                } else {

                    sigmaftilde[i,1:ktilde[i]] <- sigmaf[i-1,1:k[i-1]]
                    sigmartilde[i,1:ktilde[i]] <- sigmar[i-1,1:k[i-1]]

                    dltilde[i,j]=sample(3:30,1)

                    sigmaftilde[i,j] <- ifelse(Lf>1, var(classesf)*(dltilde[i,j]-2)/dltilde[i,j], sigmaf[i-1,j])
                    sigmartilde[i,j] <- ifelse(Lr>1, var(classesr)*(dltilde[i,j]-2)/dltilde[i,j], sigmar[i-1,j])

                    deltatilde[i,1:ktilde[i]] <- delta[i-1,1:k[i-1]]
                    deltatilde[i,j] <- tnormale(zeta, 1/(sigmaftilde[i,j]^{-1}+sigmartilde[i,j]^{-1}), deltamin, deltamax)

                    alpha <- rep(1,k[i-1])
                    alphatilde <- rep(1,ktilde[i])
                    ennetilde <- dimtilde[i,1:ktilde[i]]
                    alphaproptilde <- rep(1,ktilde[i])
                    alphaprop <- rep(1,k[i-1])
                    wtilde[i,1:ktilde[i]] <- rdirichlet(1,alphaproptilde)

                    ###calcul du rapport de vraisemblance de M-H move

                    for (m in 1:ktilde[i]) {
                        Kaf[,m] <- (wtilde[i,m]*(1/sqrt(sigmaftilde[i,m]))*dt(( startPosForwardReads -mutilde[i,m]+deltatilde[i,m]/2)/sqrt(sigmaftilde[i,m]),dltilde[i,m]))
                        Kbf[,m] <- (w[i-1,m]*(1/sqrt(sigmaf[i-1,m]))*dt(( startPosForwardReads -mu[i-1,m]+delta[i-1,m]/2)/sqrt(sigmaf[i-1,m]),dl[i-1,m]))
                    }
                    for (s in 1:nf) {
                        Y1f[s] <- log(sum(Kaf[s,1:ktilde[i]]))
                        Y2f[s] <- log(sum(Kbf[s,1:k[i-1]]))
                    }

                    for (m in 1:ktilde[i]) {
                        Kar[,m] <- (wtilde[i,m]*(1/sqrt(sigmartilde[i,m]))*dt(( startPosReverseReads -mutilde[i,m]-deltatilde[i,m]/2)/sqrt(sigmartilde[i,m]),dltilde[i,m]))
                        Kbr[,m] <- (w[i-1,m]*(1/sqrt(sigmar[i-1,m]))*dt(( startPosReverseReads -mu[i-1,m]-delta[i-1,m]/2)/sqrt(sigmar[i-1,m]),dl[i-1,m]))
                    }
                    for (s in 1:nr) {
                        Y1r[s] <- log(sum(Kar[s,1:ktilde[i]]))
                        Y2r[s] <- log(sum(Kbr[s,1:k[i-1]]))
                    }

                    Kn <- sum(Y1f) + sum(Y1r)
                    Kd <- sum(Y2f) + sum(Y2r)

                    q <-Kn - Kd
                    rap.q<-exp(q)

                    rap.vrais <- rap.q

                    rap.priormu <- (priorMuDensity(mutilde[i,1:ktilde[i]],y)/priorMuDensity(mu[i-1,1:k[i-1]],y))
                    rap.priorw <- (ddirichlet(wtilde[i,1:ktilde[i]], alphatilde)/ddirichlet(w[i-1,1:k[i-1]],alpha))
                    rap.priorenne <- dmultinom(ennetilde, nbrReads, wtilde[i,1:ktilde[i]])/dmultinom(dim[i-1,1:k[i-1]], nbrReads,w[i-1,1:k[i-1]])
                    rap.priork <- (1)
                    rap.propmu <- (1)
                    rap.propw <- (ddirichlet(w[i-1,1:k[i-1]],alphaprop)/ddirichlet(wtilde[i,1:ktilde[i]],alphaproptilde))
                    rap.prior <- rap.priormu * rap.priorw * rap.priorenne * rap.priork
                    rap.prop <-  rap.propmu  * rap.propw
                    rhomh[i] <- min(1, rap.vrais * (rap.prior)  *  (rap.prop ))

                }

                rhomh[i] <- ifelse ( is.na(rhomh[i])==FALSE, rhomh[i], 0)

                v <- runif(1)

                if (rhomh[i] >= v  ) {
                    k[i]                    <- ktilde[i]
                    maxValue                <- as.integer(k[i])
                    mu[i, 1:maxValue]       <- mutilde[i, 1:maxValue]
                    sigmaf[i, 1:maxValue]   <- sigmaftilde[i, 1:maxValue]
                    sigmar[i, 1:maxValue]   <- sigmartilde[i, 1:maxValue]
                    delta[i, 1:maxValue]    <- deltatilde[i, 1:maxValue]
                    dl[i, 1:maxValue]       <- dltilde[i, 1:maxValue]
                    w[i, 1:maxValue]        <- wtilde[i, 1:maxValue]
                    dim[i, 1:maxValue]      <- dimtilde[i, 1:maxValue]
                    a[i, 1:(maxValue + 1)]  <- atilde[i, 1:(maxValue + 1)]
                } else {
                    k[i]                    <- k[i-1]
                    maxValue                <- as.integer(k[i])
                    mu[i, 1:maxValue]       <- mu[i-1, 1:maxValue]
                    sigmaf[i, 1:maxValue]   <- sigmaf[i-1, 1:maxValue]
                    sigmar[i, 1:maxValue]   <- sigmar[i-1, 1:maxValue]
                    delta[i, 1:maxValue]    <- delta[i-1, 1:maxValue]
                    dl[i, 1:maxValue]       <- dl[i-1, 1:maxValue]
                    w[i, 1:maxValue]        <- w[i-1, 1:maxValue]
                    dim[i, 1:maxValue]      <- dim[i-1, 1:maxValue]
                    a[i, 1:(maxValue + 1)]  <- a[i-1, 1:(maxValue + 1)]
                }
            }
            } #end of else of B move and M-H move
        } ###end of moves in case k>=2
    } ###end of boucle RJMCMC

    liste <- rep(list(NULL),nbrIterations)

    for (i in 1:nbrIterations)
    {
        kVal <- k[i]
        new.list <- list(
                        k      = kVal,
                        mu     = mu[i, 1:kVal],
                        sigmaf = sigmaf[i, 1:kVal],
                        sigmar = sigmar[i, 1:kVal],
                        delta  = delta[i, 1:kVal],
                        dl     = dl[i, 1:kVal],
                        w      = w[i, 1:kVal]
        )

        liste[[i]] <- mergeNucleosomes(startPosForwardReads,
                                        startPosReverseReads,y, new.list,
                                        minInterval,
                                        maxInterval, minReads)
    }

    kmax        <- max(kmax, sapply(1:nbrIterations,
                                        function(i){liste[[i]]$k}))
    mu          <- matrix(0, nrow=nbrIterations, ncol=kmax)
    sigmaf      <- matrix(0, nrow=nbrIterations, ncol=kmax)
    sigmar      <- matrix(0, nrow=nbrIterations, ncol=kmax)
    delta       <- matrix(0, nrow=nbrIterations, ncol=kmax)
    w           <- matrix(0, nrow=nbrIterations, ncol=kmax)
    dl          <- matrix(0, nrow=nbrIterations, ncol=kmax)

    for (i in 1:nbrIterations)
    {
        kVal                <- liste[[i]]$k
        k[i]                <- kVal
        mu[i, 1:kVal]       <- liste[[i]]$mu
        sigmaf[i, 1:kVal]   <- liste[[i]]$sigmaf
        sigmar[i, 1:kVal]   <- liste[[i]]$sigmar
        delta[i, 1:kVal]    <- liste[[i]]$delta
        w[i, 1:kVal]        <- liste[[i]]$w
        dl[i, 1:kVal]       <- liste[[i]]$dl
    }

    ## Astrid : Si la fonction mode() retourne NA, le cas n'est pas gere
    ## Getting the number of nucleosomes with the highest frequency
    km <- elementWithHighestMode(k)
    kPositions  <- which(k == km)

    mu_hat     <- colMeans(mu[kPositions, 1:km, drop=FALSE])
    sigmaf_hat <- colMeans(sigmaf[kPositions, 1:km, drop=FALSE])
    sigmar_hat <- colMeans(sigmar[kPositions, 1:km, drop=FALSE])
    w_hat      <- colMeans(w[kPositions, 1:km, drop=FALSE])
    delta_hat  <- colMeans(delta[kPositions, 1:km, drop=FALSE])
    dl_hat     <- round(colMeans(dl[kPositions, 1:km, drop=FALSE]))

    # Getting 2.5% and 97.5% quantiles for each important data type
    qmu     <- t(apply(mu[,1:km, drop=FALSE], MARGIN=2, FUN=quantile,
                        probs=c(0.025, 0.975)))
    qsigmaf <- t(apply(sigmaf[,1:km, drop=FALSE], MARGIN=2, FUN=quantile,
                        probs=c(0.025, 0.975)))
    qsigmar <- t(apply(sigmar[,1:km, drop=FALSE], MARGIN=2, FUN=quantile,
                        probs=c(0.025, 0.975)))
    qdelta  <- t(apply(delta[,1:km, drop=FALSE], MARGIN=2, FUN=quantile,
                        probs=c(0.025, 0.975)))
    qdl     <- t(apply(dl[,1:km, drop=FALSE], MARGIN=2, FUN=quantile,
                        probs=c(0.025, 0.975)))
    qw      <- t(apply(w[,1:km, drop=FALSE], MARGIN=2, FUN=quantile,
                        probs=c(0.025, 0.975)))

    # Create the final list
    result <- list(
        call    = cl,
        K       = k,
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
