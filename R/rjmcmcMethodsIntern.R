#' @title Death Submove Probability
#'
#' @description Calculation of the death submove probability of a randomly
#' selected nucleosome using a truncated Poisson distribution.
#'
#' @param k a positive \code{integer}, the number of nucleosomes.
#'
#' @param lambda a positive \code{numeric}, the theorical mean
#' of the Poisson distribution.
#'
#' @param kMax a positive \code{integer}, the maximum number of nucleosomes
#' authorized. When \code{k} is equal or superior to \code{kMax}, the
#' returned value is \code{0}. Default: \code{30}.
#'
#' @return a \code{numeric} value representing the calculated death submove
#' probability. The value \code{0} when \code{k} is equal or superior to
#' \code{kMax} or when \code{k} is equal to \code{1}.
#'
#' @examples
#'
#' ## Return the death submove probability
#' rjmcmc:::Dk(k = 12L, lambda = 2L, kMax = 30L)
#'
#' ## Zero is returned when k = 1
#' rjmcmc:::Dk(k = 1L, lambda = 3L, kMax = 30L)
#'
#' ## Zerio is returned when k is superior to kMax
#' rjmcmc:::Dk(k = 31L, lambda = 2L, kMax = 30L)
#'
#' @author Rawane Samb
#' @importFrom stats dpois
#' @keywords internal
Dk <- function(k, lambda, kMax = 30) {
    ifelse((k == 1 || k > kMax), 0,
            0.5*min(1, dpois(k - 1, lambda)/dpois(k, lambda)))
}


#' @title Birth Submove Probability
#'
#' @description Calculation of the birth submove probability of adding a new
#' nucleosome using a truncated Poisson distribution.
#'
#' @param k a positive \code{integer}, the number of nucleosomes.
#'
#' @param lambda a positive \code{numeric}, the theorical mean
#' of the Poisson distribution.
#'
#' @param kMax a positive \code{integer}, the maximum number of nucleosomes
#' authorized. When \code{k} is equal or superior to \code{kMax}, the
#' returned value is \code{0}. Default: \code{30}.
#'
#' @return a \code{numeric} value. The value \code{0} when \code{k} is equal
#' or superior to \code{kMax} or when \code{k} is equal to \code{1}.
#'
#' @examples
#'
#' ## Return the birth submove probability
#' rjmcmc:::Bk(k = 14L, lambda = 1L, kMax = 30L)
#'
#' ## Zero is returned when k = 1
#' rjmcmc:::Bk(k = 1L, lambda = 3L, kMax = 30L)
#'
#' ## Zero is returned when k is superior to kMax
#' rjmcmc:::Bk(k = 31L, lambda = 2L, kMax = 30L)
#'
#' @author Rawane Samb
#' @importFrom stats dpois
#' @keywords internal
Bk <- function(k, lambda, kMax = 30) {
    ifelse((k >= kMax), 0,
           0.5 * min(1, dpois(k + 1, lambda) / dpois(k, lambda)))
}


#' @title Random deviate from a truncated normal distribution
#'
#' @description Generate a random deviate value from a normal distribution.
#' The returned value is included inside a specified range ]minValue,maxValue[
#' specified by user. The mean and variance of the normal distribution is
#' also specified by user.
#'
#' @param mu a \code{numeric} value, the mean of the normal distribution.
#'
#' @param sigma a non-negative \code{numeric}, the variance of the normal
#' distribution.
#'
#' @param minValue a \code{numeric} value, the inferior boundary of the
#' range in which the output value must be located. The output value has to be
#' superior to \code{minValue}.
#'
#' @param maxValue a \code{numeric} value, the superior boundary of the range
#' in which the output value must be located. The output value has to be
#' inferior to \code{maxValue}.
#'
#' @return a \code{numeric} which is superior to \code{minValue} and inferior
#' to \code{maxValue}.
#'
#' @author Rawane Samb
#' @importFrom stats rnorm
#' @keywords internal
tnormale <- function(mu, sigma, minValue, maxValue) {
    ## TODO : voir si on ne peut pas optimiser en créant un vecteur de valeurs
    ## Astrid : Est-ce qu'on doit ajouter une contrainte sur le nombre de
    ## boucles maximum
    repeat {
        y <- rnorm(1, mu, sd = sqrt(sigma))
        if (y > minValue & y < maxValue) break()
    }
    return(y)
}


#' @title Student Mixture Model
#'
#' @description Generation a value from a Student Mixture distribution.
#'
#' @param i a \code{integer}, a count parameter.
#'
#' @param k a positive \code{integer} value, the number of nucleosomes in the
#' analyzed region.
#'
#' @param weight a \code{vector} of positive \code{numerical} of length
#' \code{k}, the weight for each
#' nucleosome. The sum of all \code{weight} values must be equal to \code{1}.
#' The length of \code{weight} must be equal to \code{k}.
#'
#' @param mu a \code{vector} of positive \code{integer} of length \code{k},
#' the positions of all the nucleosomes in the analyzed region. The length
#' of \code{weight} must be equal to \code{k}.
#'
#' @param sigma a \code{vector} of \code{numeric} of length \code{k}, the
#' variance for each nucleosome. The length of \code{sigma} must be equal
#' to \code{k}.
#'
#' @param dfr a \code{vector} of \code{numeric} of length \code{k}, the degree
#' of freedom for each nucleosome. The length of \code{dfr} must be equal
#' to \code{k}.
#'
#' @return a \code{numerical}, the value generated from a Student Mixture
#' distribution.
#'
#' @examples
#'
#' ## Return a value generated from a student mixture
#' rjmcmc:::student.mixture(i = 1L, k = 4L, weight = c(0.1, 0.3, 0.34, 0.26),
#' mu = c(12L, 15L, 25L, 44L), sigma = c(4, 7, 6, 5), dfr = c(5L, 3L, 12L, 4L))
#'
#' @importFrom stats runif rt
#' @author Rawane Samb, Astrid Deschenes
#' @keywords internal
student.mixture <- function(i, k, weight, mu, sigma, dfr) {
    # Adding zero to the weight vector and calculating the cumulative sums
    sumWeight <- cumsum(c(0, weight))

    u <- runif(1, 0, 1)

    # Get the maximal position where the sum of weight is inferior to u
    position <- max(which(sumWeight < u))

    return(mu[position] + sqrt(sigma[position]) * rt(1, dfr[position]))
#     v <- c(0, weight)
#     u <- runif(1, 0, 1)
#     for (j in 1:k) {
#         if(sum(v[1:j]) < u & u <= sum(v[1:(j+1)])) {
#             mixte <- mu[j] + sqrt(sigma[j]) * rt(1, dfr[j])
#         }
#     }
#     return(mixte)
}


#' @title Normal Mixture Model
#'
#' @description Generation a value from a Normal Mixture distribution.
#'
#' @param i a \code{integer},  a count parameter.
#'
#' @param k a positive \code{integer} value, the number of nucleosomes in the
#' analyzed region.
#'
#' @param weight a \code{vector} of length \code{k}, the weight for each
#' nucleosome. The sum of all \code{weight} values must be equal to \code{1}.
#' The length of \code{weight} must be equal to \code{k}.
#'
#' @param mu a \code{vector} of positive \code{integer} of length \code{k},
#' the positions of all the nucleosomes in the analyzed region. The length
#' of \code{weight} must be equal to \code{k}.
#'
#' @param sigma a \code{vector} of length \code{k}, the variance for each
#' nucleosome. The length of \code{sigma} must be equal to \code{k}.
#'
#' @return a \code{numerical}, the value generated from a Normal Mixture
#' distribution.
#'
#' @examples
#'
#' ## Return a value generated from a normal mixture
#' rjmcmc:::normal.mixture(i = 1L, k = 4L, weight = c(0.2, 0.3, 0.24, 0.26),
#' mu = c(12L, 15L, 25L, 44L), sigma = c(4, 7, 6, 5))
#'
#' @importFrom stats runif
#' @author Rawane Samb, Astrid Deschenes
#' @keywords internal
normal.mixture <- function(i, k, weight, mu, sigma) {
    # Adding zero to the weight vector and calculating the cumulative sums
    sumWeight <- cumsum(c(0, weight))

    u <- runif(1, 0, 1)

    # Get the maximal position where the sum of weight is inferior to u
    position <- max(which(sumWeight < u))

    return(rnorm(1, mu[position], sd = sqrt(sigma[position])))

#   ANCIEN CODE MOINS RAPIDE
#   v <- c(0, weight)
#   u <- runif(1, 0, 1)
#     for (j in 1:k) {
#         if (sum(v[1:j]) < u & u <= sum(v[1:(j+1)]))
#         {
#             mixte <- rnorm(1, mu[j], sd = sqrt(sigma[j]))
#         }
#     }
#     return(mixte)
}


#' @title Prior density of \eqn{mu}
#'
#' @description Computes the prior density of \eqn{mu} conditionally to
#' the number of nucleosomes.
#'
#' For more information on the calculation of the prior density of \eqn{mu},
#' see Proposotion 1 and equation (11) of the cited article.
#'
#' @param mu a \code{vector} of positive \code{integer} containing the
#' positions of all nucleosomes.
#'
#' @param readPositions a \code{vector} of positive \code{integer}
#' corresponding to the
#' positions of all reads, including forward and reverse strands. The
#' values insinde \code{readPositions} must be sorted.
#'
#' @return  a \code{numeric}, the exact prior density of \code{mu} given the
#' number of nucleosomes.
#'
#' @references Samb R., Khadraoui K., Lakhal L., Belleau P. and Droit A. Using
#' informative Multinomial-Dirichlet prior in a t-mixture with
#' reversible jump estimation of nucleosome positions for genome-wide
#' profiling. Submitted (2015).
#' @examples
#'
#' ## Sorted vector of read positions, including forward and reverse
#' readPositions <- c(9909L, 9928L, 9935L, 26603L, 26616L, 26632L, 26636L,
#'              26640L, 44900L, 44902L, 44909L,  44910L, 44910L, 44918L,
#'              44924L, 44931L, 44935L, 44942L, 44946L)
#'
#' ## Position of the group of nucleosomes
#' mu <- c(10000L, 26700L, 45000L)
#'
#' ## Calculation of the exact prior density of mu
#' density <- rjmcmc:::priorMuDensity(mu, readPositions)
#'
#' @author Rawane Samb, Astrid Deschenes
#' @keywords internal
priorMuDensity <- function(mu, readPositions) {
    ## Get the number of nucleosomes
    k <- length(mu)
    ## Create a matrix used in the calculation of the priors
    basicMatrix <- matrix(0L, nrow = k, ncol = k)
    for (i in 1:k) {
        basicMatrix[i, i] <- 1L
    }
    if (k > 1) {
        for (i in 2:k) {
            basicMatrix[i , i - 1] <- -1L
        }
    }
    omega <- t(basicMatrix) %*% basicMatrix
    ## Calculating the range (R)
    R <- max(readPositions) - min(readPositions)
    ## Calculating the mean (E)
    E <- (max(readPositions) + min(readPositions))/2
    tau <- 1/R^2
    M <- rep(E, k)
    const <- (pi/(2*tau))^{-k/2}
    ## The calculation of the prior
    ## Equation 11 in the cited article
    return(const * exp(-(tau/2) * (t(mu - M) %*% omega %*% (mu - M))))
}


#' @title Element with the hightest number of occurences
#'
#' @description Returned the \code{integer} with the highest number of
#' occurences in a \code{vector}.
#' When more than one \code{integer} have the highest number of occurences,
#' \code{NA} is returned.
#'
#' @param sample a \code{numeric} \code{vector} (of positive \code{integer}
#' values). If the elements of \code{sample} are \code{numeric} but not
#' \code{integer}, the elements are truncated by \code{as.integer}.
#'
#' @return  a \code{integer} with the highest number of occurences or
#' \code{NA} when more than one \code{integer} have the highest number
#' of occurences.
#'
#' @author Rawane Samb, Astrid Deschenes
#' @keywords internal
#' @examples
#'
#' ## Return the element with the hightest number of occurence
#' data01 <- c(1L, 2L, 5L, 10L, 5L, 10L, 5L)
#' rjmcmc:::elementWithHighestMode(data01)
#'
#' data02 <- c(3L, 6L, 4L, 3L, 6L)
#' rjmcmc:::elementWithHighestMode(data02)
#'
elementWithHighestMode <- function(sample) {
    tabsample <- tabulate(sample)
    maxOccurence <- tabsample == max(tabsample)
    ifelse(sum(maxOccurence) == 1, which(maxOccurence), NA)
}


#' @title Merging two nucleosomal regions
#'
#' @description Merging two nucleosomal regions into one region with respect of
#' the minimal and maximal intervals allowed.
#'
#' @param yf a \code{vector} of positive \code{numeric}
#' corresponding to the
#' positions of all forward reads. The
#' values inside \code{yf} must be sorted.
#'
#' @param yr a \code{vector} of positive \code{numeric}
#' corresponding to the
#' positions of all reverse reads. The
#' values inside \code{yr} must be sorted.
#'
#' @param y a \code{vector} of positive \code{numeric}
#' corresponding to the
#' positions of all forward and reverse reads. The
#' values inside \code{y} must be sorted.
#'
#' @param list a \code{list} containing:
#' \itemize{
#'     \item k a \code{integer}, the number of nucleosomes.
#'     \item mu a \code{vector} of \code{numeric}, the positions of
#' the nucleosomes.
#'     \item sigmaf a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the forward reads for each nucleosome.
#'     \item sigmar a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the reverse reads for each nucleosome.
#'     \item delta TODO
#'     \item dl TODO
#'     \item w a \code{vector} of positive \code{numerical} of length
#' \code{k}, the weight for each nucleosome. The sum of all \code{w} values
#' must be equal to \code{1}.
#' }
#'
#' @param minInterval a \code{numeric}, the minimum distance between two
#' nucleosomes.
#'
#' @param maxInterval a \code{numeric}, the maximum distance between two
#' nucleosomes.
#'
#' @param minReads a positive \code{integer}, the minimum
#' number of reads in a potential canditate region.
#'
#' @return a \code{list} containing the updated values:
#' \itemize{
#'     \item k a \code{integer}, the number of nucleosomes.
#'     \item mu a \code{vector} of \code{numeric}, the positions of
#' the nucleosomes.
#'     \item sigmaf a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the forward reads for each nucleosome.
#'     \item sigmar a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the reverse reads for each nucleosome.
#'     \item delta TODO
#'     \item dl TODO
#'     \item w a \code{vector} of positive \code{numerical} of length
#' \code{k}, the weight for each nucleosome. The sum of all \code{w} values
#' must be equal to \code{1}.
#' }
#'
#' @author Rawane Samb, Astrid Deschenes
#' @keywords internal
mergeNucleosomes <- function(yf, yr, y, list,
                                minInterval, maxInterval, minReads) {
    ## Get the number of nucleosomes
    k <- list$k

    ## Nucleosome can only be merged when there is more than one
    if (k > 1) {
        ## Find the smallest distance between 2 nucleosomes
        ecart     <- diff(list$mu)
        ecart.min <- min(ecart)

        ## The distance between 2 nucleosomes is inferior to minimum distance
        ## allowed
        if (ecart.min < minInterval) {
            ## Merging nucleosomes with distance inferior to minimum distance
            ## up to the moment there is no nucleosome with inferior distance
            repeat
            {
                ## Find the first position with minimum gap
                p <- which.min(ecart)

                classes  <- y[y >= list$mu[p] & y < list$mu[p+1]]
                classesf <- yf[yf >= list$mu[p] & yf < list$mu[p+1]]
                classesr <- yr[yr >= list$mu[p] & yr < list$mu[p+1]]

                if (length(classes) > minReads){
                    mu <- mean(round(classes))
                } else {
                    mu <- mean(c(list$mu[p], list$mu[p+1]))
                }

                list$mu        <- sort(list$mu[-p])
                list$mu[p]     <- mu
                list$mu        <- sort(list$mu)
                list$sigmaf    <- list$sigmaf[-p]
                list$sigmar    <- list$sigmar[-p]
                list$delta     <- list$delta[-p]
                list$dl        <- list$dl[-p]
                list$w         <- list$w[-p]/sum(list$w[-p])

                # Downgrade the number of nucleosome
                k <- k - 1

                # Update the smallest distance between 2 nucleosomes
                if (k > 1) {
                    ecart     <- diff(list$mu)
                    ecart.min <- min(ecart)
                }
                if (k == 1 || ecart.min > minInterval) break()
            } ### end of boucle repeat

            ## Updating resulting list
            list <- list(  k       = k,
                            mu     = list$mu,
                            sigmaf = list$sigmaf,
                            sigmar = list$sigmar,
                            delta  = list$delta,
                            dl     = list$dl,
                            w      = list$w)
        } ### end of condition if (ecart.min < minInterval)

        ## Trying to split resulting nucleosomes
        list <- splitNucleosome(yf, yr, y, list, minInterval,
                                    maxInterval, minReads)
    } ### end of condition if (k > 1)

    return(list)
}


#' @title Spliting a nucleosomal region into two regions
#'
#' @description Split a nucleosomal region into two regions with respect of
#' the minimal and maximal intervals allowed.
#'
#' @param yf a \code{vector} of positive \code{numeric}
#' corresponding to the
#' positions of all forward reads. The
#' values inside \code{yf} must be sorted.
#'
#' @param yr a \code{vector} of positive \code{numeric}
#' corresponding to the positions of all reverse reads. The
#' values inside \code{yr} must be sorted.
#'
#' @param y a \code{vector} of positive \code{numeric}
#' corresponding to the
#' positions of all forward and reverse reads. The
#' values inside \code{y} must be sorted.
#'
#' @param list a \code{list} containing:
#' \itemize{
#'     \item k a \code{integer}, the number of nucleosomes.
#'     \item mu a \code{vector} of \code{numeric}, the positions of
#' the nucleosomes.
#'     \item sigmaf a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the forward reads for each nucleosome.
#'     \item sigmar a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the reverse reads for each nucleosome.
#'     \item delta TODO
#'     \item dl TODO
#'     \item w a \code{vector} of positive \code{numerical} of length
#' \code{k}, the weight for each nucleosome. The sum of all \code{w} values
#' must be equal to \code{1}.
#' }
#'
#' @param minInterval a \code{numeric}, the minimum distance between two
#' nucleosomes.
#'
#' @param maxInterval a \code{numeric}, the maximum distance between two
#' nucleosomes.
#'
#' @param minReads a positive \code{integer}, the minimum
#' number of reads in a potential canditate region.
#'
#' @return a \code{list} containing the updated values:
#' \itemize{
#'     \item k a \code{integer}, the number of nucleosomes.
#'     \item mu a \code{vector} of \code{numeric}, the positions of
#' the nucleosomes.
#'     \item sigmaf TODO
#'     \item sigmar a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the reverse reads for each nucleosome.
#'     \item delta a \code{vector} of \code{numeric} of length
#' \code{k}, the variance of the forward reads for each nucleosome.
#'     \item dl TODO
#'     \item w a \code{vector} of positive \code{numerical} of length
#' \code{k}, the weight for each nucleosome. The sum of all \code{w} values
#' must be equal to \code{1}.
#' }
#'
#' @author Rawane Samb, Astrid Deschenes
#' @keywords internal
splitNucleosome <- function(yf, yr, y, list, minInterval, maxInterval,
                                    minReads) {
    ## Get the number of nucleosomes
    k <- list$k

    if (k > 1) {

        ## Find the largest distance between 2 nucleosomes
        ## Its position and its value
        ecart       <- diff(list$mu)
        p           <- order(ecart, decreasing = TRUE)[1]
        ecart.max   <- ecart[p]

        if (ecart.max > maxInterval) {
            j <- 1
            repeat {
                classes <- y[y >= list$mu[p] & y < list$mu[p + 1]]
                classesf <- yf[yf >= list$mu[p] & yf < list$mu[p + 1]]
                classesr <- yr[yr >= list$mu[p] & yr < list$mu[p + 1]]

                if (length(classes) > minReads)
                {
                    new.mu <- sort(c(list$mu[1:k], mean(round(classes))))

                    new.sigmaf <- c(list$sigmaf[1:k],
                                    (list$sigmaf[p] + list$sigmaf[p + 1])/2)
                    new.sigmaf[p + 1] <- (list$sigmaf[p] + list$sigmaf[p + 1])/2
                    new.sigmaf[k + 1] <- list$sigmaf[k]

                    new.sigmar <- c(list$sigmar[1:k],
                                    (list$sigmar[p] + list$sigmar[p + 1])/2)
                    new.sigmar[p + 1] <- (list$sigmar[p] + list$sigmar[p + 1])/2
                    new.sigmar[k + 1] <- list$sigmar[k]

                    new.delta        <- c(list$delta[1:k],
                                   (list$delta[p] + list$delta[p + 1])/2)
                    new.delta[p + 1] <- (list$delta[p] + list$delta[p + 1])/2
                    new.delta[k + 1] <- list$delta[k]

                    new.dl <- round(c(list$dl[1:k],
                                        (list$dl[p] + list$dl[p + 1])/2))
                    new.dl[p + 1]   <- (list$dl[p] + list$dl[p + 1])/2
                    new.dl[k + 1]   <- list$dl[k]

                    new.w           <- c(list$w[1:k], (list$w[p] +
                                                            list$w[p + 1])/2)
                    new.w[p + 1]    <- (list$w[p]+list$w[p + 1])/2
                    new.w[k + 1]    <- list$w[k]

                    k       <- length(new.mu)
                    list    <- list(k       = k,
                                    mu      = new.mu,
                                    sigmaf  = new.sigmaf,
                                    sigmar  = new.sigmar,
                                    delta   = new.delta,
                                    dl      = new.dl,
                                    w       = new.w/sum(new.w))

                    ## Update the vector of distance between nucleosomes
                    ecart <- diff(list$mu)
                }
                else
                {
                    ## Update to select the next maximum distance value
                    j  <- j + 1
                }

                ## Select the next nucleosome to be potentially split
                p           <- order(ecart, decreasing = TRUE)[j]
                ecart.max   <- ecart[p]

                if ( j >= (k - 1) || ecart.max <= maxInterval) break()
            } ### end of boucle repeat
        } ### end of condition if (ecart.max > maxInterval)
    } ### end of condition if (k > 1)
    return(list)
}


#' @title Parameters validation for the \code{\link{RJMCMC}}
#' function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{RJMCMC}} function.
#'
#' @param startPosForwardReads a \code{vector} of positive \code{integer}, the
#' start position of all the forward reads.
#'
#' @param startPosReverseReads a \code{vector} of positive \code{integer}, the
#' positions of all the reverse reads. Beware that the start position of
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
#' @param verbose a \code{logical} indicating if extra information must be
#' printed or not.
#'
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @author Astrid Deschenes
#' @keywords internal
validateParameters <- function(startPosForwardReads, startPosReverseReads,
                                    nbrIterations, kmax, lambda,
                                    minInterval, maxInterval, minReads,
                                    verbose) {
    ## Validate the nbrIterations parameter
    if (!isInteger(nbrIterations) || as.integer(nbrIterations) < 1) {
        stop("nbrIterations must be a positive integer or numeric")
    }

    ## Validate the kmax parameter
    if (!isInteger(kmax) || as.integer(kmax) < 1) {
        stop("kmax must be a positive integer or numeric")
    }

    ## Validate the minReads parameter
    if (!isInteger(minReads) || as.integer(minReads) < 1) {
        stop("minReads must be a positive integer or numeric")
    }

    ## Validate the lambda parameter
    if (!isInteger(lambda) || lambda <= 0) {
        stop("lambda must be a positive numeric")
    }

    ## Validate that the startPosForwardReads has at least one read
    ## and that the values are integer
    if (!is.vector(startPosForwardReads) || !is.numeric(startPosForwardReads)
        || length(startPosForwardReads) < 1)
    {
        stop(paste0("startPosForwardReads must be a non-empty vector of ",
                    "numeric values."))
    }

    ## Validate that the startPosReverseReads has at least one read
    ## and that the values are integer
    if (!is.vector(startPosReverseReads) || !is.numeric(startPosReverseReads)
        || length(startPosReverseReads) < 1)
    {
        stop(paste0("startPosReverseReads must be a non-empty vector of ",
                    "numeric values."))
    }

    ## Validate that verbose is a logical
    if (!is.logical(verbose)) {
        stop("verbose must be a logical.")
    }

    return(0)
}


#' @title Validate if a value is an integer
#'
#' @description Validate if the value passed to the function can be casted
#' into a \code{integer} or
#' not. The value must have a length
#' of 1. The type of value can be a \code{integer} or \code{numerical}.
#'
#' @param value an object to validate.
#'
#' @return \code{TRUE} is the parameter is a integer; otherwise \code{FALSE}
#'
#' @examples
#'
#' ## Return TRUE because the input is an integer of length 1
#' rjmcmc:::isInteger(33L)
#'
#' ## Return FALSE because the length of the input is not 1
#' rjmcmc:::isInteger(c(21L, 33L))
#'
#' ## Return TRUE because the input is a numericalof length 1
#' rjmcmc:::isInteger(323.21)
#'
#' ## Return FALSE because the length of the input is not 1
#' rjmcmc:::isInteger(c(444.2, 442.1))
#'
#' @author Astrid Deschenes
#' @keywords internal
#'
isInteger <- function(value) {
    return((is.integer(value) && length(value) == 1) || (is.numeric(value) &&
        length(value) == 1))
}



#' @title Birth move in the case that only one nucleosome is present
#'
#' @description  Attempt to add a new nucleosome in the case that only one
#' nucleosome is present, \code{case k = 1}.
#'
#' @param paramValues a \code{list} containing:
#' \itemize{
#'     \item startPSF a \code{vector} of \code{numeric}, the
#' start position of all the forward reads.
#'     \item startPSR a \code{vector} of \code{numeric}, the
#' start position of all the reverse reads.
#'     \item kmax a positive \code{integer} or \code{numeric}, the maximum
#' number of nucleosomes per region.
#'     \item lambda a positive \code{numeric}, the theorical mean
#' of the Poisson distribution.
#'     \item minReads a positive \code{integer} or \code{numeric}, the minimum
#' number of reads in a potential canditate region.
#'     \item y TODO
#'     \item nr TODO
#'     \item nf TODO
#'     \item nbrReads TODO
#'     \item zeta TODO
#'     \item deltamin TODO
#'     \item detlamax TODO
#'     \item minReadPos TODO
#'     \item maxReadPos TODO
#' }
#'
#' @param kValue a \code{integer}, the number of nucleosomes.
#'
#' @param muValue a \code{vector} of \code{numeric} of length
#' \code{kValue}, the positions of the nucleosomes.
#'
#' @param sigmafValue a \code{vector} of \code{numeric} of length
#' \code{kValue}, the variance of the forward reads for each nucleosome.
#'
#' @param sigmarValue a \code{vector} of \code{numeric} of length
#' \code{kValue}, the variance of the reverse reads for each nucleosome.
#'
#' @param deltaValue
#'
#' @param wValue a \code{vector} of positive \code{numerical} of length
#' \code{kValue}, the weight for each nucleosome. The sum of all \code{wValue}
#' values must be equal to \code{1}.
#'
#' @param dlValue
#'
#' @param aValue
#'
#' @param dimValue
#'
#' @param paramValues = list(startPSF =  startPosForwardReads, startPSR = startPosReverseReads
#'    , kmax = kmax, lambda = lambda, minReads = minReads
#'    , y = y, nr = nr, nf =nf, nbrReads = nbrReads , zeta = zeta
#'    , deltamin = deltamin, deltamax = deltamax
#'    , minReadPos= minReadPos, maxReadPos = maxReadPos )
#'    kValue, muValue, sigmafValue, sigmarValue, deltaValue, wValue, dlValue, aValue, dimValue
#'
#' @return list( k=0L, mu=numeric(parValue$kmax)
#'    , sigmaf=numeric(parValue$kmax), sigmar=numeric(parValue$kmax)
#'    , delta=numeric(parValue$kmax), w=numeric(parValue$kmax)
#'    , dl=numeric(parValue$kmax), a=numeric(parValue$kmax)
#'    , dim=numeric(parValue$kmax), rho=0)
#'
#' @examples
#' birthMoveK1(paramValues , kValue, muValue, sigmafValue, sigmarValue,
#' deltaValue, wValue, dlValue, aValue, dimValue )
#'
#'
#' @author Rawane Samb, Pascal Belleau, Astrid Deschenes
#'
birthMoveK1 <- function(paramValues, kValue, muValue, sigmafValue,
                            sigmarValue, deltaValue, wValue, dlValue,
                            aValue, dimValue) {

    varTilde <- list( k=0L, mu=numeric(parValue$kmax)
                      , sigmaf=numeric(parValue$kmax), sigmar=numeric(parValue$kmax)
                      , delta=numeric(parValue$kmax), w=numeric(parValue$kmax)
                      , dl=numeric(parValue$kmax), a=numeric(parValue$kmax)
                      , dim=numeric(parValue$kmax), rho=0)

    Kn              <- rep(0, nbrIterations)
    Kaf             <- matrix(0, nrow = paramValues$nf, ncol = kmax)
    Kbf             <- matrix(0, nrow = paramValues$nf, ncol = kmax)
    Kar             <- matrix(0, nrow = paramValues$nr, ncol = kmax)
    Kbr             <- matrix(0, nrow = paramValues$nr, ncol = kmax)

    Y1f             <- rep(0, paramValues$nf)
    Y2f             <- rep(0, paramValues$nf)
    Y1r             <- rep(0, paramValues$nr)
    Y2r             <- rep(0, paramValues$nr)


    varTilde$k <- kValue + 1L
    count  <- 1L
    repeat {
        j <- sample(1:kValue, 1)
        varTilde$mu[j] <- runif(1, paramValues$minReadPos, muValue[ j])
        varTilde$mu[ 1:varTilde$k] <- sort(c(muValue[ 1:kValue],
                                             varTilde$mu[j]))

        varTilde$a[ j+1] <- runif(1, varTilde$mu[j], varTilde$mu[j+1])
        varTilde$a[ 1:(varTilde$k+1)] <- sort(c(aValue[ 1:varTilde$k],
                                                varTilde$a[ j+1]))
        varTilde$a[ 1]                <- paramValues$minReadPos
        varTilde$a[ (varTilde$k+1)]    <- paramValues$maxReadPos

        varTilde$dim[ 1] <- length(paramValues$y[varTilde$a[ 1] <= paramValues$y &
                                                     paramValues$y < varTilde$a[ 2]])
        varTilde$dim[ varTilde$k] <- length(paramValues$y[varTilde$a[ varTilde$k] <=
                                                              paramValues$y & paramValues$y <= paramValues$maxReadPos])
        if (varTilde$k > 2) {   # impossible si kValue == 1
            for (m in 2:(varTilde$k-1)) {
                varTilde$dim[m] <- length(paramValues$y[(varTilde$a[ m] <= paramValues$y &
                                                             paramValues$y < varTilde$a[ m+1])])
            }
        }
        Pr <- min(varTilde$dim[ 1:varTilde$k])

        ybar <- mean(paramValues$y[varTilde$a[ j] <= paramValues$y & paramValues$y <= varTilde$a[ j+1]])
        classesf <- paramValues$y[varTilde$a[ j] <= paramValues$y & paramValues$y <= ybar]
        classesr <- paramValues$y[ybar <= paramValues$y & paramValues$y <= varTilde$a[ j+1]]

        Lf <- length(classesf[!duplicated(classesf)])
        Lr <- length(classesr[!duplicated(classesr)])

        count <- count + 1L

        if ((Pr > 1 & Lf > 1 & Lr > 1)  ||
                count == 1000L) break()
    }

    if (count == 1000L) {
        varTilde$rho <- 0
    } else  {

        varTilde$dl[ j] <- sample(3:30, 1)

        varTilde$sigmaf[ j] <- ifelse(Lf > 1,
                                      var(classesf) * (varTilde$dl[ j] - 2)/varTilde$dl[ j],
                                      sigmafValue[ j])
        varTilde$sigmar[ j] <- ifelse(Lr > 1,
                                      var(classesr) * (varTilde$dl[ j] - 2)/varTilde$dl[ j],
                                      sigmarValue[ j])

        varTilde$sigmaf[ 1:varTilde$k] <- c(sigmafValue[ 1:kValue],
                                            varTilde$sigmaf[ j])
        varTilde$sigmar[ 1:varTilde$k] <- c(sigmarValue[ 1:kValue],
                                            varTilde$sigmar[ j])

        varTilde$delta[ j] <- tnormale(paramValues$zeta,
                                       1/(varTilde$sigmaf[ j]^{-1} + varTilde$sigmar[j]^{-1}),
                                       paramValues$deltamin, paramValues$deltamax)
        varTilde$delta[ 1:varTilde$k] <- c(varTilde$delta[ j],
                                           deltaValue[ 1:kValue])

        alpha                   <- rep(1, kValue)
        alphatilde              <- rep(1, varTilde$k)
        alphaproptilde          <- rep(1, varTilde$k)
        alphaprop               <- rep(1, kValue)
        varTilde$w[ 1:varTilde$k]  <- rdirichlet(1, alphaproptilde)
        ennetilde               <- varTilde$dim[ 1:varTilde$k]
        enne <- rmultinom(1, paramValues$nbrReads, wValue[ 1:kValue])

        #Rapport de vraisemblance

        for (m in 1:kValue) {
            Kaf[,m] <- (varTilde$w[ m]*(1/sqrt(varTilde$sigmaf[m]))*dt((paramValues$startPSF-varTilde$mu[m]+varTilde$delta[m]/2)/sqrt(varTilde$sigmaf[m]),varTilde$dl[m]))
            Kbf[,m] <- (wValue[ m]*(1/sqrt(sigmafValue[m]))*dt((paramValues$startPSF-muValue[m]+deltaValue[m]/2)/sqrt(sigmafValue[m]),dlValue[m]))
        }
        Kaf[,varTilde$k] <- (varTilde$w[ varTilde$k]*(1/sqrt(varTilde$sigmaf[varTilde$k]))*dt(( paramValues$startPSF -varTilde$mu[varTilde$k]+varTilde$delta[m]/2)/sqrt(varTilde$sigmaf[varTilde$k]),varTilde$dl[m]))
        for (s in 1:paramValues$nf) {
            Y1f[s] <- log(sum(Kaf[s, 1:varTilde$k]))
            Y2f[s] <- log(sum(Kbf[s, 1:kValue]))
        }

        for (m in 1:kValue) {
            Kar[,m] <- (varTilde$w[ m]*(1/sqrt(varTilde$sigmar[m]))*dt((paramValues$startPSR-varTilde$mu[m]-varTilde$delta[m]/2)/sqrt(varTilde$sigmar[m]),varTilde$dl[m]))
            Kbr[,m] <- (wValue[ m]*(1/sqrt(sigmarValue[m]))*dt((paramValues$startPSR-muValue[m]-deltaValue[m]/2)/sqrt(sigmarValue[m]),dlValue[m]))
        }
        Kar[,varTilde$k] <- (varTilde$w[ varTilde$k]*(1/sqrt(varTilde$sigmar[varTilde$k]))*dt(( paramValues$startPSR -varTilde$mu[varTilde$k]-varTilde$delta[m]/2)/sqrt(varTilde$sigmar[varTilde$k]),varTilde$dl[m]))
        for (s in 1:paramValues$nr) {
            Y1r[s] <- log(sum(Kar[s, 1:varTilde$k]))
            Y2r[s] <- log(sum(Kbr[s, 1:kValue]))
        }

        Kn <- sum(Y1f) + sum(Y1r)
        Kd <- sum(Y2f) + sum(Y2r)

        q <- Kn - Kd
        rap.q <- exp(q)

        rap.vrais <- rap.q

        if (j == 1) {
            qalloc <- 1/(muValue[ j] - paramValues$minReadPos) # Density of varTilde$mu[j]
        } else {
            qalloc <- 1/(muValue[ j] - muValue[ j-1])
        }

        rap.priormu   <- (priorMuDensity(varTilde$mu[ 1:varTilde$k],paramValues$y)/priorMuDensity(muValue[ 1:kValue],paramValues$y))
        rap.priorw    <- (ddirichlet(varTilde$w[ 1:varTilde$k],alphatilde)/ddirichlet(wValue[ 1:kValue],alpha) )
        rap.priorenne <- dmultinom(ennetilde, paramValues$nbrReads,varTilde$w[1:varTilde$k])/dmultinom(dimValue[1:kValue], paramValues$nbrReads,wValue[ 1:kValue])
        rap.priork    <- (dpois(varTilde$k, paramValues$lambda)/dpois(kValue, paramValues$lambda))
        rap.propmu    <- (1/(qalloc))
        rap.propw     <- (ddirichlet(wValue[ 1:kValue], alphaprop)/ddirichlet(varTilde$w[ 1:varTilde$k], alphaproptilde))

        rap.prior     <- rap.priormu * rap.priorw * rap.priorenne * rap.priork
        rap.prop      <-  rap.propmu  * rap.propw

        varTilde$rho       <- min(1,(rap.vrais) * (rap.prior) * (rap.prop ) * (Dk(varTilde$k,paramValues$lambda,paramValues$kmax)/Bk(kValue, paramValues$lambda, paramValues$kmax)))

    }

    varTilde$rho <- ifelse (is.na(varTilde$rho), 0, varTilde$rho)
    return(varTilde)
}

birthMove <- function(paramValues, kValue, muValue, sigmafValue, sigmarValue, deltaValue, wValue, dlValue, aValue, dimValue ){
    #Birth move

    varTilde <- list( k=0L, mu=numeric(parValue$kmax)
                      , sigmaf=numeric(parValue$kmax), sigmar=numeric(parValue$kmax)
                      , delta=numeric(parValue$kmax), w=numeric(parValue$kmax)
                      , dl=numeric(parValue$kmax), a=numeric(parValue$kmax)
                      , dim=numeric(parValue$kmax), rho=0)

    Kn              <- rep(0, nbrIterations)
    Kaf             <- matrix(0, nrow = paramValues$nf, ncol = kmax)
    Kbf             <- matrix(0, nrow = paramValues$nf, ncol = kmax)
    Kar             <- matrix(0, nrow = paramValues$nr, ncol = kmax)
    Kbr             <- matrix(0, nrow = paramValues$nr, ncol = kmax)

    Y1f             <- rep(0, paramValues$nf)
    Y2f             <- rep(0, paramValues$nf)
    Y1r             <- rep(0, paramValues$nr)
    Y2r             <- rep(0, paramValues$nr)


    varTilde$k <- kValue + 1L
    count <- 1L

    repeat {
        j <- sample(1:kValue, 1)

        if (j == 1) {
            varTilde$mu[j] <- runif(1, paramValues$minReadPos, muValue[j])
        } else {
            varTilde$mu[j] <- runif(1, muValue[j-1], muValue[j])
        }

        varTilde$mu[ 1:varTilde$k] <- sort(c(muValue[ 1:kValue],varTilde$mu[ j]))

        varTilde$a[j+1] <- ifelse(j<ktilde[i-1], runif(1,varTilde$mu[j],varTilde$mu[j+1]), runif(1,varTilde$mu[j],paramValues$maxReadPos))
        varTilde$a[1:(varTilde$k+1)] <- sort(c(aValue[1:varTilde$k],varTilde$a[j+1]))

        if (j == 1) {
            varTilde$a[j] <- paramValues$minReadPos
        } else {
            varTilde$a[j] <- runif(1,varTilde$mu[j-1],varTilde$mu[j])
        }
        varTilde$a[1]             <- paramValues$minReadPos
        varTilde$a[varTilde$k+1]   <- paramValues$maxReadPos

        varTilde$dim[1]         <- length(paramValues$y[varTilde$a[1]<=paramValues$y & paramValues$y<varTilde$a[2]])
        varTilde$dim[varTilde$k] <- length(paramValues$y[varTilde$a[varTilde$k]<=paramValues$y & paramValues$y<=paramValues$maxReadPos])

        if (varTilde$k > 2) {
            for (m in 2: (varTilde$k-1)) {
                varTilde$dim[m] <- length(paramValues$y[varTilde$a[m]<=paramValues$y & paramValues$y<varTilde$a[m+1]])
            }
        }
        Pr <- min(varTilde$dim[1:varTilde$k])

        ybar <- mean(paramValues$y[varTilde$a[j] <= paramValues$y & paramValues$y <= varTilde$a[j+1]])
        classesf <- paramValues$y[varTilde$a[j] <= paramValues$y & paramValues$y <= ybar]
        classesr <- paramValues$y[ybar <= paramValues$y & paramValues$y <= varTilde$a[j+1]]

        Lf <- length(classesf[!duplicated(classesf)])
        Lr <- length(classesr[!duplicated(classesr)])
        count <- count + 1L
        if ( (Pr>1 & Lf>1 & Lr>1)  ||
                 count == 1000L) break()

    }

    if (count == 1000L) {
        varTilde$rho <- 0
    } else {
        varTilde$dl[j] <- sample(3:30, 1)
        varTilde$sigmaf[j] <- ifelse(Lf > 1, var(classesf)*(varTilde$dl[j]-2)/varTilde$dl[j], sigmafValue[j])
        varTilde$sigmar[j] <- ifelse(Lr > 1, var(classesr)*(varTilde$dl[j]-2)/varTilde$dl[j], sigmarValue[j])

        if (j == 1) {
            varTilde$sigmaf[1:varTilde$k]<- c(sigmafValue[1:kValue],varTilde$sigmaf[j])
            varTilde$sigmar[1:varTilde$k]<- c(sigmarValue[1:kValue],varTilde$sigmar[j] )
        } else {
            varTilde$sigmaf[1:varTilde$k] <- c(sigmafValue[1:(j-1)],varTilde$sigmaf[j],sigmafValue[j:kValue])
            varTilde$sigmar[1:varTilde$k] <- c(sigmarValue[1:(j-1)],varTilde$sigmar[j],sigmarValue[j:kValue])
        }

        varTilde$delta[j] <- tnormale(paramValues$zeta, 1/(varTilde$sigmaf[j]^{-1}+varTilde$sigmar[j]^{-1}), paramValues$deltamin, paramValues$deltamax)

        if (j == 1) {
            varTilde$delta[ 1:varTilde$k] <- c(varTilde$delta[ j], deltaValue[ 1:kValue])
        } else if (j == kValue) {
            varTilde$delta[ 1:varTilde$k] <- c(deltaValue[ 1:kValue], varTilde$delta[ j])
        } else {
            varTilde$delta[ 1:varTilde$k] <- c(deltaValue[ 1:(j-1)], varTilde$delta[j], deltaValue[j:kValue])
        }

        alpha                 <- rep(1, kValue)
        alphatilde            <- rep(1, varTilde$k)
        alphaproptilde        <- rep(1, varTilde$k)
        alphaprop             <- rep(1, kValue)
        ennetilde             <- varTilde$dim[ 1:varTilde$k]
        enne                  <- rmultinom(1, paramValues$nbrReads, wValue[ 1:kValue])
        varTilde$w[1:varTilde$k] <- rdirichlet(1, alphaproptilde)

        ### Rapport de vraisemblance ###

        for (m in 1:kValue) {
            Kaf[,m]<- (varTilde$w[m]*(1/sqrt(varTilde$sigmaf[m]))*dt(( paramValues$startPSF -varTilde$mu[m]+varTilde$delta[m]/2)/sqrt(varTilde$sigmaf[m]),varTilde$dl[m]))
            Kbf[,m]<- (wValue[m]*(1/sqrt(sigmafValue[m]))*dt(( paramValues$startPSF -muValue[m]+deltaValue[m]/2)/sqrt(sigmafValue[m]),dlValue[m]))
        }
        Kaf[,varTilde$k]<- (varTilde$w[varTilde$k]*(1/sqrt(varTilde$sigmaf[varTilde$k]))*dt(( paramValues$startPSF -varTilde$mu[varTilde$k]+varTilde$delta[m]/2)/sqrt(varTilde$sigmaf[varTilde$k]),varTilde$dl[m]))
        for (s in 1:paramValues$nf) {
            Y1f[s] <-log(sum(Kaf[s,1:varTilde$k]))
            Y2f[s] <-log(sum(Kbf[s,1:kValue]))
        }

        for (m in 1:kValue) {
            Kar[,m]<- (varTilde$w[m]*(1/sqrt(varTilde$sigmar[m]))*dt(( paramValues$startPSR -varTilde$mu[m]-varTilde$delta[m]/2)/sqrt(varTilde$sigmar[m]),varTilde$dl[m]))
            Kbr[,m]<- (wValue[m]*(1/sqrt(sigmarValue[m]))*dt(( paramValues$startPSR -muValue[m]-deltaValue[m]/2)/sqrt(sigmarValue[m]),dlValue[m]))
        }
        Kar[,varTilde$k]<- (varTilde$w[varTilde$k]*(1/sqrt(varTilde$sigmar[varTilde$k]))*dt(( paramValues$startPSR -varTilde$mu[varTilde$k]-varTilde$delta[m]/2)/sqrt(varTilde$sigmar[varTilde$k]),varTilde$dl[m]))
        for (s in 1:paramValues$nr) {
            Y1r[s] <-log(sum(Kar[s,1:varTilde$k]))
            Y2r[s] <-log(sum(Kbr[s,1:kValue]))
        }

        Kn <- sum(Y1f)+sum(Y1r)
        Kd <- sum(Y2f)+sum(Y2r)

        q         <- Kn - Kd
        rap.q     <-exp(q)

        rap.vrais <- rap.q

        #Density of
        varTilde$mu[j]
        if (j==1) {
            qalloc <- 1/(muValue[ j] - paramValues$minReadPos)
        } else {
            qalloc <- 1/(muValue[ j] - muValue[ j-1])
        }

        rap.priormu   <- (priorMuDensity(varTilde$mu[1:varTilde$k],paramValues$y)/priorMuDensity(muValue[ 1:kValue],paramValues$y))
        rap.priorw    <- (ddirichlet(varTilde$w[1:varTilde$k], alphatilde)/ddirichlet(wValue[ 1:kValue], alpha))
        rap.priorenne <- dmultinom(ennetilde, paramValues$nbrReads,varTilde$w[1:varTilde$k])/dmultinom(dimValue[ 1:kValue], paramValues$nbrReads,wValue[ 1:kValue])
        rap.priork    <- (dpois(varTilde$k,paramValues$lambda)/dpois(kValue,paramValues$lambda))
        rap.propmu    <- (1/(qalloc))
        rap.propw     <- (ddirichlet(wValue[1:kValue],alphaprop)/ddirichlet(varTilde$w[1:varTilde$k],alphaproptilde))

        rap.prior     <- rap.priormu * rap.priorw * rap.priorenne * rap.priork
        rap.prop      <- rap.propmu  * rap.propw
        varTilde$rho       <- min(1, (rap.vrais) * (rap.prior)  *  (rap.prop ) * (Dk(varTilde$k,paramValues$lambda,paramValues$kmax)/Bk(kValue,paramValues$lambda,paramValues$kmax)))
    }



    varTilde$rho <- ifelse(is.na(varTilde$rho), 0, varTilde$rho)

    return(varTilde)
}


mhMoveK1 <- function(paramValues , kValue, muValue, sigmafValue, sigmarValue, deltaValue, wValue, dlValue, aValue, dimValue ){
    ###Metropolis-Hastings move
    varTilde <- list( k=0L, mu=numeric(parValue$kmax)
                      , sigmaf=numeric(parValue$kmax), sigmar=numeric(parValue$kmax)
                      , delta=numeric(parValue$kmax), w=numeric(parValue$kmax)
                      , dl=numeric(parValue$kmax), a=numeric(parValue$kmax)
                      , dim=numeric(parValue$kmax), rho=0)

    Kn              <- rep(0, nbrIterations)
    Kaf             <- matrix(0, nrow = paramValues$nf, ncol = kmax)
    Kbf             <- matrix(0, nrow = paramValues$nf, ncol = kmax)
    Kar             <- matrix(0, nrow = paramValues$nr, ncol = kmax)
    Kbr             <- matrix(0, nrow = paramValues$nr, ncol = kmax)

    Y1f             <- rep(0, paramValues$nf)
    Y2f             <- rep(0, paramValues$nf)
    Y1r             <- rep(0, paramValues$nr)
    Y2r             <- rep(0, paramValues$nr)

    varTilde$k <- kValue
    count     <- 1L
    repeat {
        j                        <- sample(1:kValue, 1)
        varTilde$mu[ j]            <- runif(1, muValue[ j], paramValues$maxReadPos)
        varTilde$mu[ 1:varTilde$k]  <- sort(c(varTilde$mu[ 1:varTilde$k]))

        varTilde$a[ j]   <- paramValues$minReadPos
        varTilde$a[ j+1] <- paramValues$maxReadPos

        varTilde$dim[ 1]           <- length(paramValues$y[varTilde$a[1] <=paramValues$y &  paramValues$y< varTilde$a[2]])
        varTilde$dim[ varTilde$k]   <- length(paramValues$y[varTilde$a[varTilde$k] <= paramValues$y & paramValues$y <= paramValues$maxReadPos])

        if (varTilde$k > 2) {
            for (m in 2:(varTilde$k-1)) {
                varTilde$dim[m] <-length(paramValues$y[(varTilde$a[ m] <=paramValues$y & paramValues$y < varTilde$a[ m + 1])])
            }
        }

        Pr <- min(varTilde$dim[ 1:varTilde$k])

        ybar        <- mean(paramValues$y[varTilde$a[j]<=paramValues$y & paramValues$y<=varTilde$a[j+1]])
        classesf    <- paramValues$y[varTilde$a[j]<=paramValues$y & paramValues$y<=ybar]
        classesr    <- paramValues$y[ybar <= paramValues$y & paramValues$y <= varTilde$a[ j + 1]]

        Lf      <- length(classesf[!duplicated(classesf)])
        Lr      <- length(classesr[!duplicated(classesr)])
        count   <- count + 1L

        if ((Pr > 1 & Lf > 1 & Lr > 1)  ||
                count == 1000L) break()
    }

    if (count == 1000L) {
        varTilde$rho <- 0
    } else {

        varTilde$sigmaf[ 1:varTilde$k] <- sigmafValue[ 1:kValue]
        varTilde$sigmar[ 1:varTilde$k] <- sigmarValue[ 1:kValue]

        varTilde$dl[ j] <- sample(3:30,1)

        varTilde$sigmaf[ j] <- ifelse(Lf>1, var(classesf)*(varTilde$dl[j]-2)/varTilde$dl[j], sigmafValue[j])
        varTilde$sigmar[ j] <- ifelse(Lr>1, var(classesr)*(varTilde$dl[j]-2)/varTilde$dl[j], sigmarValue[j])

        varTilde$delta[ 1:varTilde$k] <- deltaValue[1:kValue]

        varTilde$delta[ j] <- tnormale(paramValues$zeta, 1/(varTilde$sigmaf[j]^{-1}+varTilde$sigmar[j]^{-1}), paramValues$deltamin, paramValues$deltamax)

        alpha <- rep(1, kValue)
        alphatilde <- rep(1,varTilde$k)
        ennetilde <- varTilde$dim[1:varTilde$k]
        alphaproptilde <- rep(1,varTilde$k)
        alphaprop <- rep(1, kValue)
        varTilde$w[ 1:varTilde$k] <- rdirichlet(1,alphaproptilde)

        ### calcul du rapport de vraisemblance de M-H move

        for (m in 1:varTilde$k) {
            Kaf[,m] <- (varTilde$w[ m]*(1/sqrt(varTilde$sigmaf[m]))*dt(( paramValues$startPSF -varTilde$mu[m]+varTilde$delta[m]/2)/sqrt(varTilde$sigmaf[m]),varTilde$dl[m]))
            Kbf[,m] <- (wValue[ m]*(1/sqrt(sigmafValue[m]))*dt(( paramValues$startPSF -muValue[m]+deltaValue[m]/2)/sqrt(sigmafValue[m]),dlValue[m]))
        }
        for (s in 1:paramValues$nf) {
            Y1f[s] <-log(sum(Kaf[s, 1:varTilde$k]))
            Y2f[s] <-log(sum(Kbf[s, 1:kValue]))
        }

        for (m in 1:varTilde$k) {
            Kar[,m] <- (varTilde$w[ m]*(1/sqrt(varTilde$sigmar[m]))*dt((paramValues$startPSR -varTilde$mu[m]-varTilde$delta[m]/2)/sqrt(varTilde$sigmar[m]),varTilde$dl[m]))
            Kbr[,m] <- (wValue[ m]*(1/sqrt(sigmarValue[m]))*dt((paramValues$startPSR -muValue[m]-deltaValue[m]/2)/sqrt(sigmarValue[m]),dlValue[m]))
        }
        for (s in 1:paramValues$nr) {
            Y1r[s] <- log(sum(Kar[s, 1:varTilde$k]))
            Y2r[s] <- log(sum(Kbr[s, 1:kValue]))
        }

        Kn <- sum(Y1f) + sum(Y1r)
        Kd <- sum(Y2f) + sum(Y2r)

        q <- Kn - Kd
        rap.q <- exp(q)

        rap.vrais <- rap.q

        rap.priormu   <- (priorMuDensity(varTilde$mu[ 1:varTilde$k], paramValues$y)/priorMuDensity(muValue[ 1:kValue], paramValues$y))
        rap.priorw    <- (ddirichlet(varTilde$w[1:varTilde$k],alphatilde)/ddirichlet(wValue[ 1:kValue],alpha) )
        rap.priorenne <- dmultinom(ennetilde, paramValues$nbrReads,varTilde$w[1:varTilde$k])/dmultinom(dimValue[ 1:kValue], paramValues$nbrReads,wValue[1:kValue])
        rap.priork    <- 1
        rap.propmu    <- 1
        rap.propw     <- (ddirichlet(wValue[ 1:kValue],alphaprop)/ddirichlet(varTilde$w[ 1:varTilde$k], alphaproptilde))

        rap.prior     <- rap.priormu * rap.priorw * rap.priorenne * rap.priork
        rap.prop      <-  rap.propmu  * rap.propw
        varTilde$rho      <- min(1, rap.vrais * (rap.prior)  *  (rap.prop))

    }

    varTilde$rho <- ifelse(is.na(varTilde$rho), 0, varTilde$rho)
    return(varTilde)
}

mhMove <- function(paramValues , kValue, muValue, sigmafValue, sigmarValue, deltaValue, wValue, dlValue, aValue, dimValue ){
    ### Metropolis-Hastings move
    varTilde <- list( k=0L, mu=numeric(parValue$kmax)
                      , sigmaf=numeric(parValue$kmax), sigmar=numeric(parValue$kmax)
                      , delta=numeric(parValue$kmax), w=numeric(parValue$kmax)
                      , dl=numeric(parValue$kmax), a=numeric(parValue$kmax)
                      , dim=numeric(parValue$kmax), rho=0)

    Kn              <- rep(0, nbrIterations)
    Kaf             <- matrix(0, nrow = paramValues$nf, ncol = kmax)
    Kbf             <- matrix(0, nrow = paramValues$nf, ncol = kmax)
    Kar             <- matrix(0, nrow = paramValues$nr, ncol = kmax)
    Kbr             <- matrix(0, nrow = paramValues$nr, ncol = kmax)

    Y1f             <- rep(0, paramValues$nf)
    Y2f             <- rep(0, paramValues$nf)
    Y1r             <- rep(0, paramValues$nr)
    Y2r             <- rep(0, paramValues$nr)

    varTilde$k   <- kValue
    count       <- 1L
    repeat {
        j                       <- sample(2:kValue, 1)
        varTilde$mu[1:varTilde$k]  <- muValue[ 1:kValue]

        if (j==1) {
            varTilde$mu[j] <- runif(1, paramValues$minReadPos, muValue[j+1])
        } else {
            if (j==varTilde$k) {
                varTilde$mu[j] <- runif(1,muValue[j-1],paramValues$maxReadPos)
            } else {
                varTilde$mu[j] <- runif(1,muValue[j-1],muValue[j+1])
            }
        }

        varTilde$mu[1:varTilde$k] <- sort(c(varTilde$mu[ 1:varTilde$k]))

        varTilde$a[1:(varTilde$k+1)] <- sort(c(aValue[ 1:(kValue + 1)]))

        if (j==varTilde$k) {
            varTilde$a[j] <- runif(1,varTilde$mu[j],paramValues$maxReadPos)
            varTilde$a[j+1] <- paramValues$maxReadPos
        } else {
            if (j==1) {
                varTilde$a[j] <- paramValues$minReadPos
                varTilde$a[j+1] <- runif(1,varTilde$mu[j],varTilde$mu[j+1])
            } else {
                varTilde$a[j] <- runif(1,varTilde$mu[j-1],varTilde$mu[j])
                varTilde$a[j+1] <- runif(1,varTilde$mu[j],varTilde$mu[j+1])
            }
        }

        varTilde$a[1]           <- paramValues$minReadPos
        varTilde$a[varTilde$k+1] <- paramValues$maxReadPos

        varTilde$dim[1]         <- length(paramValues$y[varTilde$a[1]<=paramValues$y & paramValues$y<varTilde$a[2]])
        varTilde$dim[varTilde$k] <- length(paramValues$y[varTilde$a[varTilde$k]<=paramValues$y & paramValues$y<=paramValues$maxReadPos])
        if (varTilde$k>2) {
            for (m in 2: (varTilde$k-1)) {
                varTilde$dim[m]<-length(paramValues$y[varTilde$a[m]<=paramValues$y & paramValues$y<varTilde$a[m+1]])
            }
        }
        Pr <- min(varTilde$dim[1:varTilde$k])

        ybar      <- mean(paramValues$y[varTilde$a[j]<=paramValues$y & paramValues$y<=varTilde$a[j+1]])
        classesf  <- paramValues$y[varTilde$a[j]<=paramValues$y & paramValues$y<=ybar]
        classesr  <- paramValues$y[ybar<=paramValues$y & paramValues$y<=varTilde$a[j+1]]

        Lf    <- length(classesf[!duplicated(classesf)])
        Lr    <- length(classesr[!duplicated(classesr)])
        count <- count + 1L
        if ( (Pr > 1 & Lf > 1 & Lr > 1)
             || count == 1000L) break()
    }

    if (count == 1000L) {
        varTilde$rho <- 0
    } else {

        varTilde$sigmaf[ 1:varTilde$k] <- sigmafValue[ 1:kValue]
        varTilde$sigmar[ 1:varTilde$k] <- sigmarValue[ 1:kValue]

        varTilde$dl[j] <- sample(3:30, 1)

        varTilde$sigmaf[j] <- ifelse(Lf>1, var(classesf)*(varTilde$dl[j]-2)/varTilde$dl[j], sigmafValue[j])
        varTilde$sigmar[j] <- ifelse(Lr>1, var(classesr)*(varTilde$dl[j]-2)/varTilde$dl[j], sigmarValue[j])

        varTilde$delta[1:varTilde$k] <- deltaValue[ 1:kValue]
        varTilde$delta[j] <- tnormale(paramValues$zeta, 1/(varTilde$sigmaf[j]^{-1}+varTilde$sigmar[j]^{-1}), paramValues$deltamin, paramValues$deltamax)

        alpha <- rep(1, kValue)
        alphatilde <- rep(1,varTilde$k)
        ennetilde <- varTilde$dim[1:varTilde$k]
        alphaproptilde <- rep(1,varTilde$k)
        alphaprop <- rep(1, kValue)
        varTilde$w[1:varTilde$k] <- rdirichlet(1,alphaproptilde)

        ###calcul du rapport de vraisemblance de M-H move

        for (m in 1:varTilde$k) {
            Kaf[,m] <- (varTilde$w[m]*(1/sqrt(varTilde$sigmaf[m]))*dt(( paramValues$startPSF -varTilde$mu[m]+varTilde$delta[m]/2)/sqrt(varTilde$sigmaf[m]),varTilde$dl[m]))
            Kbf[,m] <- (wValue[m]*(1/sqrt(sigmafValue[m]))*dt(( paramValues$startPSF -muValue[m]+deltaValue[m]/2)/sqrt(sigmafValue[m]),dlValue[m]))
        }
        for (s in 1:paramValues$nf) {
            Y1f[s] <- log(sum(Kaf[s, 1:varTilde$k]))
            Y2f[s] <- log(sum(Kbf[s, 1:kValue]))
        }

        for (m in 1:varTilde$k) {
            Kar[,m] <- (varTilde$w[m]*(1/sqrt(varTilde$sigmar[m]))*dt(( paramValues$startPSR - varTilde$mu[m] - varTilde$delta[m]/2)/sqrt(varTilde$sigmar[m]),varTilde$dl[m]))
            Kbr[,m] <- (wValue[m]*(1/sqrt(sigmarValue[m]))*dt(( paramValues$startPSR - muValue[m]-deltaValue[m]/2)/sqrt(sigmarValue[m]),dlValue[m]))
        }
        for (s in 1:paramValues$nr) {
            Y1r[s] <- log(sum(Kar[s, 1:varTilde$k]))
            Y2r[s] <- log(sum(Kbr[s, 1:kValue]))
        }

        Kn <- sum(Y1f) + sum(Y1r)
        Kd <- sum(Y2f) + sum(Y2r)

        q         <- Kn - Kd
        rap.q     <- exp(q)

        rap.vrais <- rap.q

        rap.priormu   <- (priorMuDensity(varTilde$mu[1:varTilde$k],paramValues$y)/priorMuDensity(muValue[ 1:kValue],paramValues$y))
        rap.priorw    <- (ddirichlet(varTilde$w[1:varTilde$k], alphatilde)/ddirichlet(wValue[ 1:kValue],alpha))
        rap.priorenne <- dmultinom(ennetilde, paramValues$nbrReads, varTilde$w[1:varTilde$k])/dmultinom(dimValue[ 1:kValue], paramValues$nbrReads,wValue[ 1:kValue])
        rap.priork    <- (1)
        rap.propmu    <- (1)
        rap.propw     <- (ddirichlet(wValue[ 1:kValue],alphaprop)/ddirichlet(varTilde$w[1:varTilde$k],alphaproptilde))
        rap.prior     <- rap.priormu * rap.priorw * rap.priorenne * rap.priork
        rap.prop      <- rap.propmu  * rap.propw
        varTilde$rho      <- min(1, rap.vrais * (rap.prior) * (rap.prop ))

    }

    varTilde$rho <- ifelse(is.na(varTilde$rho) == FALSE, varTilde$rho, 0)
    return(varTilde)
}


deathMove <- function(paramValues , kValue, muValue, sigmafValue, sigmarValue, deltaValue, wValue, dlValue, aValue, dimValue ){
    ### Death move
    varTilde <- list( k=0L, mu=numeric(parValue$kmax)
                      , sigmaf=numeric(parValue$kmax), sigmar=numeric(parValue$kmax)
                      , delta=numeric(parValue$kmax), w=numeric(parValue$kmax)
                      , dl=numeric(parValue$kmax), a=numeric(parValue$kmax)
                      , dim=numeric(parValue$kmax), rho=0)

    Kn              <- rep(0, nbrIterations)
    Kaf             <- matrix(0, nrow = paramValues$nf, ncol = kmax)
    Kbf             <- matrix(0, nrow = paramValues$nf, ncol = kmax)
    Kar             <- matrix(0, nrow = paramValues$nr, ncol = kmax)
    Kbr             <- matrix(0, nrow = paramValues$nr, ncol = kmax)

    Y1f             <- rep(0, paramValues$nf)
    Y2f             <- rep(0, paramValues$nf)
    Y1r             <- rep(0, paramValues$nr)
    Y2r             <- rep(0, paramValues$nr)

    varTilde$k <- kValue - 1L
    count  <- 1L
    repeat {

        j <- sample(1:kValue,1)
        X <- muValue[ 1:kValue]
        varTilde$mu[1:varTilde$k] <- sort(c(X[-j]))

        Yf <- sigmafValue[1:kValue]
        varTilde$sigmaf[1:varTilde$k] <- Yf[-j]

        Yr <- sigmarValue[1:kValue]
        varTilde$sigmar[1:varTilde$k] <- Yr[-j]

        Delta <- deltaValue[1:kValue]
        varTilde$delta[1:varTilde$k] <- Delta[-j]

        Dl <- dlValue[1:kValue]
        varTilde$dl[1:varTilde$k] <- Dl[-j]

        Z <- aValue[1:(kValue+1)]
        if (j == kValue) {
            varTilde$a[1:(varTilde$k+1)] <- Z[-j]
        } else {
            varTilde$a[1:(varTilde$k+1)] <- Z[-(j+1)]
        }
        varTilde$a[1]             <- paramValues$minReadPos
        varTilde$a[varTilde$k+1]   <- paramValues$maxReadPos

        varTilde$dim[1] <- length(paramValues$y[varTilde$a[1]<=paramValues$y & paramValues$y<varTilde$a[2]])
        varTilde$dim[varTilde$k] <- length(paramValues$y[varTilde$a[varTilde$k]<=paramValues$y & paramValues$y<= paramValues$maxReadPos])
        if (varTilde$k>2) {
            for (m in 2: (varTilde$k-1)) {
                varTilde$dim[m] <- length(paramValues$y[varTilde$a[m]<=paramValues$y & paramValues$y<varTilde$a[m+1]]) }}
        Pr <- min(varTilde$dim[1:varTilde$k])

        ybar <- mean(paramValues$y[varTilde$a[j]<=paramValues$y & paramValues$y<=varTilde$a[j+1]])
        classesf <- paramValues$y[varTilde$a[j]<=paramValues$y & paramValues$y<=ybar]
        classesr <- paramValues$y[ybar<=paramValues$y & paramValues$y<=varTilde$a[j+1]]

        Lf <- length(classesf[!duplicated(classesf)])
        Lr <- length(classesr[!duplicated(classesr)])
        count <- count + 1L

        if ( (Pr>1 & Lf>1 & Lr>1)  ||
                 count == 1000L) break()
    }

    if (count == 1000L) {
        varTilde$rho <- 0
    } else {
        alpha <- rep(1, kValue)
        alphatilde <- rep(1, varTilde$k)
        alphaproptilde <- rep(1, varTilde$k)
        alphaprop <- rep(1, kValue)
        ennetilde <- varTilde$dim[1:varTilde$k]
        enne <- rmultinom(1, paramValues$nbrReads, wValue[1:kValue])
        varTilde$w[1:varTilde$k] <- rdirichlet(1, alphaproptilde)

        ### Rapport de vraisemblance ###

        for (m in 1:varTilde$k) {
            Kaf[,m] <- (varTilde$w[m]*(1/sqrt(varTilde$sigmaf[m]))*dt(( paramValues$startPSF -varTilde$mu[m]+varTilde$delta[m]/2)/sqrt(varTilde$sigmaf[m]),varTilde$dl[m]))
            Kbf[,m] <- (wValue[m]*(1/sqrt(sigmafValue[m]))*dt(( paramValues$startPSF -muValue[m]+deltaValue[m]/2)/sqrt(sigmafValue[m]),dlValue[m]))
        }
        Kbf[, kValue] <- (wValue[ kValue]*(1/sqrt(sigmafValue[ kValue]))*dt(( paramValues$startPSF -muValue[ kValue]+deltaValue[m]/2)/sqrt(sigmafValue[ kValue]),varTilde$dl[m]))

        for (s in 1:paramValues$nf) {
            Y1f[s] <- log(sum(Kaf[s, 1:varTilde$k]))
            Y2f[s] <- log(sum(Kbf[s, 1:kValue]))
        }

        for (m in 1:varTilde$k) {
            Kar[,m] <- (varTilde$w[m]*(1/sqrt(varTilde$sigmar[m]))*dt(( paramValues$startPSR -varTilde$mu[m]-varTilde$delta[m]/2)/sqrt(varTilde$sigmar[m]),varTilde$dl[m]))
            Kbr[,m] <- (wValue[m]*(1/sqrt(sigmarValue[m]))*dt(( paramValues$startPSR -muValue[m]-deltaValue[m]/2)/sqrt(sigmarValue[m]),dlValue[m]))
        }
        Kbr[, kValue] <- (wValue[ kValue]*(1/sqrt(sigmarValue[ kValue]))*dt((paramValues$startPSR -muValue[ kValue]-deltaValue[m]/2)/sqrt(sigmarValue[ kValue]),dlValue[m]))

        for (s in 1:paramValues$nr) {
            Y1r[s] <- log(sum(Kar[s, 1:varTilde$k]))
            Y2r[s] <- log(sum(Kbr[s, 1:kValue]))
        }

        Kn <- sum(Y1f) + sum(Y1r)
        Kd <- sum(Y2f) + sum(Y2r)

        q <- Kn - Kd
        rap.q <- exp(q)

        rap.vrais <- rap.q

        # Density of varTilde$mu[j]
        if (j == 1) {
            qalloc <- 1/(muValue[j+1] - paramValues$minReadPos)
        } else {
            if (j == kValue) {
                qalloc <- 1/(paramValues$maxReadPos - muValue[ j-1])
            } else {
                qalloc <- 1/(muValue[ j+1] - muValue[ j-1])
            }
        }

        rap.priormu    <- (priorMuDensity(varTilde$mu[1:varTilde$k],paramValues$y)/priorMuDensity(muValue[ 1:kValue], paramValues$y))
        rap.priorw     <- (ddirichlet(varTilde$w[1:varTilde$k],alphatilde)/ddirichlet(wValue[ 1:kValue], alpha))
        rap.priorenne  <- dmultinom(ennetilde, paramValues$nbrReads,varTilde$w[1:varTilde$k])/dmultinom(dimValue[ 1:kValue], paramValues$nbrReads,wValue[ 1:kValue])
        rap.priork     <- (dpois(varTilde$k,paramValues$lambda)/dpois(kValue,paramValues$lambda))
        rap.propmu     <- (qalloc)
        rap.propw      <- (ddirichlet(varTilde$w[ 1:varTilde$k], alphaproptilde)/ddirichlet(wValue[ 1:kValue], alphaprop))

        rap.prior      <- rap.priormu * rap.priorw * rap.priorenne * rap.priork
        rap.prop       <- rap.propmu  * rap.propw

        varTilde$rho        <- min(1,(rap.vrais) * (rap.prior)  *  (rap.prop) * (Bk(varTilde$k,paramValues$lambda,paramValues$kmax)/Dk(kValue, paramValues$lambda, paramValues$kmax)))

    }

    varTilde$rho <-  ifelse( is.na(varTilde$rho) == FALSE, varTilde$rho, 0)
    return(varTilde)
}

