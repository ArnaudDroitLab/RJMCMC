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
#' @author Rawane Samb, Astrid Louise Deschenes
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
#' @author Rawane Samb, Astrid Louise Deschenes
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
#' @author Rawane Samb, Astrid Louise Deschenes
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
#' @author Rawane Samb, Astrid Louise Deschenes
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
#' @param liste a \code{list} containing:
#' \itemize{
#'     \item k a \code{integer}, the number of nucleosomes.
#'     \item mu a \code{vector} of \code{numeric}, the positions of
#' the nucleosomes.
#'     \item sigmaf TODO
#'     \item sigmar TODO
#'     \item delta TODO
#'     \item dl TODO
#'     \item w TODO
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
#'     \item sigmar TODO
#'     \item delta TODO
#'     \item dl TODO
#'     \item w TODO
#' }
#'
#' @author Rawane Samb, Astrid Louise Deschenes
#' @keywords internal
mergeNucleosomes <- function(yf, yr, y, liste,
                                minInterval, maxInterval, minReads) {
    ## Get the number of nucleosomes
    k <- liste$k

    ## Nucleosome can only be merged when there is more than one
    if (k > 1) {

#         ecart.min <- min(sapply(1:(k-1),
#                             function(j){liste$mu[j+1] - liste$mu[j]}))

        ## Find the smallest distance between 2 nucleosomes
        ecart <- diff(liste$mu)
        ecart.min <- min(ecart)

        ## The distance between 2 nucleosomes is inferior to minimum distance
        ## allowed
        if (ecart.min < minInterval) {
            ## Merging nucleosomes with distance inferior to minimum distance
            ## up to the moment there is no nucleosome with inferior distance
            repeat
            {
#                 p <- which(sapply(1:(k-1),
#                             function(j){liste$mu[j+1] - liste$mu[j]}) ==
#                             ecart.min)[1]

                ## Find the first position with minimum gap
                p <- which.min(ecart)

                classes  <- y[y >= liste$mu[p] & y < liste$mu[p+1]]
                classesf <- yf[yf >= liste$mu[p] & yf < liste$mu[p+1]]
                classesr <- yr[yr >= liste$mu[p] & yr < liste$mu[p+1]]

                if (length(classes) > minReads){
                    mu <- mean(round(classes))
                } else {
                    mu <- mean(c(liste$mu[p], liste$mu[p+1]))
                }

                liste$mu        <- sort(liste$mu[-p])
                liste$mu[p]     <- mu
                liste$mu        <- sort(liste$mu)
                liste$sigmaf    <- liste$sigmaf[-p]
                liste$sigmar    <- liste$sigmar[-p]
                liste$delta     <- liste$delta[-p]
                liste$dl        <- liste$dl[-p]
                liste$w         <- liste$w[-p]/sum(liste$w[-p])

                # Downgrade the number of nucleosome
                k <- k - 1

                # Update the smallest distance between 2 nucleosomes
                if (k > 1) {
#                     ecart.min <- min(sapply(1:(k-1),
#                                     function(i){liste$mu[i+1]-liste$mu[i]}))
                    ecart <- diff(liste$mu)
                    ecart.min <- min(ecart)
                }
                if (k == 1 || ecart.min > minInterval) break()
            } ### end of boucle repeat

            ## Updating resulting list
            liste <- list(  k      = k,
                            mu     = liste$mu,
                            sigmaf = liste$sigmaf,
                            sigmar = liste$sigmar,
                            delta  = liste$delta,
                            dl     = liste$dl,
                            w      = liste$w)
        } ### end of condition if (ecart.min < minInterval)
#         else {
#             liste <- liste
#         }

        ## Trying to split resulting nucleosomes
        liste <- splitNucleosome(yf, yr, y, liste, minInterval,
                                    maxInterval, minReads)
    } ### end of condition if (k > 1)

    return(liste)
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
#' @param liste a \code{list} containing:
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
#' @author Rawane Samb, Astrid Louise Deschenes
#' @keywords internal
splitNucleosome <- function(yf, yr, y, liste, minInterval, maxInterval,
                                    minReads) {
    ## Get the number of nucleosomes
    k <- liste$k

    if (k > 1) {
#         ecart.max <- max(sapply(1:(k-1),
#                                 function(j){liste$mu[j+1]-liste$mu[j]}))

        ## Find the largest distance between 2 nucleosomes
        ## Its position and its value
        ecart       <- diff(liste$mu)
        p           <- order(ecart, decreasing = TRUE)[1]
        ecart.max   <- ecart[p]

        if (ecart.max > maxInterval) {
            j <- 1
            repeat {
#                 p <- which(sapply(1:(k-1),
#                                 function(j){
#                                     liste$mu[j+1]-liste$mu[j]
#                                 }) == ecart.max)
#                 p <- which(sapply(1:(k-1),
#                                   function(j){
#                                       liste$mu[j+1]-liste$mu[j]
#                                   }) == ecart.max)[1]

                classes <- y[y>=liste$mu[p] & y<liste$mu[p+1]]
                classesf <- yf[yf>=liste$mu[p] & yf<liste$mu[p+1]]
                classesr <- yr[yr>=liste$mu[p] & yr<liste$mu[p+1]]
#                 j <- 1
                if (length(classes) > minReads)
                {
                    new.mu <- sort(c(liste$mu[1:k],mean(round(classes))))
                    new.sigmaf <- c(liste$sigmaf[1:k],
                                        (liste$sigmaf[p]+liste$sigmaf[p+1])/2)
                    new.sigmaf[p+1] <- (liste$sigmaf[p]+liste$sigmaf[p+1])/2
                    new.sigmaf[k+1] <- liste$sigmaf[k]
                    new.sigmar <- c(liste$sigmar[1:k],
                                        (liste$sigmar[p]+liste$sigmar[p+1])/2)
                    new.sigmar[p+1] <- (liste$sigmar[p]+liste$sigmar[p+1])/2
                    new.sigmar[k+1] <- liste$sigmar[k]
                    new.delta <- c(liste$delta[1:k],
                                   (liste$delta[p]+liste$delta[p+1])/2)
                    new.delta[p+1] <- (liste$delta[p]+liste$delta[p+1])/2
                    new.delta[k+1] <- liste$delta[k]
                    new.dl <- round(c(liste$dl[1:k],
                                        (liste$dl[p]+liste$dl[p+1])/2))
                    new.dl[p+1] <- (liste$dl[p]+liste$dl[p+1])/2
                    new.dl[k+1] <- liste$dl[k]
                    new.w <- c(liste$w[1:k], (liste$w[p]+liste$w[p+1])/2)
                    new.w[p+1] <- (liste$w[p]+liste$w[p+1])/2
                    new.w[k+1] <- liste$w[k]
                    k <- length(new.mu)
                    liste <- list(  k       = k,
                                    mu      = new.mu,
                                    sigmaf  = new.sigmaf,
                                    sigmar  = new.sigmar,
                                    delta   = new.delta,
                                    dl      = new.dl,
                                    w       = new.w/sum(new.w))
#                     ecart.max <- max(sapply(1:(k-1),
#                                     function(j){liste$mu[j+1]-liste$mu[j]}))

#                     j           <- 1

                    ## Update the vector of distance between nucleosomes
                    ecart <- diff(liste$mu)
                }
                else
                {
#                     liste <- liste
#                     ecart.max <- sort(sapply(1:(k-1),
#                                     function(j){
#                                         liste$mu[j+1]-liste$mu[j]
#                                     }))[k - 1 - j]

                    ## Update to select the next maximum distance value
                    j  <- j + 1
                }

                ## Select the next nucleosome to be potentially split
                p           <- order(ecart, decreasing = TRUE)[j]
                ecart.max   <- ecart[p]

                if ( j >= (k - 1) || ecart.max <= maxInterval) break()
#                 if ( j == (k - 1) || ecart.max <= maxInterval) break()
            } ### end of boucle repeat
        } ### end of condition if (ecart.max > maxInterval)
    } ### end of condition if (k > 1)
    return(liste)
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
#' zero. The maximum value of \code{nbrIterations} is \code{100000}.
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
#' @return \code{0} indicating that all parameters validations have been
#' successful.
#'
#' @author Astrid Louise Deschenes
#' @keywords internal
validateParameters <- function(startPosForwardReads, startPosReverseReads,
                                    nbrIterations, kmax, lambda,
                                    minInterval, maxInterval, minReads) {
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
        || length(startPosForwardReads) < 1 || !all(startPosForwardReads > 0))
    {
        stop(paste0("startPosForwardReads must be a non-empty vector of ",
                    "non-negative numeric values."))
    }

    ## Validate that the startPosReverseReads has at least one read
    ## and that the values are integer
    if (!is.vector(startPosReverseReads) || !is.numeric(startPosReverseReads)
        || length(startPosReverseReads) < 1 || !all(startPosReverseReads > 0))
    {
        stop(paste0("startPosReverseReads must be a non-empty vector of ",
                    "non-negative numeric values."))
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
#' @author Astrid Louise Deschenes
#' @keywords internal
#'
isInteger <- function(value) {
    return((is.integer(value) && length(value) == 1) || (is.numeric(value) &&
        length(value) == 1))
}

