#' @title Death Submove Probability
#'
#' @description Calculation of the death submove
#'      probability of a randomly selected nucleosome using
#'      a truncated Poisson distribution.
#'
#' @param k a positive \code{integer}, the number of nucleosomes.
#'
#' @param lambda a \code{numeric}, the theorical mean
#'      of the Poisson distribution.
#'
#' @param kMax a positive \code{numeric}, the maximum number of nucleosomes
#'      authorized. When
#'      \code{k} is equal or superior to \code{kMax}, the
#'      returned value is \code{0}.
#'
#' @return a \code{numeric} value representing the calculated death submove
#'      probability. The value \code{0} when
#'      \code{k} is equal or superior to \code{kMax} or
#'      when \code{k} is equal to \code{1}.
#' @examples
#'
#' ## Return the death submove probability
#' rjmcmc:::Dk(k = 12L, lambda = 4.4, kMax = 30L)
#'
#' ## Zero is returned when k = 1
#' rjmcmc:::Dk(k = 1L, lambda = 4.4, kMax = 30L)
#'
#' ## Zerio is returned when k is superior to kMax
#' rjmcmc:::Dk(k = 31L, lambda = 4.4, kMax = 30L)
#'
#' @author Rawane Samb
#' @importFrom stats dpois
#' @keywords internal
Dk <- function(k, lambda, kMax) {
    ifelse((k == 1 || k > kMax), 0,
            0.5*min(1, dpois(k-1,lambda)/dpois(k, lambda)))
}


#' @title Birth Submove Probability
#'
#' @description Calculation of the birth submove
#'      probability of adding a new nucleosome using a
#'      truncated Poisson distribution.
#'
#' @param k a positive \code{integer}, the number of
#'      nucleosomes.
#'
#' @param lambda a \code{numeric}, the theorical mean
#'      of the Poisson distribution.
#'
#' @param kMax a positive \code{integer}, the maximum number of nucleosomes
#'      authorized. When
#'      \code{k} is equal or superior to \code{kMax}, the
#'      returned value is \code{0}.
#'
#' @return a \code{numeric} value. The value \code{0} when
#'      \code{k} is equal or superior to \code{kMax} or
#'      when \code{k} is equal to \code{1}.
#' @examples
#'
#' ## Return the birth submove probability
#' rjmcmc:::Bk(k = 14L, lambda = 2.4, kMax = 22L)
#'
#' ## Zero is returned when k = 1
#' rjmcmc:::Bk(k = 1L, lambda = 3.4, kMax = 20L)
#'
#' ## Zerio is returned when k is superior to kMax
#' rjmcmc:::Bk(k = 31L, lambda = 4.4, kMax = 30L)
#'
#' @author Rawane Samb
#' @importFrom stats dpois
#' @keywords internal
Bk <- function(k, lambda, kMax) {
    ifelse((k == 1 || k > kMax), 0,
           0.5 * min(1, dpois(k + 1, lambda) / dpois(k, lambda)))
}


#' @title Truncated random deviate from a normal distribution
#'
#' @description Generate a random deviate value from a normal
#'    distribution. The returned value is included inside a
#'    sp
#'    ecified range ]minValue,maxValue[ specified by user. The mean and
#'    variance of the normal distribution is also specified by
#'    user.
#'
#' @param mu a \code{numeric} mean of the normal distribution.
#'
#' @param sigma a \code{numeric} variance of the normal distribution.
#'
#' @param minValue a \code{numeric} inferior boundary of the range
#'      in which the output value must be located. The
#'      output value has to be superior to \code{minValue}.
#'
#' @param maxValue a \code{numeric} superior boundary of the range
#'      in which the output value must be located. The
#'      output value has to be inferior to \code{maxValue}.
#'
#' @return a \code{numeric} superior to \code{minValue} and inferior
#'      to \code{maxValue}.
#'
#' @author Rawane Samb
#' @importFrom stats rnorm
#' @keywords internal
tnormale <- function(mu, sigma, minValue, maxValue)
{
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
#' @description TODO
#'
#' @param i a \code{vector} TODO
#'
#' @param k a \code{integer} number of nucleosomes in a region.
#'
#' @param w a \code{vector} weight for the nucleosome occupancy.
#'
#' @param mu a \code{vector} mean of the Student-t distribution.
#'
#' @param sigma a \code{vector} TODO
#'
#' @param dfr a \code{vector} of \code{numeric} containing the degree
#'      of freedom.
#'
#' @return \code{0} TODO
#'
#' @author Rawane Samb
#' @keywords internal
student.mixture <- function(i, k, w, mu, sigma, dfr)
{
    ## TODO : valider si k doit petre un entier
    ## TODO : i n'est pas utilisé
    v <- c(0, w)
    u <- runif(1, 0, 1)
    for (j in 1:k) {
        if(sum(v[1:j]) < u & u <= sum(v[1:(j+1)])) {
            mixte <- mu[j] + sqrt(sigma[j]) * rt(1, dfr[j])
        }
    }
    return(mixte)
}


#' @title Normal Mixture Model
#'
#' @description TODO
#'
#' @param i a \code{vector} TODO
#'
#' @param k a \code{vector} TODO
#'
#' @param w a \code{vector} TODO
#'
#' @param mu a \code{vector} TODO
#'
#' @param sigma a \code{vector} TODO
#'
#' @return \code{0} TODO
#'
#' @author Rawane Samb
#' @keywords internal
normal.mixture <- function(i, k, w, mu, sigma)
{
    v <- c(0, w)
    u <- runif(1, 0, 1)
    for (j in 1:k) {
        if (sum(v[1:j]) < u & u <= sum(v[1:(j+1)]))
        {
            mixte <- rnorm(1, mu[j], sd = sqrt(sigma[j]))
        }
    }
    return(mixte)
}


#' @title Prior density of \eqn{mu}
#'
#' @description Computes the prior density of \eqn{mu} conditionally to
#'  the number of nucleosomes.
#'
#'  For more information on the calculation of the prior density of \eqn{mu},
#'  see Proposotion 1 of the cited article.
#'
#' @param mu a \code{vector} of positive \code{integer} containing the
#'      positions of all nucleosomes.
#'
#' @param reads a \code{vector} of \code{TODO} corresponding to the read
#'      data, including forward and reverse strands.
#'
#' @return  the exact prior density of \code{mu} given the
#'      number of nucleosomes.
#'
#' @references Samb R., Khadraoui K., Lakhal L., Belleau P. and Droit A. Using
#'      informative Multinomial-Dirichlet prior in a t-mixture with
#'      reversible jump estimation of nucleosome positions for genome-wide
#'      profiling. Submitted (2015).
#' @author Rawane Samb, Astrid Louise Deschenes
#' @keywords internal
priormu <- function(mu, reads)
{
    k <- length(mu)
    ## Create a matrix used in the calculation of the priors
    basicMatrix <- matrix(0L, nrow = k, ncol = k)
    for (i in 1:k) {
        for (j in 1:k) {
            if (j == i) {
                basicMatrix[i, j] <- 1L
            } else if (j == i - 1) {
                basicMatrix[i, j] <- -1L
            }
        }
    }
    omega <- t(basicMatrix) %*% basicMatrix
    # Calculating the range (R)
    R <- max(reads) - min(reads)
    # Calculating the mean (E)
    E <- (max(reads) + min(reads))/2
    tau <- 1/R^2
    M <- rep(E, k)
    const <- (pi/(2*tau))^{-k/2}
    prior <- const * exp(-(tau/2) * (t(mu - M) %*% omega %*% (mu - M)))
    return(prior)
}



#' @title Element with the hightest number of occurences
#'
#' @description \code{mode} takes the integer-valued vector \code{sample} and
#'      returned the \code{integer} with the highest number of occurences.
#'      When more than one \code{integer} have the highest number of
#'      occurences, \code{NA} is returned.
#'
#' @param sample a \code{numeric} \code{vector} (of positive \code{integer}
#'      values). If the elements of \code{sample} are \code{numeric} but not
#'      \code{integer}, the elements are truncated by \code{as.integer}.

#' @return  a \code{integer} with the highest number of occurences or
#'      \code{NA} when more than one \code{integer} have the highest number
#'      of occurences.
#'
#' @author Rawane Samb, Astrid Louise Deschenes
#' @keywords internal
#' @examples
#'
#' ## Return the element with the hightest number of occurence
#' data01 <- c(1L, 2L, 5L, 10L, 5L, 10L, 5L)
#' rjmcmc:::mode(data01)
#'
#' data02 <- c(3L, 6L, 4L, 3L, 6L)
#' rjmcmc:::mode(data02)
#'
mode <- function(sample) {
    tabsample <- tabulate(sample)
    #samplemode <- which(tabsample == max(tabsample))[1]
    #if(sum(tabsample == max(tabsample)) > 1) {
    #    samplemode <- NA
    #}
    ## CODE MODIFIER PAR ASTRID
    maxOccurence <- tabsample == max(tabsample)
    ifelse(sum(maxOccurence) == 1, which(maxOccurence), NA)
}


#' @title Merging two nucleosomal regions
#'
#' @description TODO
#'
#' @param yf TODO
#'
#' @param yr TODO
#'
#' @param y TODO
#'
#' @param liste TODO
#'
#' @param ecartmin TODO
#'
#' @param ecartmax TODO
#'
#' @param minReads TODO
#'
#' @return \code{0} TODO
#'
#' @author Rawane Samb
#' @keywords internal
merge <- function(yf, yr, y, liste, ecartmin, ecartmax, minReads)
{
    k <- length(liste$mu)
    if (k > 1) {
        ecart.min <- min(sapply(1:(k-1),
                            function(j){liste$mu[j+1] - liste$mu[j]}))
        if (ecart.min < ecartmin)
        {
            repeat
            {
                p <- which(sapply(1:(k-1),
                            function(j){liste$mu[j+1] - liste$mu[j]}) ==
                            ecart.min)[1]

                classes <- y[y >= liste$mu[p] & y < liste$mu[p+1]]
                classesf <- yf[yf >= liste$mu[p] & yf < liste$mu[p+1]]
                classesr <- yr[yr >= liste$mu[p] & yr < liste$mu[p+1]]

                if (length(classes) > minReads){
                    mu <- mean(round(classes))
                } else {
                    mu <- mean(c(liste$mu[p], liste$mu[p+1]))
                }

                liste$mu <- sort(liste$mu[-p])
                liste$mu[p] <- mu
                liste$mu <- sort(liste$mu)
                liste$sigmaf <- liste$sigmaf[-p]
                liste$sigmar <- liste$sigmar[-p]
                liste$delta <- liste$delta[-p]
                liste$dl <- liste$dl[-p]
                liste$w <- liste$w[-p]/sum(liste$w[-p])
                k <- k-1
                if (k > 1) {
                    ecart.min <- min(sapply(1:(k-1),
                                    function(i){liste$mu[i+1]-liste$mu[i]}))}
                if (k == 1 || ecart.min > ecartmin) break()
            } ### end of boucle repeat
            liste <- list(k = k,
                            mu = liste$mu,
                            sigmaf = liste$sigmaf,
                            sigmar = liste$sigmar,
                            delta = liste$delta,
                            dl = liste$dl,
                            w = liste$w)
        } ### end of condition if (ecart.min < ecartmin)
        else {
            liste <- liste
        }
        liste <- split(yf, yr, y, liste, ecartmin, ecartmax, minReads)
    } ### condition else if (k > 1)
    return(liste)
}


#' @title Spliting a nucleosomal region into two regions
#'
#' @description TODO
#'
#' @param yf TODO
#'
#' @param yr TODO
#'
#' @param y TODO
#'
#' @param liste TODO
#'
#' @param ecartmin TODO
#'
#' @param ecartmax TODO
#'
#' @param minReads TODO
#'
#' @return \code{0} TODO
#'
#' @author Rawane Samb
#' @keywords internal
split <- function(yf, yr, y, liste, ecartmin, ecartmax, minReads)
{
    ## Astrid : on devrait changer le nom car porte a confusion avec
    ## les fonctions split() existantes
    k <- length(liste$mu)
    if (k>1) {
        ecart.max <- max(sapply(1:(k-1),
                                function(j){liste$mu[j+1]-liste$mu[j]}))
        if (ecart.max > ecartmax) {
            j <- 1
            repeat {
                p <- which(sapply(1:(k-1),
                                function(j){
                                    liste$mu[j+1]-liste$mu[j]
                                }) == ecart.max)

                classes <- y[y>=liste$mu[p] & y<liste$mu[p+1]]
                classesf <- yf[yf>=liste$mu[p] & yf<liste$mu[p+1]]
                classesr <- yr[yr>=liste$mu[p] & yr<liste$mu[p+1]]
                j <- 1
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
                    liste <- list(k = k,
                                    mu = new.mu,
                                    sigmaf = new.sigmaf,
                                    sigmar = new.sigmar,
                                    delta = new.delta,
                                    dl = new.dl,
                                    w = new.w/sum(new.w))
                    ecart.max <- max(sapply(1:(k-1),
                                    function(j){liste$mu[j+1]-liste$mu[j]}))
                }
                else
                {
                    liste <- liste
                    ecart.max <- sort(sapply(1:(k-1),
                                    function(j){
                                        liste$mu[j+1]-liste$mu[j]
                                    }))[k - 1 - j]
                    j <- j + 1
                }
                if ( j == (k - 1) || ecart.max <= ecartmax) break()
            } ### end of boucle repeat
        } ### end of condition if (ecart.max > ecartmax)
    } ### end of condition if (k>1)
    return(liste)
}
