###################################################
# Created by Astrid Louise Deschenes
# 2015-06-12
###################################################

###################################################
## Test the rjmcmcMethodsIntern.R functions
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "rjmcmc" )
}

### }}}

###########################################################
## Dk() function
###########################################################

test.Dk_result_k_inferior_to_kmax <- function() {
    obs <- mapply(c(8, 6, 5, 4, 3, 1, 0, -1),
                    FUN = function(x) {rjmcmc:::Dk(9, 8, x)})
    message <- paste0(" Dk_result_k_inferior_to_kmax() ",
                    "- A kmax inferior to k did not always returned 0")
    checkIdentical(obs, rep(0, 8), msg = message)
}

test.Dk_result_with_various_values_of_k <- function() {
    obs <- mapply(c(9, 8, 6, 5, 4, 3, 1, 0),
                    FUN = function(x) {rjmcmc:::Dk(x, 8, 10)})
    exp <- c(0.5000, 0.5000, 0.3750, 0.3125, 0.2500, 0.1875,
                    0.0000, 0.0000)
    message <- paste0(" Dk_result_with_various_values_of_k() ",
                    "- Not all tested data generated expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

test.Dk_result_with_various_values_of_lambda <- function() {
    obs <- mapply(c(9, 8, 6, 35, 24, 13, 1),
                    FUN = function(x) {rjmcmc:::Dk(7, x, 10)})
    exp <- c(0.38888888888888, 0.4375000, 0.5000000, 0.1000000,
                    0.14583333333, 0.2692307692, 0.5000000)
    message <- paste0(" Dk_result_with_various_values_of_lambda() ",
                    "- Not all tested data generated expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

###########################################################
## Bk() function
###########################################################

test.Bk_result_k_inferior_to_kmax <- function() {
    obs <- mapply(c(8, 6, 5, 4, 3, 1, 0, -1),
                  FUN = function(x) {rjmcmc:::Bk(9, 8, x)})
    message <- paste0(" Bk_result_k_inferior_to_kmax() ",
                    "- A kmax inferior to k did not always returned 0")
    checkIdentical(obs, rep(0, 8), msg = message)
}

test.Bk_result_with_various_values_of_k <- function() {
    obs <- mapply(c(9, 8, 6, 5, 4, 3, 1, 0),
                    FUN = function(x) {rjmcmc:::Bk(x, 8, 10)})
    exp <- c(0.4000000, 0.444444444444444, 0.5000000, 0.5000000, 0.5000000,
                    0.5000000, 0.5000000, 0.5000000)
    message <- paste0(" Bk_result_with_various_values_of_k() ",
                    "- Not all tested data generated expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

test.Bk_result_with_various_values_of_lambda <- function() {
    obs <- mapply(c(9, 8, 6, 2, 24, 3, 1),
                  FUN = function(x) {rjmcmc:::Bk(7, x, 10)})
    exp <- c(0.5000, 0.5000, 0.3750, 0.1250, 0.5000, 0.1875, 0.0625)
    message <- paste0(" Bk_result_with_various_values_of_lambda() ",
                    "- Not all tested data generated expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

###########################################################
## tnormale() function
###########################################################

test.tnormale_result_with_various_values_of_mu <- function() {
    set.seed(101)
    obs <- mapply(c(9, 8, 6, 2, 24, 3, 1),
                  FUN = function(x) {rjmcmc:::tnormale(x, 7, x-1, x+1)})
    exp <- c(8.1373885, 8.5671418, 6.8222181, 1.7017330, 23.4093112,
             2.3737946,  0.48847584)
    message <- paste0(" tnormale_result_with_various_values_ofm_mu() ",
                    "- Not all tested data with various mu generated",
                    " expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

test.tnormale_result_with_various_values_of_a <- function() {
    set.seed(101)
    obs <- mapply(c(1, 2, 3, 4, 5, 6, 7),
                  FUN = function(x) {rjmcmc:::tnormale(8, 10090, x, 10)})
    exp <- c(2.2091408741, 2.8786961586, 8.5411309921, 4.5581059972,
             7.8565459593, 6.2785794538, 9.8004753529)
    message <- paste0(" tnormale_result_with_various_values_of_a() ",
                    "- Not all tested data with various a generated",
                    " expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

test.tnormale_result_with_various_values_of_b <- function() {
    set.seed(101)
    obs <- mapply(c(11, 12, 13, 14, 15, 16, 17),
                  FUN = function(x) {rjmcmc:::tnormale(10, 10090, 10, x)})
    exp <- c(10.5411309921, 11.8004753529, 12.4993162124, 10.5176897766,
             11.2364951059, 15.3437143272, 15.1289317804)
    message <- paste0(" tnormale_result_with_various_values_of_b() ",
                    "- Not all tested data with various b generated",
                    " expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

test.tnormale_result_with_various_values_of_lambda <- function() {
    set.seed(101)
    obs <- mapply(c(9, 8, 6, 2, 24, 3, 1),
                  FUN = function(x) {rjmcmc:::tnormale(8, x, 7, 10)})
    exp <- c(7.0218905285, 9.5625980973, 8.5250712962, 8.4394940419,
             7.4477169047, 9.5883395894, 7.7767406354)
    message <- paste0(" tnormale_result_with_various_values_of_lambda() ",
                    "- Not all tested data with various lambda generated",
                    " expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

#########################################################
## elementWithHighestMode() function
#########################################################

test.elementWithHighestMode_results_with_different_vectors <- function() {
    set.seed(101)
    obs <- lapply(list(test1 = c(1, 2, 3, 3, 21, 22),
                    test2 = c(2, 4, 2, 4, 4, 0),
                    test3 = c(12, 13, 13, 12, 1),
                    test4 = c(23, 22, 21, 20, 22, 22, 21, 21, 22)),
                  FUN = function(x) {rjmcmc:::elementWithHighestMode(x)})
    exp <- list(test1 = 3L, test2 = 4L, test3 = NA, test4 = 22L)
    message <- paste0(" elementWithHighestMode_results_with_different_vectors() ",
                    "- Not all tested vectors generated",
                    " expected values.")
    checkEquals(obs, exp, msg = message)
}

#########################################################
## priorMuDensity() function
#########################################################

test.priorMuDensity_with_1000_reads <- function() {
    set.seed(101)
    k <- 3L
    nbrReads <- 500L
    mu <- c(10000L, 26700L, 45000L)
    sigma <- rep(400L, k)
    delta  <- rep(147, k)
    weight <- c(0.3, 0.2, 0.5)
    readsForward <- sapply(1:nbrReads, rjmcmc:::normal.mixture, k = k,
                        w = weight, mu = mu - delta/2, sigma = sigma)
    readsReverse <- sapply(1:nbrReads, rjmcmc:::normal.mixture, k = k,
                        w = weight, mu = mu - delta/2, sigma = sigma)
    reads <- sort(c(readsForward, readsReverse))
    obs <- rjmcmc:::priorMuDensity(mu, reads)
    exp <- 8.0883349761e-15
    message <- paste0(" priorMuDensity_with_1000_reads() ",
                    "- The result is not the expected value.")
    checkEqualsNumeric(obs, exp, msg = message)
}

test.priorMuDensity_results_with_various_values_of_mu <- function() {
    set.seed(101)
    obs <- mapply(list(A=c(10200, 10300), B=c(10108, 10206, 10222),
                    C=c(10333, 10455, 10899)),
                    FUN = function(x) { rjmcmc:::priorMuDensity(x,
                                            10000:11000) })
    exp <- c(6.0557145969e-7, 4.6807063829e-10, 4.5053092518e-10)
    message <- paste0(" priorMuDensity_results_with_various_values_of_mu() ",
                    "- Not all tested data with various mu generated",
                    " expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

#########################################################
## normal.mixture() function
#########################################################

test.normal_mixture_results_with_various_values_of_mu <- function() {
    set.seed(101)
    k <- 6
    weight <- c(0.1, 0.15, 0.15, 0.2, 0.4)
    sigma <- c(4,4,4,4,4,4)
    mu <- list(A=c(4,6,7,8,34,44), B=c(2,16,17,21,24,34),
            C=c(102, 103, 106, 200, 201, 222))
    obs <- mapply(mu,
                    FUN = function(x) { rjmcmc:::normal.mixture(k = k,
                    weight = weight, mu = x, sigma = sigma)} )
    exp <- c(3.58414331129, 22.65011231209, 199.13928138368)
    message <- paste0(" normal_mixture_results_with_various_values_of_mu() ",
                      "- Not all tested data with various mu generated",
                      " expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

#########################################################
## student.mixture() function
#########################################################

test.student_mixture_results_with_various_values_of_mu <- function() {
    set.seed(101)
    k <- 6
    weight <- c(0.1, 0.15, 0.15, 0.2, 0.4)
    sigma <- c(4,4,4,4,4,4)
    dfr <- sigma
    mu <- list(A=c(4,6,7,8,34,44), B=c(2,16,17,21,24,34),
               C=c(102, 103, 106, 200, 201, 222))
    obs <- mapply(mu,
                  FUN = function(x) { rjmcmc:::student.mixture(k = k,
                                        weight = weight, mu = x,
                                        sigma = sigma, dfr = dfr)} )
    exp <- c(3.61671572308, 17.43930668620, 201.78192965619)
    message <- paste0(" student_mixture_results_with_various_values_of_mu() ",
                      "- Not all tested data with various mu generated",
                      " expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

#########################################################
## isInteger() function
#########################################################

## Test the result when a integer is passed to the function
test.isInteger_with_integer <- function() {
    obs <- tryCatch(rjmcmc:::isInteger(value = 444L),
                    error=conditionMessage)
    message <- paste0(" isInteger_with_integer() ",
                      "- A integer did not generated TRUE.")
    checkEquals(obs, TRUE, msg = message)
}

## Test the result when a numeric is passed to the function
test.isInteger_with_numeric <- function() {
    obs <- tryCatch(rjmcmc:::isInteger(value = 444.2),
                    error=conditionMessage)
    message <- paste0(" isInteger_with_numeric() ",
                      "- A numeric did not generated TRUE.")
    checkEquals(obs, TRUE, msg = message)
}


## Test the result when a string is passed to the function
test.isInteger_with_string <- function() {
    obs <- tryCatch(rjmcmc:::isInteger(value = "444"),
                    error=conditionMessage)
    message <- paste0(" isInteger_with_string() ",
                      "- A string did not generated FALSE.")
    checkEquals(obs, FALSE, msg = message)
}

## Test the result when a list is passed to the function
test.isInteger_with_list_of_integers <- function() {
    obs <- tryCatch(rjmcmc:::isInteger(value = list(a=444L, b=333L)),
                    error=conditionMessage)
    message <- paste0(" isInteger_with_list_of_integers() ",
                      "- A list did not generated FALSE.")
    checkEquals(obs, FALSE, msg = message)
}

## Test the result when a vector is passed to the function
test.isInteger_with_vector_of_integers <- function() {
    obs <- tryCatch(rjmcmc:::isInteger(value = c(444L, 333L)),
                    error=conditionMessage)
    message <- paste0(" isInteger_with_vector_of_integers() ",
                      "- A vector of integers did not generated FALSE.")
    checkEquals(obs, FALSE, msg = message)
}
