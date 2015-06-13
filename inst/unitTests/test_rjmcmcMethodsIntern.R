###################################################
# Created by Astrid Louise Deschenes
# 2015-06-12
###################################################

###################################################
## Test the findConsensusPeakRegions.R functions
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "rjmcmc" )
}

### }}}

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
                    0.5000000, 0.0000000, 0.5000000)
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

# test.tnormale_result_with_various_values_of_mu <- function() {
#     set.seed(101)
#     obs <- mapply(c(9, 8, 6, 2, 24, 3, 1),
#                   FUN = function(x) {rjmcmc:::tnormale(x, 7, 7, 10)})
#     exp <- c(7.0218905285, 9.5625980973, 8.5250712962, 8.4394940419,
#              7.4477169047, 9.5883395894, 7.7767406354)
#     message <- paste0(" tnormale_result_with_various_values_of_lambda() ",
#                       "- Not all tested data with various lambda generated",
#                       " expected values.")
#     checkEqualsNumeric(obs, exp, msg = message)
# }


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







