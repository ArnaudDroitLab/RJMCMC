###################################################
# Created by Astrid Louise Deschenes
# 2015-06-30
###################################################

###################################################
## Test the rjmcmcMethod.R functions
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "rjmcmc" )
}

### }}}

data(reads_demo)


###########################################################
## RJMCMC() function
###########################################################

test.rjmcmc_good_result_01 <- function() {
    set.seed(101)
    obs <- RJMCMC(startPosForwardReads = reads_demo$readsForward,
                        startPosReverseReads = reads_demo$readsReverse,
                        nbrIterations = 100, lambda = 2, kmax = 30,
                        minInterval = 146, maxInterval = 292, minReads = 5)
    exp.k       <- 2
    exp.mu      <- c(72669.922485424002, 72898.810781247870)
    exp.sigmaf  <- c(21509.940563849304, 4736.790424950710)
    exp.sigmar  <- c(12509.376406312418, 5617.742233704413)
    exp.delta   <- c(144.376399340025, 142.465706790323)
    exp.dl      <- c(13, 3)
    exp.w       <- c(0.617525448736, 0.382474551264)
    exp.qmu     <- matrix(c(72669.922485424002, 0.000000000000,
                            72898.810781247870, 72898.810781247870),  nrow=2)
    exp.qsigmaf <- matrix(c(21509.940563849304, 0.000000000000,
                            21509.940563849304, 16141.114669140097),  nrow=2)
    exp.qsigmar <- matrix(c(12509.376406312418, 0.000000000000,
                            12509.376406312418, 10923.572448627263),  nrow=2)
    exp.qdelta  <- matrix(c(142.465706790323, 0.000000000000,
                            144.376399340025, 145.079185774579),  nrow=2)
    exp.qdl     <- matrix(c(3.00000000000, 0.000000000000,
                        13.000000000000, 9.000000000000),  nrow=2)
    exp.qw      <- matrix(c(0.386473566240, 0.000000000000, 1.000000000000,
                            0.382474551264),  nrow=4)

    message     <- paste0(" rjmcmc_good_result_01() ",
                      "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$mu, exp.mu, msg = message)
    checkEqualsNumeric(obs$sigmaf, exp.sigmaf, msg = message)
    checkEqualsNumeric(obs$sigmar, exp.sigmar, msg = message)
    checkEqualsNumeric(obs$delta, exp.delta, msg = message)
    checkEqualsNumeric(obs$dl, exp.dl, msg = message)
    checkEqualsNumeric(obs$w, exp.w, msg = message)
    checkEqualsNumeric(obs$qmu, exp.qmu, msg = message)
    checkEqualsNumeric(obs$qsigmaf, exp.qsigmaf, msg = message)
    checkEqualsNumeric(obs$qsigmar, exp.qsigmar, msg = message)
    checkEqualsNumeric(obs$qdelta, exp.qdelta, msg = message)
    checkEqualsNumeric(obs$qdl, exp.qdl, msg = message)
    checkEqualsNumeric(obs$qw, exp.qw, msg = message)
}

test.rjmcmc_good_result_02 <- function() {
    set.seed(101)
    obs <- RJMCMC(startPosForwardReads = reads_demo$readsForward,
                    startPosReverseReads = reads_demo$readsReverse,
                    nbrIterations = 200, lambda = 3, kmax = 30,
                    minInterval = 146, maxInterval = 292, minReads = 5)
    exp.k       <- 2
    exp.mu      <- c(72670.198935220964, 72898.810781247870)
    exp.sigmaf  <- c(21143.627354636745, 4736.790424950710)
    exp.sigmar  <- c(12296.706280648166, 5617.742233704413)
    exp.delta   <- c(144.343830717019, 143.998297780212)
    exp.dl      <- c(3, 3)
    exp.w       <- c(0.638561449326, 0.361438550674)
    exp.qmu     <- matrix(c(72451.14614368614, 0.000000000000,
                            72898.810781247870, 72898.810781247870),  nrow=2)
    exp.qsigmaf <- matrix(c(4580.202311729198, 0.000000000000,
                            21509.940563849304, 4736.790424950710),  nrow=2)
    exp.qsigmar <- matrix(c(5432.031746639804, 0.000000000000,
                            12509.376406312418, 5617.742233704413),  nrow=2)
    exp.qdelta  <- matrix(c(142.465706790323, 0.0000000000000,
                            145.843850738089, 144.943480560172),  nrow=2)
    exp.qdl     <- matrix(c(3.00000000000, 0.000000000000, 11.050000000000,
                            7.000000000000),  nrow=2)
    exp.qw      <- matrix(c(0.282811859549, 0.000000000000, 1.000000000000,
                            0.382474551264),  nrow=2)

    message     <- paste0(" rjmcmc_good_result_02() ",
                      "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$mu, exp.mu, msg = message)
    checkEqualsNumeric(obs$sigmaf, exp.sigmaf, msg = message)
    checkEqualsNumeric(obs$sigmar, exp.sigmar, msg = message)
    checkEqualsNumeric(obs$delta, exp.delta, msg = message)
    checkEqualsNumeric(obs$dl, exp.dl, msg = message)
    checkEqualsNumeric(obs$w, exp.w, msg = message)
    checkEqualsNumeric(obs$qmu, exp.qmu, msg = message)
    checkEqualsNumeric(obs$qsigmaf, exp.qsigmaf, msg = message)
    checkEqualsNumeric(obs$qsigmar, exp.qsigmar, msg = message)
    checkEqualsNumeric(obs$qdelta, exp.qdelta, msg = message)
    checkEqualsNumeric(obs$qdl, exp.qdl, msg = message)
    checkEqualsNumeric(obs$qw, exp.qw, msg = message)
}

test.rjmcmc_good_result_03 <- function() {
    set.seed(101)
    obs <- RJMCMC(startPosForwardReads = reads_demo$readsForward,
                  startPosReverseReads = reads_demo$readsReverse,
                  nbrIterations = 110, lambda = 3, kmax = 30,
                  minInterval = 100, maxInterval = 200, minReads = 335)
    exp.k       <- 2
    exp.mu      <- c(72670.440093554484, 72898.810781247870)
    exp.sigmaf  <- c(20824.077533834301, 4736.790424950710)
    exp.sigmar  <- c(12111.185532728285, 5617.742233704413)
    exp.delta   <- c(144.315419790566, 143.959205935876)
    exp.dl      <- c(4, 3)
    exp.w       <- c(0.632718493260, 0.367281506740)
    exp.qmu     <- matrix(c(72669.922485424002, 0.000000000000,
                            72898.810781247870, 72898.810781247870), nrow=2)
    exp.qsigmaf <- matrix(c(19.565623379218, 0.000000000000,
                            21509.940563849304, 4736.790424950710), nrow=2)
    exp.qsigmar <- matrix(c(32.729034009596, 0.000000000000,
                            12509.376406312418, 5617.742233704413), nrow=2)
    exp.qdelta  <- matrix(c(142.465706790323, 0.0000000000000,
                            144.376399340025, 144.043110382256), nrow=2)
    exp.qdl     <- matrix(c(3.00000000000, 0.000000000000, 13.000000000000,
                            3.000000000000),  nrow=2)
    exp.qw      <- matrix(c(0.529802297936, 0.000000000000, 1.000000000000,
                            0.470197702064),  nrow=2)

    message     <- paste0(" rjmcmc_good_result_02() ",
                          "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$mu, exp.mu, msg = message)
    checkEqualsNumeric(obs$sigmaf, exp.sigmaf, msg = message)
    checkEqualsNumeric(obs$sigmar, exp.sigmar, msg = message)
    checkEqualsNumeric(obs$delta, exp.delta, msg = message)
    checkEqualsNumeric(obs$dl, exp.dl, msg = message)
    checkEqualsNumeric(obs$w, exp.w, msg = message)
    checkEqualsNumeric(obs$qmu, exp.qmu, msg = message)
    checkEqualsNumeric(obs$qsigmaf, exp.qsigmaf, msg = message)
    checkEqualsNumeric(obs$qsigmar, exp.qsigmar, msg = message)
    checkEqualsNumeric(obs$qdelta, exp.qdelta, msg = message)
    checkEqualsNumeric(obs$qdl, exp.qdl, msg = message)
    checkEqualsNumeric(obs$qw, exp.qw, msg = message)
}
