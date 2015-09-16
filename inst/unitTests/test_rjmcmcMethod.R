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
    exp.qsigmaf <- matrix(c(21509.940563849305, 0.000000000000,
                            21509.940563849305, 4736.79042495071),  nrow=2)
    exp.qsigmar <- matrix(c(12509.376406312418, 0.000000000000,
                            12509.376406312418, 5617.742233704413),  nrow=2)
    exp.qdelta  <- matrix(c(142.4657067903234, 0.0000000000000,
                            144.3763993400250, 142.4657067903234),  nrow=2)
    exp.qdl     <- matrix(c(6.80000000000, 0.000000000000,
                        13.000000000000, 3.000000000000),  nrow=2)
    exp.qw      <- matrix(c(0.6175254487359728, 0.000000000000, 1.000000000000,
                            0.3824745512640272),  nrow=2)

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
    exp.mu      <- c(72669.9224854240, 72898.8107812479)
    exp.sigmaf  <- c(21509.94056384930, 4736.79042495071)
    exp.sigmar  <- c(12509.37640631242, 5617.74223370441)
    exp.delta   <- c(144.376399340025, 142.465706790323)
    exp.dl      <- c(13, 3)
    exp.w       <- c(0.617525448735973, 0.382474551264027)
    exp.qmu     <- matrix(c(72460.7128374504, 0.000000000000,
                            72898.8107812479, 72898.8107812479),  nrow=2)
    exp.qsigmaf <- matrix(c(21509.94056384930, 0.000000000000,
                            21509.94056384930, 4736.79042495071),  nrow=2)
    exp.qsigmar <- matrix(c(12509.37640631242, 0.000000000000,
                            12509.37640631242, 5617.74223370441),  nrow=2)
    exp.qdelta  <- matrix(c(142.465706790323, 0.0000000000000,
                            145.861078485865, 144.376399340025),  nrow=2)
    exp.qdl     <- matrix(c(11.00000000000, 0.000000000000, 19.000000000000,
                            3.000000000000),  nrow=2)
    exp.qw      <- matrix(c(0.2667356970318264, 0.0000000000000000,
                            1.0000000000000000,
                            0.3824745512640272),  nrow=2)

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
    exp.mu      <- c(72669.9224854240, 72898.8107812479)
    exp.sigmaf  <- c(21509.94056384930, 4736.79042495071)
    exp.sigmar  <- c(12509.37640631242, 5617.74223370441)
    exp.delta   <- c(144.376399340025, 142.465706790323)
    exp.dl      <- c(13, 3)
    exp.w       <- c(0.617525448735973, 0.382474551264027)
    exp.qmu     <- matrix(c(72669.9224854240, 0.000000000000,
                            72898.8107812479, 72898.8107812479), nrow=2)
    exp.qsigmaf <- matrix(c(21509.94056384930, 0.000000000000,
                            21509.94056384930, 4736.79042495071), nrow=2)
    exp.qsigmar <- matrix(c(12509.37640631242, 0.000000000000,
                            12509.37640631242, 5617.74223370441), nrow=2)
    exp.qdelta  <- matrix(c(142.465706790323, 0.0000000000000,
                            144.376399340025, 142.465706790323), nrow=2)
    exp.qdl     <- matrix(c(11.00000000000, 0.000000000000, 13.000000000000,
                            3.000000000000),  nrow=2)
    exp.qw      <- matrix(c(0.617525448735973, 0.000000000000, 1.000000000000,
                            0.382474551264027),  nrow=2)

    message     <- paste0(" rjmcmc_good_result_03() ",
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
