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
data(reads_demo_02)


###########################################################
## RJMCMC() function
###########################################################

test.rjmcmc_good_result_01 <- function() {
    set.seed(101)
    obs <- RJMCMC(startPosForwardReads = reads_demo$readsForward,
                        startPosReverseReads = reads_demo$readsReverse,
                        nbrIterations = 100, lambda = 2, kMax = 30,
                        minInterval = 146, maxInterval = 292, minReads = 5)
    exp.k       <- 2
    exp.mu      <- c(72669.922485424002, 72898.810781247870)
    exp.sigmaf  <- c(21509.940563849304, 4736.790424950710)
    exp.sigmar  <- c(12509.376406312418, 5617.742233704413)
    exp.delta   <- c(144.376399340025, 142.465706790323)
    exp.dl      <- c(13, 3)
    exp.w       <- c(0.617525448736, 0.382474551264)
    exp.qmu     <- matrix(c(72669.9224854240, 72883.3070521507,
                            72898.8107812479, 72898.8107812479),  nrow = 2)
    exp.qsigmaf <- matrix(c(21509.94056384930, 4736.79042495071,
                            21509.94056384930, 4736.79042495071),  nrow = 2)
    exp.qsigmar <- matrix(c(12509.37640631242, 5617.74223370441,
                            12509.37640631242, 5617.74223370441),  nrow = 2)
    exp.qdelta  <- matrix(c(142.465706790323, 142.465706790323,
                            144.376399340025, 142.465706790323),  nrow = 2)
    exp.qdl     <- matrix(c(6.80000000000, 3.000000000000,
                            13.000000000000, 3.000000000000),  nrow = 2)
    exp.qw      <- matrix(c(0.617525448735973, 0.331484238246196, 1.000000000000,
                            0.382474551264027),  nrow = 2)

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
                    nbrIterations = 200, lambda = 3, kMax = 30,
                    minInterval = 146, maxInterval = 292, minReads = 5)
    exp.k       <- 2
    exp.mu      <- c(72669.9224854240, 72898.8107812479)
    exp.sigmaf  <- c(21509.94056384930, 4736.79042495071)
    exp.sigmar  <- c(12509.37640631242, 5617.74223370441)
    exp.delta   <- c(144.376399340025, 142.465706790323)
    exp.dl      <- c(13, 3)
    exp.w       <- c(0.617525448735973, 0.382474551264027)
    exp.qmu     <- matrix(c(72460.7128374504, 72669.9224854240,
                            72898.8107812479, 72898.8107812479),  nrow = 2)
    exp.qsigmaf <- matrix(c(21509.94056384930, 4736.79042495071,
                            21509.94056384930, 4736.79042495071),  nrow = 2)
    exp.qsigmar <- matrix(c(12509.37640631242, 5617.74223370441,
                            12509.37640631242, 5617.74223370441),  nrow = 2)
    exp.qdelta  <- matrix(c(142.465706790323, 142.465706790323,
                            145.861078485865, 144.376399340025),  nrow = 2)
    exp.qdl     <- matrix(c(11.00000000000, 3.000000000000, 19.000000000000,
                            3.000000000000),  nrow = 2)
    exp.qw      <- matrix(c(0.266735697031826, 0.253034010399643,
                            1.000000000000000,
                            0.382474551264027),  nrow = 2)

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
                  nbrIterations = 110, lambda = 3, kMax = 30,
                  minInterval = 100, maxInterval = 200, minReads = 335)
    exp.k       <- 2
    exp.mu      <- c(72669.9224854240, 72898.8107812479)
    exp.sigmaf  <- c(21509.94056384930, 4736.79042495071)
    exp.sigmar  <- c(12509.37640631242, 5617.74223370441)
    exp.delta   <- c(144.376399340025, 142.465706790323)
    exp.dl      <- c(13, 3)
    exp.w       <- c(0.617525448735973, 0.382474551264027)
    exp.qmu     <- matrix(c(72669.9224854240, 72898.8107812479,
                            72898.8107812479, 72898.8107812479), nrow = 2)
    exp.qsigmaf <- matrix(c(21509.94056384930, 4736.79042495071,
                            21509.94056384930, 4736.79042495071), nrow = 2)
    exp.qsigmar <- matrix(c(12509.37640631242, 5617.74223370441,
                            12509.37640631242, 5617.74223370441), nrow = 2)
    exp.qdelta  <- matrix(c(142.465706790323, 142.465706790323,
                            144.376399340025, 142.465706790323), nrow = 2)
    exp.qdl     <- matrix(c(11.00000000000, 3.000000000000, 13.000000000000,
                            3.000000000000),  nrow = 2)
    exp.qw      <- matrix(c(0.617525448735973, 0.382474551264027, 1.000000000000,
                            0.382474551264027),  nrow = 2)

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

test.rjmcmc_good_result_04 <- function() {
    set.seed(331)
    obs <- RJMCMC(startPosForwardReads = reads_demo_02$readsForward,
                  startPosReverseReads = reads_demo_02$readsReverse,
                  nbrIterations = 210, lambda = 3, kMax = 30,
                  minInterval = 100, maxInterval = 200, minReads = 10)
    exp.k       <- 7
    exp.mu      <- c(18830.2584268975, 18994.0909090909, 19130.9529411765,
                     19299.7221025075, 19410.9770121171, 19524.9898072680,
                     19661.6084682850)
    exp.sigmaf  <- c(19360.298739977090, 10226.311392692363, 5280.179285721085,
                     448.243057393201, 396.751032236276, 396.751032236276,
                     396.751032236276)
    exp.sigmar  <- c(8246.736691816854, 4623.276659771355, 2442.451348288455,
                     286.591428155266, 233.932761372705, 233.932761372705,
                     233.932761372705)
    exp.delta   <- c(151.637089900930, 148.427811353354, 147.120483906917,
                     143.613396554890, 143.726384903608, 143.726384903608,
                     143.726384903608)
    exp.dl      <- c(3, 3, 7, 3, 3, 3, 3)
    exp.w       <- c(0.2808308442231622, 0.2452808306768512, 0.2049774378308375,
                     0.0569431136923228, 0.0706559245256088, 0.0706559245256088,
                     0.0706559245256088)
    exp.qmu     <- matrix(c(18830.2584268975, 18904.7226562500, 18995.6777777778,
                            19121.6120802148, 19288.3044692737, 19387.5833333333,
                            19551.7515240235, 18830.2584268975, 18994.0909090909,
                            19130.9529411765, 19299.7221025075, 19430.0992779783,
                            19551.7515240235, 19661.6084682850), nrow = 2)
    exp.qsigmaf <- matrix(c(19360.298739977090, 10226.311392692363,
                            5192.730579796065, 396.751032236276,
                            396.751032236276, 396.751032236276,
                            396.751032236276, 19360.298739977090,
                            14688.813010253683, 7077.307093476294,
                            1319.680196024897, 396.751032236276,
                            396.751032236276, 396.751032236276), nrow = 2)
    exp.qsigmar <- matrix(c(8246.736691816854, 4623.276659771355,
                            2351.450468767910, 233.932761372705,
                            233.932761372705, 233.932761372705,
                            233.932761372705, 8246.736691816854,
                            8198.509840053977, 6156.291229289655,
                            2174.002689949019, 233.932761372705,
                            233.932761372705, 233.932761372705), nrow = 2)
    exp.qdelta  <- matrix(c(151.637089900930, 148.427811353354,
                            145.766380169507, 143.104948985660,
                            143.104948985660, 143.104948985660,
                            143.104948985660, 151.637089900930,
                            150.307000182081, 149.869135532029,
                            144.556807719391, 143.900213132404,
                            143.900213132404, 143.900213132404), nrow = 2)
    exp.qdl     <- matrix(c(3.000000000000, 3.000000000000, 3.000000000000,
                            3.000000000000, 3.000000000000, 3.000000000000,
                            3.000000000000, 3.000000000000, 9.000000000000,
                            9.000000000000, 3.000000000000, 3.000000000000,
                            3.000000000000, 3.000000000000), nrow = 2)
    exp.qw      <- matrix(c(0.1602552234736732, 0.1726393303598764,
                            0.1519735262157934, 0.0033610588369209,
                            0.0033610588369209, 0.0033610588369209,
                            0.0033610588369209, 0.3479111596593648,
                            0.3939779874925808, 0.2623494692460391,
                            0.1094172630972199, 0.1397737902033675,
                            0.1397737902033675, 0.1397737902033675), nrow = 2)

    message     <- paste0(" test.rjmcmc_good_result_04() ",
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

