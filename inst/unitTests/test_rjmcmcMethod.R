###################################################
# Created by Astrid Deschenes
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
    obs <- rjmcmc(startPosForwardReads = reads_demo$readsForward,
                        startPosReverseReads = reads_demo$readsReverse,
                        nbrIterations = 100, lambda = 2, kMax = 30,
                        minInterval = 146, maxInterval = 292, minReads = 5)
    exp.k       <- 2
    exp.mu      <- c(72669.922485424002, 72904.348062342819)
    exp.sigmaf  <- c(21509.940563849304, 4624.834390666467)
    exp.sigmar  <- c(12509.376406312418, 5630.799813272683)
    exp.delta   <- c(144.376399340025, 142.578131511464)
    exp.df      <- c(13, 4)
    exp.w       <- c(0.620173536972, 0.379826463028)
    exp.qmu     <- matrix(c(72669.922485424002, 72898.810781247870,
                            72669.922485424002, 72898.810781247870),  nrow = 2)
    exp.qsigmaf <- matrix(c(21509.940563849304, 4736.790424950710,
                            21509.940563849304, 4736.790424950710),  nrow = 2)
    exp.qsigmar <- matrix(c(12509.376406312418, 5617.742233704413,
                            12509.376406312418, 5617.742233704413),  nrow = 2)
    exp.qdelta  <- matrix(c(144.376399340025, 142.465706790323,
                            144.376399340025, 142.465706790323),  nrow = 2)
    exp.qdf     <- matrix(c(13.000000000000, 3.000000000000,
                            13.000000000000, 3.000000000000),  nrow = 2)
    exp.qw      <- matrix(c(0.617525448736, 0.382474551264, 0.617525448736,
                            0.382474551264),  nrow = 2)

    message     <- paste0(" rjmcmc_good_result_01() ",
                      "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$mu, exp.mu, msg = message)
    checkEqualsNumeric(obs$sigmaf, exp.sigmaf, msg = message)
    checkEqualsNumeric(obs$sigmar, exp.sigmar, msg = message)
    checkEqualsNumeric(obs$delta, exp.delta, msg = message)
    checkEqualsNumeric(obs$df, exp.df, msg = message)
    checkEqualsNumeric(obs$w, exp.w, msg = message)
    checkEqualsNumeric(obs$qmu, exp.qmu, msg = message)
    checkEqualsNumeric(obs$qsigmaf, exp.qsigmaf, msg = message)
    checkEqualsNumeric(obs$qsigmar, exp.qsigmar, msg = message)
    checkEqualsNumeric(obs$qdelta, exp.qdelta, msg = message)
    checkEqualsNumeric(obs$qdf, exp.qdf, msg = message)
    checkEqualsNumeric(obs$qw, exp.qw, msg = message)
}

test.rjmcmc_good_result_02 <- function() {
    set.seed(101)
    obs <- rjmcmc(startPosForwardReads = reads_demo$readsForward,
                    startPosReverseReads = reads_demo$readsReverse,
                    nbrIterations = 200, lambda = 3, kMax = 30,
                    minInterval = 146, maxInterval = 292, minReads = 5)
    exp.k       <- 2
    exp.mu      <- c(72669.9224854240, 72898.8107812479)
    exp.sigmaf  <- c(21509.94056384930, 4736.79042495071)
    exp.sigmar  <- c(12509.37640631242, 5617.74223370441)
    exp.delta   <- c(144.376399340025, 142.465706790323)
    exp.df      <- c(13, 3)
    exp.w       <- c(0.617525448735973, 0.382474551264027)
    exp.qmu     <- matrix(c(72669.922485424002, 72898.810781247870,
                            72669.922485424002, 72898.810781247870),  nrow = 2)
    exp.qsigmaf <- matrix(c(21509.940563849304, 4736.790424950710,
                            21509.940563849304, 4736.790424950710),  nrow = 2)
    exp.qsigmar <- matrix(c(12509.376406312418, 5617.742233704413,
                            12509.376406312418, 5617.742233704413),  nrow = 2)
    exp.qdelta  <- matrix(c(144.376399340025, 142.465706790323,
                            144.376399340025, 142.465706790323),  nrow = 2)
    exp.qdf     <- matrix(c(13.00000000000, 3.000000000000, 13.000000000000,
                            3.000000000000),  nrow = 2)
    exp.qw      <- matrix(c(0.617525448736, 0.382474551264,
                            0.617525448736, 0.382474551264),  nrow = 2)

    message     <- paste0(" rjmcmc_good_result_02() ",
                      "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$mu, exp.mu, msg = message)
    checkEqualsNumeric(obs$sigmaf, exp.sigmaf, msg = message)
    checkEqualsNumeric(obs$sigmar, exp.sigmar, msg = message)
    checkEqualsNumeric(obs$delta, exp.delta, msg = message)
    checkEqualsNumeric(obs$df, exp.df, msg = message)
    checkEqualsNumeric(obs$w, exp.w, msg = message)
    checkEqualsNumeric(obs$qmu, exp.qmu, msg = message)
    checkEqualsNumeric(obs$qsigmaf, exp.qsigmaf, msg = message)
    checkEqualsNumeric(obs$qsigmar, exp.qsigmar, msg = message)
    checkEqualsNumeric(obs$qdelta, exp.qdelta, msg = message)
    checkEqualsNumeric(obs$qdf, exp.qdf, msg = message)
    checkEqualsNumeric(obs$qw, exp.qw, msg = message)
}

test.rjmcmc_good_result_03 <- function() {
    set.seed(101)
    obs <- rjmcmc(startPosForwardReads = reads_demo$readsForward,
                  startPosReverseReads = reads_demo$readsReverse,
                  nbrIterations = 110, lambda = 3, kMax = 30,
                  minInterval = 100, maxInterval = 200, minReads = 335)
    exp.k       <- 2
    exp.mu      <- c(72669.9224854240, 72898.8107812479)
    exp.sigmaf  <- c(21509.94056384930, 4736.79042495071)
    exp.sigmar  <- c(12509.37640631242, 5617.74223370441)
    exp.delta   <- c(144.376399340025, 142.465706790323)
    exp.df      <- c(13, 3)
    exp.w       <- c(0.617525448735973, 0.382474551264027)
    exp.qmu     <- matrix(c(72669.922485424002, 72898.810781247870,
                            72669.922485424002, 72898.810781247870), nrow = 2)
    exp.qsigmaf <- matrix(c(21509.940563849304, 4736.790424950710,
                            21509.940563849304, 4736.790424950710), nrow = 2)
    exp.qsigmar <- matrix(c(12509.376406312418, 5617.742233704413,
                            12509.376406312418, 5617.742233704413), nrow = 2)
    exp.qdelta  <- matrix(c(144.376399340025, 142.465706790323,
                            144.376399340025, 142.465706790323), nrow = 2)
    exp.qdf     <- matrix(c(13.00000000000, 3.000000000000, 13.000000000000,
                            3.000000000000),  nrow = 2)
    exp.qw      <- matrix(c(0.617525448736, 0.382474551264, 0.617525448736,
                            0.382474551264),  nrow = 2)

    message     <- paste0(" rjmcmc_good_result_03() ",
                          "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$mu, exp.mu, msg = message)
    checkEqualsNumeric(obs$sigmaf, exp.sigmaf, msg = message)
    checkEqualsNumeric(obs$sigmar, exp.sigmar, msg = message)
    checkEqualsNumeric(obs$delta, exp.delta, msg = message)
    checkEqualsNumeric(obs$df, exp.df, msg = message)
    checkEqualsNumeric(obs$w, exp.w, msg = message)
    checkEqualsNumeric(obs$qmu, exp.qmu, msg = message)
    checkEqualsNumeric(obs$qsigmaf, exp.qsigmaf, msg = message)
    checkEqualsNumeric(obs$qsigmar, exp.qsigmar, msg = message)
    checkEqualsNumeric(obs$qdelta, exp.qdelta, msg = message)
    checkEqualsNumeric(obs$qdf, exp.qdf, msg = message)
    checkEqualsNumeric(obs$qw, exp.qw, msg = message)
}

test.rjmcmc_good_result_04 <- function() {
    set.seed(331)
    obs <- rjmcmc(startPosForwardReads = reads_demo_02$readsForward,
                  startPosReverseReads = reads_demo_02$readsReverse,
                  nbrIterations = 210, lambda = 3, kMax = 30,
                  minInterval = 100, maxInterval = 200, minReads = 10)
    exp.k       <- 6
    exp.mu      <- c(18830.258426897533, 19299.722102507523, 19381.662157933497,
                     19438.522564048224, 19474.666149448600, 19661.608468284954)
    exp.sigmaf  <- c(19360.298739977090, 1092.324045407637, 674.355821083463,
                     384.223266693146, 295.997780079298, 396.751032236276)
    exp.sigmar  <- c(8246.736691816854, 999.816627725857, 331.571444139750,
                     79.62427776446, 306.565419274503, 233.932761372705)
    exp.delta   <- c(151.637089900930, 145.218532805779, 146.316731025536,
                     151.310459710703, 145.543765409165, 143.900213132404)
    exp.df      <- c(3, 3, 3, 5, 11, 3)
    exp.w       <- c(0.311173734828, 0.221181496934, 0.175512530366,
                     0.050440609347, 0.181847532684, 0.059844095841)
    exp.qmu     <- matrix(c(18830.258426897533, 19299.722102507523, 19381.662157933497,
                            19438.522564048224, 19451.890925142368, 19661.608468284954,
                            18830.258426897533, 19299.722102507523, 19381.662157933497,
                            19438.522564048224, 19551.751524023533, 19661.608468284954),
                          nrow = 2)
    exp.qsigmaf <- matrix(c(19360.298739977090, 1092.324045407637,
                            674.355821083463,   384.223266693146,
                            159.149766899767,   396.751032236276,
                            19360.298739977090, 1092.324045407637,
                            674.355821083463,   384.223266693146,
                            759.175670840787,   396.751032236276),
                          nrow = 2)
    exp.qsigmar <- matrix(c(8246.736691816854,  999.816627725857,
                            331.571444139750,   79.624277764465,
                            218.518977049937,   233.932761372705,
                            8246.736691816854,  999.816627725857,
                            331.571444139750,   79.624277764465,
                            604.568762188420,   233.932761372705),
                          nrow = 2)
    exp.qdelta  <- matrix(c(151.637089900930, 145.218532805779,
                            146.316731025536, 151.310459710703,
                            143.104948985660, 143.900213132404,
                            151.637089900930, 145.218532805779,
                            146.316731025536, 151.310459710703,
                            146.264324807019, 143.900213132404),
                          nrow = 2)
    exp.qdf     <- matrix(c(3.000000000000, 3.000000000000, 3.000000000000,
                            3.000000000000, 3.000000000000, 3.000000000000,
                            3.000000000000, 3.000000000000, 3.000000000000,
                            13.000000000000, 13.000000000000, 3.000000000000),
                          nrow = 2)
    exp.qw      <- matrix(c(0.176647544632, 0.207703263307,
                            0.171681803342, 0.016202988344,
                            0.047681376154, 0.032004302325,
                            0.350920109204, 0.266800133824,
                            0.188478067985, 0.166321788126,
                            0.221487533477, 0.154071089279), nrow = 2)

    message     <- paste0(" test.rjmcmc_good_result_04() ",
                          "- RJMCMC did not generated expected values")

    checkEqualsNumeric(obs$k, exp.k, msg = message)
    checkEqualsNumeric(obs$mu, exp.mu, msg = message)
    checkEqualsNumeric(obs$sigmaf, exp.sigmaf, msg = message)
    checkEqualsNumeric(obs$sigmar, exp.sigmar, msg = message)
    checkEqualsNumeric(obs$delta, exp.delta, msg = message)
    checkEqualsNumeric(obs$df, exp.df, msg = message)
    checkEqualsNumeric(obs$w, exp.w, msg = message)
    checkEqualsNumeric(obs$qmu, exp.qmu, msg = message)
    checkEqualsNumeric(obs$qsigmaf, exp.qsigmaf, msg = message)
    checkEqualsNumeric(obs$qsigmar, exp.qsigmar, msg = message)
    checkEqualsNumeric(obs$qdelta, exp.qdelta, msg = message)
    checkEqualsNumeric(obs$qdf, exp.qdf, msg = message)
    checkEqualsNumeric(obs$qw, exp.qw, msg = message)
}

