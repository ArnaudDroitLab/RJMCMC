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
    exp.mu      <- c(72669.922485424002, 72898.810781247870)
    exp.sigmaf  <- c(21509.940563849304, 4736.790424950710)
    exp.sigmar  <- c(12509.376406312418, 5617.742233704413)
    exp.delta   <- c(144.376399340025, 142.465706790323)
    exp.df      <- c(13, 3)
    exp.w       <- c(0.617525448736, 0.382474551264)
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
    exp.k       <- 4
    exp.mu      <- c(18830.2584268975, 18994.0909090909, 19130.9529411765,
                     19299.7221025075, 19410.9770121171, 19524.9898072680,
                     19661.6084682850)
    exp.sigmaf  <- c(19360.298739977090, 1141.230922288118, 347.461377747992,
                     396.751032236276)
    exp.sigmar  <- c(8246.736691816854, 1128.665134955233, 337.427609866408,
                     233.932761372705)
    exp.delta   <- c(151.637089900930, 145.390326636372, 145.196502385400,
                     143.865885183624)
    exp.df      <- c(3, 6, 11, 3)
    exp.w       <- c(0.378890356569, 0.336939424563, 0.221839702904,
                     0.062330515963)
    exp.qmu     <- matrix(c(18830.258426897533, 19201.761590246519, 19403.403846153848,
                            19661.608468284954, 18830.258426897533, 19337.633333333335,
                            19553.238141579601, 19661.608468284954), nrow = 2)
    exp.qsigmaf <- matrix(c(19360.298739977090, 674.355821083463,
                            21.507805325987,   396.751032236276,
                            19360.298739977090,    2946.586823298956,
                            759.175670840787,    396.751032236276), nrow = 2)
    exp.qsigmar <- matrix(c(8246.736691816854, 331.571444139750,
                            77.862248759566, 233.932761372705,
                            8246.736691816854,  4932.573126036619,
                            604.568762188420,  233.932761372705), nrow = 2)
    exp.qdelta  <- matrix(c(151.637089900930, 145.218532805779,
                            143.104948985660, 143.104948985660,
                            151.637089900930, 146.580095508417,
                            146.264324807019, 143.900213132404), nrow = 2)
    exp.qdf     <- matrix(c(3.000000000000, 3.000000000000, 3.000000000000,
                            3.000000000000, 3.000000000000, 27.000000000000,
                            17.950000000000, 3.000000000000), nrow = 2)
    exp.qw      <- matrix(c(0.242004820526, 0.255755908990,
                            0.432106313995, 0.601840945406,
                            0.272729202963, 0.160678776966), nrow = 2)

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

