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
    exp.k       <- 3
    exp.mu      <- c(72583.158931887257, 72682.48124999999, 72898.810781247870)
    exp.sigmaf  <- c(21509.940563849304, 13217.346965399556, 4924.753366949810)
    exp.sigmar  <- c(12509.376406312418, 8458.825226487852, 4408.274046663286)
    exp.delta   <- c(144.555002556739, 144.84904296083, 145.143083364930)
    exp.dl      <- c(4, 4, 3)
    exp.w       <- c(0.47755744993, 0.333333333333, 0.189109216735)
    exp.qmu     <- matrix(c(72583.158931887257, 0.000000000000, 0.000000000000,
                            72898.810781247870, 72736.151927083338,
                            72898.810781247870),  nrow=3)
    exp.qsigmaf <- matrix(c(21509.940563849304, 0.000000000000, 0.000000000000,
                            21509.940563849304, 13217.346965399556,
                            4924.753366949810),  nrow=3)
    exp.qsigmar <- matrix(c(12509.376406312418, 0.000000000000, 0.000000000000,
                            12509.376406312418, 8458.825226487852,
                            4408.274046663286),  nrow=3)
    exp.qdelta  <- matrix(c(142.465706790323, 0.000000000000, 0.000000000000,
                           144.555002556739, 145.258601627752,
                           145.477755436756),  nrow=3)
    exp.qdl     <- matrix(c(3.00000000000, 0.000000000000, 0.000000000000,
                        15.000000000000, 9.000000000000,
                        3.000000000000),  nrow=3)
    exp.qw      <- matrix(c(0.46612915689, 0.000000000000, 0.000000000000,
                        1.000000000000, 0.338513883737,
                        0.243157370856),  nrow=3)

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
    exp.k       <- 3
    exp.mu      <- c(72583.15893188726, 72772.61459313818, 73014.46314048163)
    exp.sigmaf  <- c(21509.94056384930, 10956.01225808368, 402.08395231807)
    exp.sigmar  <- c(12509.37640631242, 6426.88752640830, 344.39864650418)
    exp.delta   <- c(144.55500255674, 145.67654791164, 146.79809326653)
    exp.dl      <- c(3, 15, 26)
    exp.w       <- c(0.61251561781, 0.333333333333, 0.05415104886)
    exp.qmu     <- matrix(c(72583.15893188726, 0.000000000000, 0.000000000000,
                        72898.81078124787, 72783.21695627819,
                        73024.95025541677),  nrow=3)
    exp.qsigmaf <- matrix(c(60.95168449911, 0.000000000000, 0.000000000000,
                            21509.94056384930, 13217.34696539956,
                            4924.75336694981),  nrow=3)
    exp.qsigmar <- matrix(c(79.66429332594, 0.000000000000, 0.000000000000,
                            12509.37640631242, 8893.79951911038,
                            4408.27404666329),  nrow=3)
    exp.qdelta  <- matrix(c(142.465706790323, 0.0000000000000, 0.000000000000,
                           144.555002556739, 145.779156053732,
                           147.790003631531),  nrow=3)
    exp.qdl     <- matrix(c(3.00000000000, 0.000000000000, 0.000000000000,
                        15.025000000000, 16.000000000000,
                        29.000000000000),  nrow=3)
    exp.qw      <- matrix(c(0.09371828799, 0.000000000000, 0.000000000000,
                       1.00000000000, 0.49987309733,
                       0.57294837867),  nrow=3)

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
