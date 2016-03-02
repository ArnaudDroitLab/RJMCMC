###################################################
# Created by Astrid Deschenes
# 2015-06-12
###################################################

###################################################
## Test the rjmcmcMethodsIntern.R functions
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "RJMCMC" )
}

### }}}


file_002 <- dir(system.file("extdata", package = "RJMCMC"),
                pattern = "yeastRes_Chr1_Seg_002.rds",
                full.names = TRUE)

data_002 <- readRDS(file_002)

data(RJMCMC_result)
data(reads_demo)

###########################################################
## Dk() function
###########################################################

test.Dk_result_k_inferior_to_kmax <- function() {
    obs <- mapply(c(8, 6, 5, 4, 3, 1, 0, -1),
                    FUN = function(x) {RJMCMC:::Dk(9, 8, x)})
    message <- paste0(" Dk_result_k_inferior_to_kmax() ",
                    "- A kmax inferior to k did not always returned 0")
    checkIdentical(obs, rep(0, 8), msg = message)
}

test.Dk_result_with_various_values_of_k <- function() {
    obs <- mapply(c(9, 8, 6, 5, 4, 3, 1, 0),
                    FUN = function(x) {RJMCMC:::Dk(x, 8, 10)})
    exp <- c(0.5000, 0.5000, 0.3750, 0.3125, 0.2500, 0.1875,
                    0.0000, 0.0000)
    message <- paste0(" Dk_result_with_various_values_of_k() ",
                    "- Not all tested data generated expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

test.Dk_result_with_various_values_of_lambda <- function() {
    obs <- mapply(c(9, 8, 6, 35, 24, 13, 1),
                    FUN = function(x) {RJMCMC:::Dk(7, x, 10)})
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
                  FUN = function(x) {RJMCMC:::Bk(9, 8, x)})
    message <- paste0(" Bk_result_k_inferior_to_kmax() ",
                    "- A kmax inferior to k did not always returned 0")
    checkIdentical(obs, rep(0, 8), msg = message)
}

test.Bk_result_with_various_values_of_k <- function() {
    obs <- mapply(c(9, 8, 6, 5, 4, 3, 1, 0),
                    FUN = function(x) {RJMCMC:::Bk(x, 8, 10)})
    exp <- c(0.4000000, 0.444444444444444, 0.5000000, 0.5000000, 0.5000000,
                    0.5000000, 0.5000000, 0.5000000)
    message <- paste0(" Bk_result_with_various_values_of_k() ",
                    "- Not all tested data generated expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}

test.Bk_result_with_various_values_of_lambda <- function() {
    obs <- mapply(c(9, 8, 6, 2, 24, 3, 1),
                  FUN = function(x) {RJMCMC:::Bk(7, x, 10)})
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
                  FUN = function(x) {RJMCMC:::tnormale(x, 7, x-1, x+1)})
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
                  FUN = function(x) {RJMCMC:::tnormale(8, 10090, x, 10)})
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
                  FUN = function(x) {RJMCMC:::tnormale(10, 10090, 10, x)})
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
                  FUN = function(x) {RJMCMC:::tnormale(8, x, 7, 10)})
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
                  FUN = function(x) {RJMCMC:::elementWithHighestMode(x)})
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
    readsForward <- sapply(1:nbrReads, RJMCMC:::normal.mixture, k = k,
                        w = weight, mu = mu - delta/2, sigma = sigma)
    readsReverse <- sapply(1:nbrReads, RJMCMC:::normal.mixture, k = k,
                        w = weight, mu = mu - delta/2, sigma = sigma)
    reads <- sort(c(readsForward, readsReverse))
    obs <- RJMCMC:::priorMuDensity(mu, reads)
    exp <- 8.0883349761e-15
    message <- paste0(" priorMuDensity_with_1000_reads() ",
                    "- The result is not the expected value.")
    checkEqualsNumeric(obs, exp, msg = message)
}

test.priorMuDensity_results_with_various_values_of_mu <- function() {
    set.seed(101)
    obs <- mapply(list(A=c(10200, 10300), B=c(10108, 10206, 10222),
                    C=c(10333, 10455, 10899)),
                    FUN = function(x) { RJMCMC:::priorMuDensity(x,
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
                    FUN = function(x) { RJMCMC:::normal.mixture(k = k,
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
                  FUN = function(x) { RJMCMC:::student.mixture(k = k,
                                        weight = weight, mu = x,
                                        sigma = sigma, dfr = dfr)} )
    exp <- c(3.61671572308, 17.43930668620, 201.78192965619)
    message <- paste0(" student_mixture_results_with_various_values_of_mu() ",
                      "- Not all tested data with various mu generated",
                      " expected values.")
    checkEqualsNumeric(obs, exp, msg = message)
}


#########################################################
## validatePrepMergeParameters() function
#########################################################

## Test the result when startPosForwardReads is NA
test.validatePrepMergeParameters_startPosForwardReads_NA <- function() {
    obs <- tryCatch(RJMCMC:::validatePrepMergeParameters(
                        startPosForwardReads = NA,
                        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                        72429.24, 72426.08),
                        resultRJMCMC = data_002, extendingSize = 11,
                        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("startPosForwardReads must be a non-empty vector of ",
                    "numeric values.")
    message <- paste0(" test.validatePrepMergeParameters_startPosForwardReads_NA() ",
                      "- NA for startPosForwardReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosForwardReads is empty
test.validatePrepMergeParameters_startPosForwardReads_empty <- function() {
    obs <- tryCatch(RJMCMC:::validatePrepMergeParameters(
        startPosForwardReads = c(),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("startPosForwardReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validatePrepMergeParameters_startPosForwardReads_empty() ",
                      "- empty startPosForwardReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosForwardReads is not number
test.validatePrepMergeParameters_startPosForwardReads_not_number <- function() {
    obs <- tryCatch(RJMCMC:::validatePrepMergeParameters(
        startPosForwardReads = c("A", "B"),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("startPosForwardReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validatePrepMergeParameters_startPosForwardReads_not_number() ",
                      "- not number startPosForwardReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosReverseReads is NA
test.validatePrepMergeParameters_startPosReverseReads_NA <- function() {
    obs <- tryCatch(RJMCMC:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = NA,
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("startPosReverseReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validatePrepMergeParameters_startPosReverseReads_NA() ",
                      "- NA for startPosReverseReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosReverseReads is empty
test.validatePrepMergeParameters_startPosReverseReads_empty <- function() {
    seqinfo <- GenomeInfoDb::Seqinfo(c("chr1"), c(1000000), NA, "mock1")
    obs <- tryCatch(RJMCMC:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(),
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("startPosReverseReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validatePrepMergeParameters_startPosReverseReads_empty() ",
                      "- empty startPosReverseReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosReverseReads is not number
test.validatePrepMergeParameters_startPosReverseReads_not_number <- function() {
    obs <- tryCatch(RJMCMC:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c("A", "B"),
        resultRJMCMC = data_002, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("startPosReverseReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validatePrepMergeParameters_startPosReverseReads_not_number() ",
                      "- not number startPosReverseReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when resultRJMCMC is NA
test.validatePrepMergeParameters_resultRJMCMC_NA <- function() {
    obs <- tryCatch(RJMCMC:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = NA, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("resultRJMCMC must be an object of class",
                  "\'rjmcmcNucleosomes\' or \'rjmcmcNucleosomesMerge\'.")
    message <- paste0(" test.validatePrepMergeParameters_startPosReverseReads_not_number() ",
                      "- NA resultRJMCMC did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when resultRJMCMC is a number
test.validatePrepMergeParameters_resultRJMCMC_number <- function() {
    obs <- tryCatch(RJMCMC:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = 33, extendingSize = 11,
        chrLength = 1000000), error=conditionMessage)
    exp <- paste0("resultRJMCMC must be an object of class",
                  "\'rjmcmcNucleosomes\' or \'rjmcmcNucleosomesMerge\'.")
    message <- paste0(" test.validatePrepMergeParameters_resultRJMCMC_number() ",
                      "- number resultRJMCMC did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nbBase is a string
test.validatePrepMergeParameters_nbBase_string <- function() {
    obs <- tryCatch(RJMCMC:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = "ALLO",
        chrLength = 1000000), error=conditionMessage)
    exp <- "extendingSize must be a positive integer or numeric"
    message <- paste0(" test.validatePrepMergeParameters_nbBase_number() ",
                      "- string nbBase did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nbBase is an array
test.validatePrepMergeParameters_nbBase_array <- function() {
    obs <- tryCatch(RJMCMC:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = c(10, 11),
        chrLength = 1000000), error=conditionMessage)
    exp <- "extendingSize must be a positive integer or numeric"
    message <- paste0(" test.validatePrepMergeParameters_nbBase_string() ",
                      "- array nbBase did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrLength is a string
test.validatePrepMergeParameters_chrLength_string <- function() {
    obs <- tryCatch(RJMCMC:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = 74,
        chrLength = "5000"), error=conditionMessage)
    exp <- "chrLength must be a positive integer or numeric"
    message <- paste0(" test.validatePrepMergeParameters_chrLength_string() ",
                      "- string chrLength did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when chrLength is an array
test.validatePrepMergeParameters_chrLength_array <- function() {
    obs <- tryCatch(RJMCMC:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = 74,
        chrLength = c(100, 200)), error=conditionMessage)
    exp <- "chrLength must be a positive integer or numeric"
    message <- paste0(" test.validatePrepMergeParameters_chrLength_string() ",
                      "- array chrLength did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when all parameters are valid
test.validatePrepMergeParameters_all_valid <- function() {
    obs <- RJMCMC:::validatePrepMergeParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        resultRJMCMC = data_002, extendingSize = 74,
        chrLength = 200000)
    exp <- 0
    message <- paste0(" test.validatePrepMergeParameters_all_valid() ",
                      "- All valid parameters did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}


#########################################################
## validateParameters() function
#########################################################

## Test the result when nbrIterations is NA
test.validateParameters_nbrIterations_NA <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                    72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                    72429.24, 72426.08),
        nbrIterations = NA,
        kMax = 4, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "nbrIterations must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_nbrIterations_NA() ",
                      "- NA for nbrIterations did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nbrIterations is zero
test.validateParameters_nbrIterations_zero <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 0,
        kMax = 4, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "nbrIterations must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_nbrIterations_zero() ",
                      "- Zero for nbrIterations did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when nbrIterations is negative
test.validateParameters_nbrIterations_negative <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = -1,
        kMax = 4, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "nbrIterations must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_nbrIterations_zero() ",
                      "- Negative value for nbrIterations did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when kMax is NA
test.validateParameters_kMax_NA <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = NA, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "kMax must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_kMax_NA() ",
                      "- NA value for kMax did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when kMax is zero
test.validateParameters_kMax_zero <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters (
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 0, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "kMax must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_kMax_zero() ",
                      "- Zero value for kMax did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when kMax is negative
test.validateParameters_kMax_negative <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters (
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = -1, lambda = 1, minReads = 5, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "kMax must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_kMax_negative() ",
                      "- Negative value for kMax did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when minReads is NA
test.validateParameters_minReads_NA <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = NA, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "minReads must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_minReads_NA() ",
                      "- NA value for minReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when minReads is zero
test.validateParameters_minReads_zero <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 0, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "minReads must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_minReads_zero() ",
                      "- Zero value for minReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when minReads is negative
test.validateParameters_minReads_negative <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 3, lambda = 1, minReads = -1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "minReads must be a positive integer or numeric"
    message <- paste0(" test.validateParameters_minReads_negative() ",
                      "- Negative value for minReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when lambda is NA
test.validateParameters_lambda_NA <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = NA, minReads = 2, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "lambda must be a positive numeric"
    message <- paste0(" test.validateParameters_minReads_NA() ",
                      "- NA value for lambda did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when lambda is zero
test.validateParameters_lambda_zero <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 0, minReads = 3, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "lambda must be a positive numeric"
    message <- paste0(" test.validateParameters_minReads_zero() ",
                      "- Zero value for lambda did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when lambda is negative
test.validateParameters_lambda_negative <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = -1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- "lambda must be a positive numeric"
    message <- paste0(" test.validateParameters_minReads_negative() ",
                      "- Negative value for lambda did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosForwardReads is NA
test.validateParameters_startPosForwardReads_NA <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = NA,
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- paste0("startPosForwardReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validateParameters_startPosForwardReads_NA() ",
                      "- NA value for startPosForwardReads did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosForwardReads is empty array
test.validateParameters_startPosForwardReads_empty <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- paste0("startPosForwardReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validateParameters_startPosForwardReads_empty() ",
                        "- Empty array for startPosForwardReads did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosReverseReads is NA
test.validateParameters_startPosReverseReads_NA <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = NA,
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- paste0("startPosReverseReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validateParameters_startPosReverseReads_NA() ",
                        "- NA value for startPosReverseReads did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when startPosReverseReads is empty array
test.validateParameters_startPosReverseReads_empty_array <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = FALSE), error=conditionMessage)
    exp <- paste0("startPosReverseReads must be a non-empty vector of ",
                  "numeric values.")
    message <- paste0(" test.validateParameters_startPosReverseReads_empty_array() ",
                        "- Empty array for startPosReverseReads did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when adaptIterationsToReads is string
test.validateParameters_adaptIterationsToReads_string <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = "allo"), error=conditionMessage)
    exp <- "adaptIterationsToReads must be a logical."
    message <- paste0(" test.validateParameters_adaptIterationsToReads_string() ",
                        "- String for adaptIterationsToReads did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when adaptIterationsToReads is number
test.validateParameters_adaptIterationsToReads_number <- function() {
    obs <- tryCatch(RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = 33), error=conditionMessage)
    exp <- "adaptIterationsToReads must be a logical."
    message <- paste0(" test.validateParameters_adaptIterationsToReads_number() ",
                        "- Number value for adaptIterationsToReads did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when all parameters are valid
test.validateParameters_all_valid <- function() {
    obs <- RJMCMC:::validateParameters(
        startPosForwardReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        startPosReverseReads = c(72424.14, 72431.49, 72428.21,
                                 72429.24, 72426.08),
        nbrIterations = 2,
        kMax = 10, lambda = 1, minReads = 1, minInterval = 146,
        maxInterval = 292,
        adaptIterationsToReads = TRUE)
    exp <- 0
    message <- paste0(" test.validateParameters_all_valid() ",
                      "- All valid parameters did not  ",
                      "generated expected message.")
    checkEquals(obs, exp, msg = message)
}


#########################################################
## validateRDSFilesParameters() function
#########################################################

## Test the result when RDSFiles is NA
test.validateRDSFilesParameters_RDSFiles_NA <- function() {
    obs <- tryCatch(RJMCMC:::validateRDSFilesParameters(
        RDSFiles = NA), error=conditionMessage)
    exp <- "RDSFiles must be a list of valid RDS files"
    message <- paste0(" test.validateRDSFilesParameters_RDSFiles_NA() ",
                        "- NA for RDSFiles did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}

## Test the result when RDSFiles is empty array
test.validateRDSFilesParameters_RDSFiles_empty_array <- function() {
    obs <- tryCatch(RJMCMC:::validateRDSFilesParameters(
        RDSFiles = c()), error=conditionMessage)
    exp <- "RDSFiles must be a list of valid RDS files"
    message <- paste0(" test.validateRDSFilesParameters_RDSFiles_empty_array() ",
                        "- Empty array for RDSFiles did not  ",
                        "generated expected message.")
    checkEquals(obs, exp, msg = message)
}


#########################################################
## postMerge() function
#########################################################


## Test the result of postMerge() function
test.postMerge_good_01 <- function() {
    obs <- RJMCMC:::postMerge(startPosForwardReads = reads_demo$readsForward,
                                startPosReverseReads = reads_demo$readsReverse,
                                resultRJMCMC = RJMCMC_result,
                                extendingSize = 10, chrLength = 80000)
    exp <- c(72434.766272478853, 72544.048047704578, 73146.590899701128)
    message <- paste0(" test.postMerge_good_01() ",
                      "- postMerge() did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}


## Test the result of postMerge() function
test.postMerge_good_02 <- function() {
    obs <- RJMCMC:::postMerge(startPosForwardReads = reads_demo$readsForward,
                              startPosReverseReads = reads_demo$readsReverse,
                              resultRJMCMC = RJMCMC_result,
                              extendingSize = 100, chrLength = 80000)
    exp <- c(72452.452375092398, 73146.590899701128)
    message <- paste0(" test.postMerge_good_02() ",
                      "- postMerge() did not generated expected message.")
    checkEquals(obs, exp, msg = message)
}
