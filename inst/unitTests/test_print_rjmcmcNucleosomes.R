###################################################
# Created by Astrid Deschenes
# 2015-03-08
###################################################

###################################################
## Test the print.rjmcmcNucleosomes.R functions
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "RJMCMC" )
}

### }}}

data(RJMCMC_result)

###########################################################
## print.rjmcmcNucleosomes() function
###########################################################

test.print_rjmcmcNucleosomes_test_returned_value <- function() {
    result <- print(RJMCMC_result)

    message     <- paste0(" test.print_rjmcmcNucleosomes_test_returned_value() ",
                          "- print method did not returned expected value")

    checkEquals(RJMCMC_result, result, msg = message)
}
