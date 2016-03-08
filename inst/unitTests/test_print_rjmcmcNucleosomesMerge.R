###################################################
# Created by Astrid Deschenes
# 2015-03-08
###################################################

###################################################
## Test the print.rjmcmcNucleosomesMerge.R function
###################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "RJMCMC" )
}

### }}}


###########################################################
## print.rjmcmcNucleosomesMerge() function
###########################################################

test.print_rjmcmcNucleosomesMerge_test_returned_value <- function() {

    directoryWithRDSFiles <- system.file("extdata", package = "RJMCMC")

    resultInit <- mergeAllRDSFilesFromDirectory(directoryWithRDSFiles)

    result <- print(resultInit)

    message <- paste0(" test.print_rjmcmcNucleosomesMerge_test_returned_value() ",
                          "- print method did not returned expected value")

    checkEquals(resultInit, result, msg = message)
}
