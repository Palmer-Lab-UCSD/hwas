
######################################################################
# AI CONTRIBUTION
# CLAUDE Opus 4.7 suggested using an environment for this global
# acculator variable.  Previously, I was using a list, which did
# not update correctly.
RESULTS <- new.env(parent = emptyenv())
RESULTS$failed <- vector(mode="character", length=0)
RESULTS$warning <- vector(mode="character", length=0)
RESULTS$success <- vector(mode="character", length=0)
# END AI CONTRIBUTION
######################################################################



TEST_NULL <- function(x) {
    caller <- sys.call(-1)

    if (is.null(x)) {
        RESULTS$success <- append(RESULTS$success, caller)
        return(TRUE)
    }

    print(paste("FAILED: ", caller, "NULL value not returned."))

    RESULTS$failed <- append(RESULTS$failed, caller)
    return(FALSE)
}

TEST_NOTNULL <- function(x) {
    caller <- sys.call(-1)

    if (!is.null(x)) {
        RESULTS$success <- append(RESULTS$success, caller)
        return(TRUE)
    }

    print(paste("FAILED: ", caller, "NULL value returned."))

    RESULTS$failed <- append(RESULTS$failed, caller)
    return(FALSE)
}


TEST_NEQ <- function(test_val, truth_val) {
    caller <- sys.call(-1)

    if (test_val != truth_val) {
        RESULTS$success <- append(RESULTS$success, caller)
        return(TRUE)
    }

    print(paste("FAILED:", caller,
           "values are equal when not equal expected"))
    
    RESULTS$failed <- append(RESULTS$failed, caller)
    return(FALSE)
}


TEST_EQ <- function(test_val, truth_val) {
    caller <- sys.call(-1)

    if (test_val == truth_val) {
        RESULTS$success <- append(RESULTS$success, caller)
        return(TRUE)
    }

    print(paste("FAILED:", 
          sys.call(-1), 
          "values are not equal when equality expected"))
    
    RESULTS$failed <- append(RESULTS$failed, caller)
    return(FALSE)
}


ERROR_REPORT <- function() {
    cat("===============================================\n")
    cat("Test results\n")
    cat(sprintf("Failed%20d\n", length(RESULTS$failed)))
    for (j in seq(length(RESULTS$failed))) 
        cat(sprintf("%4s\n", RESULTS$failed[j]))

    cat(sprintf("Warning%20d\n", length(RESULTS$warning)))
    for (j in seq(length(RESULTS$warning))) 
        cat(sprintf("%4s\n", RESULTS$warning[j]))

    cat(sprintf("Success%20d\n", length(RESULTS$success)))

    cat("-----------------------------------------------\n")
    if (length(RESULTS$failed) == 0)
        cat("PASS\n")
    else
        cat("FAILED\n")

    cat("===============================================\n")
    if (length(RESULTS$failed) != 0)
        stop()
}
