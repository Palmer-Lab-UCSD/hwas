
library(hwas)

RESULTS <- list(failed = vector(mode="character", length=0),
                warning = vector(mode="character", length=0),
                success = vector(mode="character", length=0))



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


set_up <- function() {
    bcf_filename <- system.file("exdata", "geno_test_data.bcf", 
                                package = "hwas",
                                mustWork = TRUE)


    samples <- readLines(system.file("exdata", "samples",
                                     package = "hwas",
                                     mustWork = TRUE))

    loci <- read.table(system.file("exdata", "pos.tsv",
                                   package = "hwas",
                                   mustWork = TRUE),
                        row.names = NULL,
                        header = TRUE)

    tmp <- 42
    return(structure(tmp,
                     samples = samples,
                     chrom = loci$chrom,
                     pos = loci$pos,
                     bcf = bcf_filename,
                     class = "data"))
}


test_bopen <- function() {
    test_data <- set_up()
    bid <- hwas::bopen(attr(test_data, "bcf"), "r")

    TEST_NOTNULL(bid);
    
    status <- hwas::bclose(bid)
    
    TEST_EQ(0, status)
}


if (!interactive()) {
    test_bopen()
    
    ERROR_REPORT()
}


