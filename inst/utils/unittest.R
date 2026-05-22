# In these unittest tools I've tried to reproduce Google Test
# nomenclature.  
# unittest environment holds the status --failed, warning and 
# success -- accumulators, expect statements, and test suite.  
# The expect statements are individual tests that update the 
# accumulators with the relevant status of the test.  The
# test suite is the collection of all user tests defined tests
# and the main function.  The main runction runs all tests
# and produces the output report.
#
# Note: the .unittests is not apart of the API, but an internal
# environment object that is manipulated


unittests <- new.env(parent = baseenv())
Expect <- new.env(parent = unittests)

###########################################################
## UNITTESTS ENVIRONMENT
###########################################################

evalq({
    .MAGIC_CONST <- "__UNIT_TEST_RESULTS__"
    .PRINT_BARRIER <- "===============================================\n"
    .PRINT_SMALL_BAR <- "-----------------------------------------------\n"
    # user defined test suites and associated tests
    # registered with unittests$TEST
    .test_suites <- list()

    # book keeping test results
    .failed  <- vector(mode="character", length = 0)
    .warning <- vector(mode="character", length = 0)
    .success <- vector(mode="character", length = 0)

    # function to register test suites and associated test functions
    TEST <- function(test_suite, test_name, func) {
        # by the evalq function, the enclosing environment is
        # unittest, therefore the environment within the function
        # scope is local and penv is the unittest env.
        penv <- parent.env(environment())
        if (!(test_suite %in% names(.test_suites)))
            assign(".test_suites",
                   append(.test_suites, 
                          stats::setNames(list(list()), 
                                          test_suite)),
                   envir = penv)

        ts <- get(".test_suites", envir = penv)
        element <- structure(func, 
                             test_suite = test_suite,
                             test_name  = test_name,
                             class      = "test")
        ts[[test_suite]] <- append(ts[[test_suite]], 
                                   stats::setNames(list(element), test_name))
        assign(".test_suites", ts, envir=penv)
    }
}, envir = unittests)


evalq({
    .print_pass_fail <- function(x) {
        cat(.PRINT_SMALL_BAR)
        cat(sprintf("%s\n", x))
        cat(.PRINT_SMALL_BAR)
}}, envir = unittests)


evalq({
    .print_report <- function() {
        penv <- parent.env(environment())


        cat(sprintf("%s\n", .MAGIC_CONST))
        cat(.PRINT_BARRIER)
        cat("TEST RESULTS\n")

        n_failed    <- length(get(".failed",    penv))
        n_warn      <- length(get(".warning",      penv))
        n_success   <- length(get(".success",   penv))
        n_tests     <- n_failed + n_warn + n_success

        if (n_tests == 0) {
            get(".print_pass_fail", penv)("FAILED")
            return(FALSE)
        }
        cat(sprintf("%-20s%d\n", "Failed", n_failed))
        for (msg in get(".failed", penv))
            cat(sprintf("%4s\n", msg))

        cat(sprintf("%-20s%d\n", "Warning", n_warn))
        for (msg in get(".warning", penv))
            cat(sprintf("%4s\n", msg))

        cat(sprintf("%-20s%d\n", "Success", n_success))

        if (n_failed == 0)
            get(".print_pass_fail", penv)("PASS")
        else 
            get(".print_pass_fail", penv)("FAILED")
}}, envir = unittests)

evalq({
    main <- function() {
        penv <- parent.env(environment())

        # run all user defined unit test
        for (test_suite_i in get(".test_suites", penv))
            for (test_i in test_suite_i)
                test_i()

        .print_report()
    }
}, envir = unittests)



###########################################################
## EXPECT 
###########################################################

evalq({true <- function(x) {
    caller <- sys.function(-1)
    penv <- parent.env(environment())
    ppenv <- parent.env(penv)

    s <- sprintf("%s:%s", 
                 attr(caller, "test_suite"),
                 attr(caller, "test_name"))

    if (x) {
        assign(".success",
               append(get(".success", ppenv), s),
               envir = ppenv)
        return(TRUE)
    }
    assign(".failed",
           append(get(".failed", ppenv), s),
           envir = ppenv)
    return(FALSE)
}}, envir = Expect)


evalq({false <- function(x) {
    caller <- sys.function(-1)
    penv <- parent.env(environment())
    ppenv <- parent.env(penv)

    s <- sprintf("%s:%s", 
                 attr(caller, "test_suite"),
                 attr(caller, "test_name"))

    if (!x) {
        assign(".success",
               append(get(".success", ppenv), s),
               envir = ppenv)
        return(TRUE)
    }
    assign(".failed",
           append(get(".failed", ppenv), s),
           envir = ppenv)
    return(FALSE)
}}, envir = Expect)


evalq({null <- function(x) {
    caller <- sys.function(-1)
    penv <- parent.env(environment())
    ppenv <- parent.env(penv)

    s <- sprintf("%s:%s", 
                 attr(caller, "test_suite"),
                 attr(caller, "test_name"))

    if (is.null(x)) {
        assign(".success",
               append(get(".success", ppenv), s),
               envir = ppenv)
        return(TRUE)
    }
    assign(".failed",
           append(get(".failed", ppenv), s),
           envir = ppenv)
    return(FALSE)
}}, envir = Expect)


evalq({not_null <- function(x) {
    caller <- sys.function(-1)
    penv <- parent.env(environment())
    ppenv <- parent.env(penv)

    s <- sprintf("%s:%s", 
                 attr(caller, "test_suite"),
                 attr(caller, "test_name"))

    if (!is.null(x)) {
        assign(".success",
               append(get(".success", ppenv), s),
               envir = ppenv)
        return(TRUE)
    }
    assign(".failed",
           append(get(".failed", ppenv), s),
           envir = ppenv)
    return(FALSE)
}}, envir = Expect)


evalq({neq <- function(test_val, truth_val) {
    caller <- sys.function(-1)
    penv <- parent.env(environment())
    ppenv <- parent.env(penv)

    s <- sprintf("%s:%s", 
                 attr(caller, "test_suite"),
                 attr(caller, "test_name"))

    if (test_val != truth_val) {
        assign(".success",
               append(get(".success", ppenv), s),
               envir = ppenv)
        return(TRUE)
    }
    cat(paste("Expected neq, instead\ntest_val:",
              test_val,
              "\n, equal to truth val:\n",
              truth_val,
              "\n"),
        file = stderr())
    assign(".failed",
           append(get(".failed", ppenv), s),
           envir = ppenv)
    return(FALSE)
}}, envir = Expect)


evalq({eq <- function(test_val, truth_val) {
    caller <- sys.function(-1)
    penv <- parent.env(environment())
    ppenv <- parent.env(penv)

    s <- sprintf("%s:%s", 
                 attr(caller, "test_suite"),
                 attr(caller, "test_name"))

    if (test_val == truth_val) {
        assign(".success",
               append(get(".success", ppenv), s),
               envir = ppenv)
        return(TRUE)
    }
    cat(paste("Expected eq, instead test_val,",
              test_val,
              ",not equal to,",
              truth_val,
              ".\n"),
        file = stderr())
    assign(".failed",
           append(get(".failed", ppenv), s),
           envir = ppenv)
    return(FALSE)
}}, envir = Expect)


