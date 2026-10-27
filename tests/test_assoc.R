

library(hwas)

source(system.file("utils/unittest.R", package = "hwas"))


unittests$TEST("TestAssoc", 
               "test_input_validation", 
               function() {

    stop("Not implemented error")
})



if (!interactive()) {
    print("Running as not interactive")
    unittests$main()
}

