
library(hwas)

source(system.file("utils/unittest.R", package = "hwas"))

unittests$TEST("InitUtils",
               "test_is_valid_dirname", 
               function() {
    Expect$true( hwas:::is_valid_dirname("grms"))
    Expect$true( hwas:::is_valid_dirname("trait/values"))
    Expect$true( hwas:::is_valid_dirname("trait_values"))
    Expect$true( hwas:::is_valid_dirname("trait-values"))
    Expect$false(hwas:::is_valid_dirname("trait.values"))
    Expect$false(hwas:::is_valid_dirname("trait?values"))
    Expect$false(hwas:::is_valid_dirname("trait\values"))
})


unittests$TEST("InitUtils",
               "test_err",
               function() {
    # the ultimate question and associated answer is from the book
    # hitch hiker's guide to the galaxy.
    ultimate_question <- "the answer to life, the universe,
        and everything"
    answer_to_ultimate_question <- 42

    m <- hwas:::err(answer_to_ultimate_question, ultimate_question)
    Expect$eq(m, answer_to_ultimate_question)
    Expect$eq(attr(m, "msg"), ultimate_question)
    Expect$eq(attr(m, "class"), "error")
})


unittests$TEST("InitUtils",
               "test_mk_dirname",
               function() {

    Expect$false(hwas:::mk_dirname(-1))
    Expect$false(hwas:::mk_dirname(c(1,2,4)))
    Expect$false(hwas:::mk_dirname(54.2))
    Expect$false(hwas:::mk_dirname(1000))
    Expect$false(hwas:::mk_dirname(0))


    year_str <- format(Sys.time(), "%Y")
    Expect$eq(hwas:::mk_dirname(1),
              paste(year_str, "001", sep="-"))

    Expect$eq(hwas:::mk_dirname(99),
              paste(year_str, "099", sep="-"))
})                  


unittests$TEST("InitUtils",
               "test_is_error",
               function() {
    v <- 42
    s <- "the answer to life, the universe, everything"
    error <- hwas:::err(v, s)

    Expect$true(hwas:::is.err(error))

    Expect$false(hwas:::is.err(v))
    Expect$false(hwas:::is.err(structure(v, class = "error")))
    Expect$false(hwas:::is.err(structure(v, msg = s)))
})


unittests$TEST("InitUtils",
               "test_err",
               function() {
    v <- character(0)
    s <- "hi!"
    error <- hwas:::err(v, s) 

    Expect$eq(class(error), "error")
    Expect$eq(length(error), 0)
    
    Expect$true(hwas:::is.err(error))

    Expect$true("msg" %in% names(attributes(error)))
    Expect$eq(attr(error, "msg"), s)

    Expect$true(is.character(error))
})


unittests$TEST("InitUtils",
               "test_get_chromname",
               function() {
    Expect$eq(hwas:::get_chromname("/path/to/chr01-genotypes.bcf"),
              "chr01")

    s <- "/path/to/cr01-genotypes.bcf"
    Expect$true(hwas:::is.err(hwas:::get_chromname(s)))
})


setUp <- function() {
    tmp_dirname <- tempdir()


    return(tmp_dirname)
}

cleanUp <- function() {

}

unittests$TEST("InitDirConstruction",
               "test_make_dirs",
               function() {
    settings <- setUp()
    on.exit(cleanUp(settings))

    curr_dir <- getwd()
    on.exit(setwd(curr_dir))

    tmpdir <- tempdir()
    setwd(tmpdir)

    dir.create("dummy_dir", mode="750")
    Expect$false(hwas:::make_analysis_directories("dummy_dir"))


    hwas:::make_analysis_directories("bmi")
})

if (!interactive())
    unittests$main()
