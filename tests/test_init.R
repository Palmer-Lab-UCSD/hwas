
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
    Expect$true(hwas:::is.err(m))
    Expect$eq(attr(m, "class"), "error")

    Expect$eq(m, answer_to_ultimate_question)
    Expect$true("msg" %in% names(attributes(m)))
    Expect$eq(attr(m, "msg"), ultimate_question)
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
               "test_get_chromname",
               function() {
    s <- "/path/to/chr01-genotypes.bcf"
    Expect$eq(hwas:::get_chromname(s), "chr01")

    s <- "/path/to/ch01-genotypes.bcf"
    Expect$true(hwas:::is.err(hwas:::get_chromname(s)))

    s <- "/path/to/achr01-genotypes.bcf"
    Expect$true(hwas:::is.err(hwas:::get_chromname(s)))

    s <- "/path/to/chr01a-genotypes.bcf"
    Expect$eq(hwas:::get_chromname(s), "chr01")

    s <- "/path/to/chr01a-genotypes"
    Expect$true(hwas:::is.err(hwas:::get_chromname(s)))
})


unittests$TEST("InitUtils",
               "test_samp_config",
               function() {
    tdir <- tempdir()

    sfname <- "sample_bad"
    sfpath <- file.path(tdir, sfname)
    file.create(sfpath)
    cfg <- hwas:::config_samples(tdir, sfname)
    print(cfg)
    Expect$eq(cfg$file, sfname)
    Expect$false(cfg$inclusion)


    dir.create(file.path(tdir, "tmp"))
    sfname <- "tmp/samples_good.tsv"
    sfpath <- file.path(tdir, sfname)
    file.create(sfpath)
    cfg <- hwas:::config_samples(tdir, sfname)
    print(cfg)
    Expect$eq(cfg$file, sfname)
    Expect$true(cfg$inclusion)
})


# setUp <- function() {
#     tmp_dirname <- tempdir()
# 
# 
#     return(tmp_dirname)
# }
# 
# cleanUp <- function() {
# 
# }
# 
# unittests$TEST("InitDirConstruction",
#                "test_make_dirs",
#                function() {
#     settings <- setUp()
#     on.exit(cleanUp(settings))
# 
#     curr_dir <- getwd()
#     on.exit(setwd(curr_dir))
# 
#     tmpdir <- tempdir()
#     setwd(tmpdir)
# 
#     dir.create("dummy_dir", mode="750")
#     Expect$false(hwas:::make_analysis_directories("dummy_dir"))
# 
# 
#     hwas:::make_analysis_directories("bmi")
# })

if (!interactive())
    unittests$main()
