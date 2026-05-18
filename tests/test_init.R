
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
    Expect$true(hwas:::is_err(m))
    Expect$eq(attr(m, "class"), "error")

    Expect$eq(m, answer_to_ultimate_question)
    Expect$true("msg" %in% names(attributes(m)))
    Expect$eq(attr(m, "msg"), ultimate_question)
})


unittests$TEST("InitUtils",
               "test_is_error",
               function() {
    v <- 42
    s <- "the answer to life, the universe, everything"
    error <- hwas:::err(v, s)

    Expect$true(hwas:::is_err(error))

    Expect$false(hwas:::is_err(v))
    Expect$false(hwas:::is_err(structure(v, class = "error")))
    Expect$false(hwas:::is_err(structure(v, msg = s)))
})



unittests$TEST("InitUtils",
               "test_get_chromname",
               function() {
    s <- "/path/to/chr01-genotypes.bcf"
    Expect$eq(hwas:::get_chromname(s), "chr01")

    s <- "/path/to/ch01-genotypes.bcf"
    Expect$true(hwas:::is_err(hwas:::get_chromname(s)))

    s <- "/path/to/achr01-genotypes.bcf"
    Expect$true(hwas:::is_err(hwas:::get_chromname(s)))

    s <- "/path/to/chr01a-genotypes.bcf"
    Expect$eq(hwas:::get_chromname(s), "chr01")

    s <- "/path/to/chr01a-genotypes"
    Expect$true(hwas:::is_err(hwas:::get_chromname(s)))
})


# tempdir/
#   sample_dir/
#       chr3_filter/
#           samples.exclude_samples
#       chr18_filter/
#           samples.exclude_samples
 
# list () of filenames
setUpSamples <- function() {
    data_dir <- system.file("exdata/", package = "hwas")

    tdir <- tempdir()
    samp_dir <- "samp_dir"
    chr3_dir <- "b_chr3_filter"
    chr18_dir <- "a_chr18_filter"
    
    dir.create(file.path(tdir, samp_dir))
    dir.create(file.path(tdir, samp_dir, chr3_dir))
    dir.create(file.path(tdir, samp_dir, chr18_dir))

    # Note that temp_dir serves the purpose of geno_rootdir
    out <- list(temp_dir= tdir,
                samp_dir = samp_dir,
                filenames = c(file.path(chr3_dir, "samps.exclude_samples"),
                              file.path(chr18_dir, "samps.exclude_samples")))

    out$filenames <- sort(out$filenames)

    for (tofile in out$filenames)
        file.copy(from = file.path(data_dir, "samps.exclude_samples"),
                  to = file.path(out$temp_dir, out$samp_dir, tofile),
                  overwrite = FALSE)

    return(out)
}

cleanUpSamples <- function(tmp_dirname) {
    if (regexpr("^/tmp/Rtmp[a-zA-Z0-9]+$", tmp_dirname) > 0)
        unlink(tmp_dirname, recursive = TRUE)
    else
        cat(sprintf("Directory, %s, is not valid tmp directory\n",
                    tmp_dirname))
}

unittests$TEST("InitConfig",
               "test_samp_config",
               function() {
    tmp_cfg <- setUpSamples()
    on.exit(cleanUpSamples(tmp_cfg$temp_dir))

    cfg <- hwas:::config_samples(tmp_cfg$temp_dir, "wrong_sampdir")
    Expect$true(hwas:::is_err(cfg))


    cfg <- hwas:::config_samples(tmp_cfg$temp_dir, tmp_cfg$samp_dir)
    Expect$false(hwas:::is_err(cfg))

    Expect$eq(cfg$dir, tmp_cfg$samp_dir)
    Expect$eq(length(cfg$filenames), length(tmp_cfg$filenames))

    for (i in seq(length(tmp_cfg$filenames)))
        Expect$eq(cfg$filenames[i], tmp_cfg$filenames[i])
})


setUpGenotypes <- function() {
    data_dir <- system.file("exdata/", package = "hwas")
    tdir <- tempdir()

    print(tdir)

    for (i in seq(10)) {
        chromdir <- sprintf("chr%d_asdfsad", i)
        dir.create(file.path(tdir, chromdir))
        file.copy(file.path(data_dir, "geno_test_data.bcf"),
                  file.path(tdir, chrom, "chrm_genos.bcf"))
    }

    return(dirname)
}
unittests$TEST("InitConfig",
               "test_chrom",
               function() {

 })


# TODO unit tests from chrom some configuration deduction


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
