
library(hwas)

source(system.file("utils/unittest.R", package = "hwas"))


cleanUp <- function(tmp_dirname) {
    if (regexpr("^/tmp/Rtmp[a-zA-Z0-9]+$", tmp_dirname) > 0)
        unlink(tmp_dirname, recursive = TRUE)
    else
        cat(sprintf("Directory, %s, is not valid tmp directory\n",
                    tmp_dirname))
}


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

    tdir <- tempdir(check=TRUE)
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


unittests$TEST("InitConfig",
               "test_samp_config",
               function() {
    tmp_cfg <- setUpSamples()
    on.exit(cleanUp(tmp_cfg$temp_dir))

    cfg <- hwas:::config_samples(tmp_cfg$temp_dir, "wrong_sampdir")
    Expect$true(hwas:::is_err(cfg))


    cfg <- hwas:::config_samples(tmp_cfg$temp_dir, tmp_cfg$samp_dir)
    Expect$false(hwas:::is_err(cfg))

    Expect$eq(cfg$dir, tmp_cfg$samp_dir)
    Expect$eq(length(cfg$filenames), length(tmp_cfg$filenames))

    for (i in seq(length(tmp_cfg$filenames)))
        Expect$eq(cfg$filenames[i], tmp_cfg$filenames[i])
})


# /tmp/Rtmpaijfdad/
#      |--- filtered_snps/
#           |--- annot_chr1_randomstr/
#                |--- pos.exclude_snps
#                |--- geno_test_data.bcf
#                |--- geno_test_data.bcf.csi
setUpGenotypes <- function() {
    data_dir <- system.file("exdata", package = "hwas")

    # tempdir has an odd behavior that only one temp directory is
    # made per session.  Consequently, my cleanUp is deleting it.
    # By setting check = TRUE, I can delete the pervious tempdir
    # and create a replacement one here.
    tdir <- tempdir(check = TRUE)
    genotypes_dir <- "filt_snps"
    
    if (!dir.create(file.path(tdir, genotypes_dir))) {
        cat(sprintf("Could not create directory %s\n", 
                    file.path(tdir, genotypes_dir)), 
            file = stderr())
        return(list())
    }

    out <- list(round_dir = tdir,
                dir = genotypes_dir,
                chr = list())

    for (i in seq(10)) {
        chrom <- list(name = sprintf("chr%d", i),
                      dir = sprintf("auto_chr%d_aadsf", i),
                      pos = "pos.exclude_snps",
                      bcf = "geno_test_data_compressed.bcf",
                      index = "geno_test_data_compressed.bcf.csi")

        dir.create(file.path(out$round_dir, out$dir, chrom$dir))

        file.copy(from = file.path(data_dir, chrom$pos),
                  to = file.path(out$round_dir,out$dir, chrom$dir),
                  overwrite = FALSE)
        file.copy(from = file.path(data_dir, chrom$bcf),
                  to = file.path(out$round_dir,out$dir, chrom$dir),
                  overwrite = FALSE)
        file.copy(from = file.path(data_dir, chrom$index),
                  to = file.path(out$round_dir,out$dir, chrom$dir),
                  overwrite = FALSE)

        out$chr[[i]] <- chrom
    }

    return(out)
}


unittests$TEST("InitConfig",
               "test_chrom",
               function() {
    tmp_cfg <- setUpGenotypes()
    on.exit(cleanUp(tmp_cfg$round_dir))

    out <- hwas:::config_chroms(tmp_cfg$round_dir, tmp_cfg$dir)

    print(length(tmp_cfg$chr))
    print(length(out))
    Expect$eq(length(tmp_cfg$chr), length(out))

    for (truth_chr in tmp_cfg$chr) {

        for (out_chr in out)
            if (truth_chr$name == out_chr$name) {
                Expect$eq(truth_chr$dir,   out_chr$dir)
                Expect$eq(truth_chr$pos,   out_chr$pos)
                Expect$eq(truth_chr$bcf,   out_chr$bcf)
                Expect$eq(truth_chr$index, out_chr$index)
                break
            }
    }
 })


setUpCfg <- function() {
}

unittests$TEST("InitConfig",
               "test_mk_config",
               function() {

    tmp_cfg <- setUpGenotypes()
    on.exit(cleanUp(tmp_cfg$round_dir))
    
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
