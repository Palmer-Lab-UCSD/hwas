# Initialize HWAS run
#
# The init function initializes an HWAS run by making the expected
# directories, copy scripts that are used in the pipeline, and set
# up configuration files.  The other functions in this file support
# the init function and are not exported.
# 
# The structure of an hwas run is as follows:
#
# |- analysis_dir/
#    |--- config.yaml
#    |--- .scripts/
#    |    |--- compute_grm.R
#    |    |--- compute_herit.R
#    |    |--- process_pos.R
#    |    |--- process_traits.R
#    |    |--- unique_samples.R
#    |
#    |--- preprocess_data/
#    |    |--- *traits.tsv*
#    |    |--- *samples*
#    |
#    |--- postprocess_data/
#    |    |--- **covariates.tsv**
#    |    |--- **phenotype.tsv**
#    |    |--- **samples**
#    |    |--- pos/
#    |         |--- **chr01.tsv**
#    |         |--- **chr02.tsv**
#    |         |--- ...
#    |         |--- **chr20.tsv**
#    |
#    |--- results/
#         |--- **heritability**
#         |--- grms/
#         |    |--- **chr01.RData**
#         |    |--- **chr02.RData**
#         |    |--- ...
#         |    |--- **chr20.RData**
#         |
#         |--- lod/
#         |    |--- **chr01.bed**
#         |    |--- ...
#         |    |--- **chr20.bed**
#         |
#         |--- blup/
#              |--- **chr01.tsv**
#              |--- ...
#              |--- **chr20.tsv**
#
# NOTE ON DIR STRUCTURE
# - *traits.tsv* a file with the sample id, covariates, and trait
#       measurements.  The user is responsible for making this file,
#       and it is considered raw data.
#
# - **covariates.tsv**, and any other name between ** are produced
#       by pipeline.  
#
#
# ASSUMPTIONS
#
# * Configuring genotype files
#   - user specified path to genotypes contains all genotype files 
#       desired
#   - genotypes are in the bcf file format with suffix .bcf
#   - bcf file names must have the chromosome string in format
#       chr[0-9]+[^a-zA-Z]+.
#


# TODO :
#   * I need to find position files per chromosome

library(yaml)

is_valid_dirname <- function(d) {
    ismatch <- grep("[^a-zA-Z0-9-_/]", d, perl=TRUE)
    return(length(ismatch) == 0)
}


is.err <- function(error) {
    return(class(error) == "error"
           && "msg" %in% names(attributes(error)))
}


err <- function(val, msg) {
    return(structure(val,
                     msg = msg,
                     class = "error"))
}


mk_dirname <- function(i) {
    if (length(i) != 1 
        || !is.numeric(i) 
        || floor(i) == 0
        || i %% floor(i) != 0)
        return(err(FALSE, "Index needs to be a natural number"))

    if (i <= 0 || i >= 1000)
        return(err(FALSE, "Index out of bounds"))

    return(sprintf("%s-%03d", format(Sys.time(), "%Y"), i))
}


get_chromname <- function(bcf_filename) {
    bcf_filename <- basename(bcf_filename)
    reg <- regexpr("^(?<chrom>chr[0-9]+)[^a-zA-Z]", 
                   bcf_filename, 
                   perl=TRUE)
    if (reg != 1)
        return(err(character(0), 
                   "Couldn't extract autosome name from path"))

    start_idx <- attr(reg, "capture.start")[, "chrom"]
    end_idx <- start_idx + attr(reg, "capture.length")[,"chrom"] - 1

    return(substr(bcf_filename, start_idx, end_idx))
}



# @title: make the directories and files
make_skeleton <- function(phenotype,
                          genotype_rootdir,
                          pos_file,
                          samples_file,
                          pos_include,
                          samples_include) {

    # Cases:
    #   dir.exists      dir.create      statement
    #   false           false           true && true    -> true
    #   false           true            true && false   -> false
    #   true            short circuit   false           -> false
    #   true            short circuit   false           -> false 
    if (!dir.exists(pheno) && !dir.create(pheno, mode="0750"))
        return(err(FALSE, "Couldn't create directory\n"))

    j <- 1
    analysis_dirname <- mk_dirname(j)
    if (!analysis_dirname)
        return(analysis_dirname)
    analysis_dirname <- file.path(phenotype, analysis_dirname)

    while (dir.exists(analysis_dirname)) {
        j <- j+1
        analysis_dirname <- mk_analysis_set_dirname(j)

        if (!analysis_dirname)
            return(analysis_dirname)

        analysis_dirname <- file.path(phenotype, analysis_dirname)
    }

    if (!dir.create(analysis_dirname, mode="750"))
        return(err(FALSE,
                   sprintf("Couldn't create dir %s\n",
                           analysis_dirname)))

    if (!dir.create(file.path(analysis_dirname, "grms"), 
                    mode="750"))
        return(err(FALSE, "Couldn't create dir grms\n"))

    if (!dir.create(file.path(analysis_dirname, "preprocess_data"), 
                    mode="750"))
        return(err(FALSE, "Couldn't create dir preprocess_data\n"))

    if (!dir.create(file.path(analysis_dirname, "postprocess_data"), 
                    mode="750"))
        return(err(FALSE, "Couldn't create dir postprocess_data\n"))

    if (!dir.create(file.path(analysis_dirname, "pos"), 
                    mode="750"))
        return(err(FALSE, "Couldn't create dir pos\n"))

    if (!dir.create(file.path(analysis_dirname, "scripts"), 
                    mode="750"))
        return(err(FALSE, "Couldn't create dir scripts\n"))

    if (!dir.create(file.path(analysis_dirname, "results"), 
                    mode="750"))
        return(err(FALSE, "Couldn't create dir results\n"))
    
    return(TRUE)
}




mk_config <- function(cfg_fname, 
                      phenotype,
                      genotype_rootdir,
                      pos_file,
                      samples_file,
                      pos_include = TRUE,
                      samples_include = TRUE) {

    cfg <- yaml::yaml.load_file(cfg_fname)

    cfg$setup$rootdir <- getwd()
    cfg$phenotype$name <- phenotype   


    cfg$genotypes$dir <- genotype_rootdir
    cfg$genotypes$pos_include <- NULL
    cfg$genotypes$pos_exclude <- NULL
    if (pos_include)
        cfg$genotypes$pos_include <- pos_file
    else
        cfg$genotypes$pos_exclude <- pos_file

    cfg$genotypes$samples_include <- NULL
    cfg$genotypes$samples_exclude <- NULL

    if (samples_include)
        cfg$genotypes$samples_include <- samples_file
    else
        cfg$genotypes$samples_exclude <- samples_file


    bcf_files <- grepv("chr[0-9]+.*\\.bcf$", 
                       list.files(cfg$genotypes$dir,
                                  pattern=".bcf",
                                  recursive = TRUE),
                       perl=TRUE)

    if (length(bcf_files) == 0)
        return(err(character(0), "No autosome genotype files with 
                     suffix .bcf were found."))

    cfg$genotypes$chrom <- vector(mode=list, length=length(bcf_files))

    # Note that bcf_files should include a directory prefix if
    # it was found in a subdirectory of cfg$genotypes$dir
    i <- 1
    for (bcf in bcf_files) {
        chrom <- get_chromname(basename(bcf))
        if (is.null(chrom))
            return(chrom)

        bcf_idx <- paste0(bcf, ".csi")
        if (file.exists(file.path(cfg$genotypes$dir, bcf_idx)))
            bcf_idx <- basename(bcf_idx)
        else
            bcf_idx <- NULL

        cfg$genotypes$chrom[[i]] <- list(name: chrom,
                                         dir: dirname(bcf),
                                         bcf: basename(bcf),
                                         index: bcf_index)
        i <- i + 1
    }

    return(cfg)
}


init <- function(phenotype,
                 genotype_rootdir,
                 pos_file,
                 samples_file,
                 pos_include = TRUE,
                 samples_include = TRUE) {

    if (length(list.files(".")) != 0)
        stop("Make a new directory with no contents.\n")

    pkg_path <- system.files(".", package="hwas")
    cfg <- mk_config(file.path(pkg_path, "templates", "config.yaml"),
                     genotype_rootdir,
                     pos_file,
                     samples_file,
                     pos_include = TRUE,
                     samples_include = TRUE)

    # Remember: we make the directory in a temporary directory, and
    # move upon successful setup 
    tmp_dirname <- make.skeleton(phenotype,
                                 genotype_rootdir,
                                 pos_file,
                                 samples_file,
                                 pos_include,
                                 samples_include)

    if (is.null(tmp_dirname) && is.err(tmp_dirname))
        stop(attr(tmp_dirname, "msg"))
    else if (is.null(tmp_dirname))
        stop("couldn't initialize project")




    ## copy scripts
    scripts_path <- system.files("scripts", package = "hwas")
    for (w in list.files(scripts_path, "*.R"))
        file.copy(file.path(scripts_path, w), ".",
                  copy.mode=TRUE,
                  overwrite=FALSE)


    for (elem in list.files(tmp)) {
        status <- file.copy(file.path(tmp, elem), 
                            ".", 
                            overwrite=FALSE, 
                            recursive=TRUE)
        
        if (!status)
            stop("failed copying hwas skeleton to current directory")
    }

}
