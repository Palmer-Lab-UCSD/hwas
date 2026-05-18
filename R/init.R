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
#    |--- scripts/
#    |    |--- compute_grm.R
#    |    |--- compute_herit.R
#    |    |--- process_pos.R
#    |    |--- process_traits.R
#    |    |--- unique_samples.R
#    |
#    |--- preprocess_data/
#    |    |--- *traits.tsv*
#    |    |--- samples/
#    |    |    |--- **chr1.exclude_samples**
#    |    |    |--- **chr2.exclude_samples**
#    |    |    |--- ...
#    |    |    |--- **chr2.exclude_samples**
#    |    |--- pos/
#    |         |--- **chr1_pos.exclude_snps**
#    |         |--- **chr2.exclude_snps**
#    |         |--- ...
#    |         |--- **chr2.exclude_snps**
#    |
#    |--- postproces_data/
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
#       ([^a-z]|^)chr[0-9]+.*\\.bcf
#


# TODO :
#   * I need to find position files per chromosome

library(yaml)


###########################################################
## SUFFIXES
###########################################################

SUFFIXES <- new.env(parent = baseenv())

assign("samples", 
       structure("exclude_samples",
                 inclusion = FALSE,
                 class = "suffix_samples"),
       envir = SUFFIXES)
assign("pos",
       structure("exclude_snps",
                 inclusion = FALSE,
                 class = "suffix_pos"),
       envir = SUFFIXES)
assign("geno", "bcf", envir = SUFFIXES)
assign("index", "csi", envir = SUFFIXES)


###########################################################
## REGEX PATTERNS
###########################################################

REGEX <- new.env(parent = baseenv())

assign("invalid_dirchars",
       "[^a-zA-Z0-9-_/]",
        envir = REGEX)
assign("chrom_from_bcfname",
       "([^a-z]|^)(?<chrom>chr[0-9]+).*\\.bcf", 
       envir = REGEX)
assign("find_bcf",
       "([^a-z]|^)chr[0-9]+.*\\.bcf$",
       envir = REGEX)
assign("geno_suffix",
       sprintf("\\.%s$", get("geno", SUFFIXES)),
       envir = REGEX)
assign("sample_filename",
       sprintf("([^a-zA-Z0-9]|^)chr([0-9]){1,2}[^0-9]*\\.%s$", 
               get("samples", SUFFIXES)),
       envir = REGEX)
assign("pos_suffix",
       sprintf("\\.%s$", get("pos", SUFFIXES)),
       envir = REGEX)


###########################################################
## UTILS
###########################################################

is_valid_dirname <- function(d) {
    ismatch <- grep(get("invalid_dirchars", REGEX),
                    d, perl=TRUE)
    return(length(ismatch) == 0)
}


err <- function(val, msg) {
    return(structure(val,
                     msg = msg,
                     class = "error"))
}


is_err <- function(error) {
    return(class(error) == "error"
           && "msg" %in% names(attributes(error)))
}


###########################################################
## CONFIG TOOLS
###########################################################

get_chromname <- function(bcf_filename) {
    bcf_filename <- basename(bcf_filename)
    reg <- regexpr(get("chrom_from_bcfname", REGEX),
                   bcf_filename, 
                   perl=TRUE)
    if (reg != 1)
        return(err(character(0), 
                   "Couldn't extract autosome name from path"))

    start_idx <- attr(reg, "capture.start")[, "chrom"]
    end_idx <- start_idx + attr(reg, "capture.length")[,"chrom"] - 1

    return(substr(bcf_filename, start_idx, end_idx))
}


#' @description
#'  Find requisite sample files assuming that they satisfy the
#'  following assumptions:
#'      * Given the geno_root and sampdir, all exluce sample files
#'          can be found in any subdirectory ${geno_root}/${sampdir}
#'      * The set of child directories in ${geno_root}/${sampdir}
#'          consist only those that contain exclude sample files for
#'          that round of imputation.
#'      * The sample files must contain the chromosome it applies to
#'          in either the relative path, i.e. 
#'          ${geno_root}/${sampdir}/<RELATIVE_PATH> or file name 
#'      * The file extension must be ".exclude_samples"
#' 
#' An example file directory is:
#' 
#'  /path/to/data/
#'           |--- filtered_samples/
#'           |    |--- annot_chr1_filt/
#'           |    |    |--- bad_samples.exclude_samples
#'           |    |    |--- ...
#'           |    |
#'           |    |--- annot_chr2_filt/
#'           |    |    |--- bad_samples.exclude_samples
#'           |    |    |--- ...
#'           |    |
#'           |    |--- ...
#'           |    |
#'           |    |--- annot_chr20_filt/
#'           |         |--- bad_samples.exclude_samples
#'           |         |--- ...
#'           |
#'           |--- genotypes/
#'
#' @return list(dir, filenames)
config_samples <- function(geno_root, sampdir) {

    samp_path <- file.path(geno_root, sampdir)
    if(!dir.exists(samp_path))
        return(err(list(), "Sample directory does not exist"))

    samp_files <- list.files(samp_path, 
                             pattern = get("samples", SUFFIXES),
                             recursive = TRUE)
    samp_files <- grep(get("sample_filename", REGEX),
                       samp_files,
                       perl = TRUE,
                       value = TRUE)
 
    if (length(samp_files) <= 0)
        return(err(list(), "No sample exclude files found"))

    samp_files <- sort(samp_files)

    return(list(dir = sampdir, filenames = samp_files))
}


#' @return list(list(name =, dir=, pos=, pos_exclude=, bcf=, index=),
#'              list(name =, dir=, pos=, pos_exclude=, bcf=, index=),
#'              ...)
config_chroms <- function(genotype_rootdir) {

    bcf_files <- grep(get("find_bcf", REGEX),
                      list.files(genotype_rootdir,
                                  pattern=get("geno_suffix", REGEX),
                                  recursive = TRUE),
                      perl=TRUE,
                      value=TRUE)

    if (length(bcf_files) == 0)
        return(err(list(), "No autosome genotype files with 
                            suffix .bcf were found."))

    chrom_cfg <- vector(mode=list, length=length(bcf_files))

    # Note that bcf_files should include a directory prefix if
    # it was found in a subdirectory of genotype_rootdir
    i <- 1
    for (bcf in bcf_files) {

        if(!is_bcf(bcf))
            return(err(list(), 
                       sprintf("Invalide bcf: %s\n", 
                               file.path(genotype_rootdir, bcf))))

        tmp <- get_chromname(basename(bcf))
        if (is_err(tmp))
            return(tmp)

        chrom_i <- list(name = tmp, 
                        dir = dirname(bcf),
                        pos = NULL,
                        bcf = bcf,
                        index = NULL)

        pos_file <- grepv(get("pos_suffix", REGEX), 
                          list.files(chrom_path, 
                                     pattern=get("pos",SUFFIXES),
                                     recursive = TRUE),
                          perl = TRUE)
        
        if (length(pos_file) == 1)
            chrom_i[["pos"]] <- list(filename = pos_file,
                                     exclusion = TRUE)

        bcf_idx <- sprintf("%s.%s", bcf, get("index", SUFFIXES))
        if (file.exists(file.path(genotype_rootdir, bcf_idx)))
            chrom_i[["index"]] <- bcf_idx

        chrom_cfg[[i]] <- chrom_i
        i <- i + 1
    }

    return(chrom_cfg)
}


mk_config <- function(phenotype,
                      genotype_rootdir,
                      samples_file,
                      pos_include = TRUE) {

    pkg_path <- system.files(".", package="hwas")
    cfg_fname <- file.path(pkg_path, "templates", "config.yaml")
    cfg <- yaml::yaml.load_file(cfg_fname)

    cfg$setup$rootdir <- getwd()
    cfg$phenotype$name <- phenotype   

    if (!dir.exists(genotype_rootdir))
        return(err(list(),
                   sprintf("Genotype directory, %s, is not found\n", 
                           genotype_rootdir)))
    cfg$genotypes$dir <- genotype_rootdir

    samp_cfg <- config_samples(cfg$genotypes$dir, samples)
    if (is.err(samp_cfg))
        return(samp_cfg)

    cfg$genotypes$samples <- samp_cfg

    cfg$genotypes$chroms <- config_chroms(cfg$genotypes$dir)

    cat(yaml::as.yaml(cfg))
    return(cfg)
}

# TODO make skeleton of directories and unit tests

# @title: make the directories and files
# make_skeleton <- function(phenotype,
#                           genotype_rootdir,
#                           pos_file,
#                           samples_file,
#                           pos_include,
#                           samples_include) {
# 
#     # Cases:
#     #   dir.exists      dir.create      statement
#     #   false           false           true && true    -> true
#     #   false           true            true && false   -> false
#     #   true            short circuit   false           -> false
#     #   true            short circuit   false           -> false 
#     if (!dir.exists(pheno) && !dir.create(pheno, mode="0750"))
#         return(err(FALSE, "Couldn't create directory\n"))
# 
#     j <- 1
#     analysis_dirname <- mk_dirname(j)
#     if (!analysis_dirname)
#         return(analysis_dirname)
#     analysis_dirname <- file.path(phenotype, analysis_dirname)
# 
#     while (dir.exists(analysis_dirname)) {
#         j <- j+1
#         analysis_dirname <- mk_analysis_set_dirname(j)
# 
#         if (!analysis_dirname)
#             return(analysis_dirname)
# 
#         analysis_dirname <- file.path(phenotype, analysis_dirname)
#     }
# 
#     if (!dir.create(analysis_dirname, mode="750"))
#         return(err(FALSE,
#                    sprintf("Couldn't create dir %s\n",
#                            analysis_dirname)))
# 
#     if (!dir.create(file.path(analysis_dirname, "grms"), 
#                     mode="750"))
#         return(err(FALSE, "Couldn't create dir grms\n"))
# 
#     if (!dir.create(file.path(analysis_dirname, "preprocess_data"), 
#                     mode="750"))
#         return(err(FALSE, "Couldn't create dir preprocess_data\n"))
# 
#     if (!dir.create(file.path(analysis_dirname, "postprocess_data"), 
#                     mode="750"))
#         return(err(FALSE, "Couldn't create dir postprocess_data\n"))
# 
#     if (!dir.create(file.path(analysis_dirname, "pos"), 
#                     mode="750"))
#         return(err(FALSE, "Couldn't create dir pos\n"))
# 
#     if (!dir.create(file.path(analysis_dirname, "scripts"), 
#                     mode="750"))
#         return(err(FALSE, "Couldn't create dir scripts\n"))
# 
#     if (!dir.create(file.path(analysis_dirname, "results"), 
#                     mode="750"))
#         return(err(FALSE, "Couldn't create dir results\n"))
#     
#     return(TRUE)
# }

###########################################################
## USER ACCESSIBLE FEATURES
###########################################################

# TODO finish init implementation and unit tests
init <- function(phenotype,
                 genotype_rootdir,
                 samples_file,
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
    # tmp_dirname <- make.skeleton(phenotype,
    #                              genotype_rootdir,
    #                              pos_file,
    #                              samples_file,
    #                              pos_include,
    #                              samples_include)

    # if (is.null(tmp_dirname) && is.err(tmp_dirname))
    #     stop(attr(tmp_dirname, "msg"))
    # else if (is.null(tmp_dirname))
    #     stop("couldn't initialize project")




    # ## copy scripts
    # scripts_path <- system.files("scripts", package = "hwas")
    # for (w in list.files(scripts_path, "*.R"))
    #     file.copy(file.path(scripts_path, w), ".",
    #               copy.mode=TRUE,
    #               overwrite=FALSE)


    # for (elem in list.files(tmp)) {
    #     status <- file.copy(file.path(tmp, elem), 
    #                         ".", 
    #                         overwrite=FALSE, 
    #                         recursive=TRUE)
    #     
    #     if (!status)
    #         stop("failed copying hwas skeleton to current directory")
    # }

}
