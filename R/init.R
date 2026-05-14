# Initialize HWAS run
#
# The init function initializes an HWAS run by making the expected
# directories, copy scripts that are used in the pipeline, and set
# up configuration files.  The other functions in this file support
# the init function and are not exported.
# 
# The structure of an hwas run is as follows:
#
# |- phenotype_name/
#    |--- 2026-001/
#    |    |--- config.yaml
#    |    |--- scripts/
#    |    |    |--- compute_grm.R
#    |    |    |--- compute_herit.R
#    |    |    |--- process_pos.R
#    |    |    |--- process_traits.R
#    |    |    |--- unique_samples.R
#    |    | 
#    |    |--- grms/
#    |    |    |--- **chr01.RData**
#    |    |    |--- **chr02.RData**
#    |    |    |--- ...
#    |    |    |--- **chr20.RData**
#    |    |
#    |    |--- preprocess_data/
#    |    |    |--- *traits.tsv*
#    |    |    |--- *samples*
#    |    |
#    |    |--- postprocess_data/
#    |    |    |--- **covariates.tsv**
#    |    |    |--- **phenotype.tsv**
#    |    |    |--- **samples**
#    |    |
#    |    |--- pos/
#    |    |    |--- **chr01.tsv**
#    |    |    |--- **chr02.tsv**
#    |    |    |--- ...
#    |    |    |--- **chr20.tsv**
#    |    |
#    |    |--- results/
#    |    |    |--- **heritability**
#    |    |    |--- lod/
#    |    |    |    |--- **chr01.bed**
#    |    |    |    |--- ...
#    |    |    |    |--- **chr20.bed**
#    |    |    |
#    |    |    |--- blup/
#    |    |    |    |--- **chr01.tsv**
#    |    |    |    |--- ...
#    |    |    |    |--- **chr20.tsv**
#    |
#    |--- YYYY-III/
#
# NOTE ON DIR STRUCTURE
# - YYYY-III is the year in format YYYY the analysis was completed
#   and III is the index, e.g. 002.
# 
# - *traits.tsv* a file with the sample id, covariates, and trait
#       measurements.  The user is responsible for making this file,
#       and it is considered raw data.
#
# - **covariates.tsv**, and any other name between ** are produced
#       by pipeline.  

is_valid_dirname <- function(d) {
    ismatch <- grep("[^a-zA-Z0-9-_/]", d, perl=TRUE)
    return(length(ismatch) == 0)
}


err <- function(val, msg) {
    return(structure(val,
                     msg = msg,
                     class = "error"))
}


mk_analysis_set_dirname <- function(i) {
    if (i <= 0 || i >= 1000)
        return(err(FALSE, "Index out of bounds"))

    return(sprintf("%s-%3d", format(Sys.time(), "%Y"), i))
}


make_analysis_directories <- function(pheno) {

    # Cases:
    #   dir.exists      dir.create      statement
    #   false           false           true && true    -> true
    #   false           true            true && false   -> false
    #   true            short circuit   false           -> false
    #   true            short circuit   false           -> false 
    if (!dir.exists(pheno) && !dir.create(pheno, mode="0750"))
        return(err(FALSE, "Couldn't create directory\n"))

    setwd(phenotype)

    j <- 1
    analysis_dirname <- mk_analysis_set_dirname(j)
    while (dir.exists(analysis_dir)) {
        j <- j+1
        analysis_dirname <- mk_analysis_set_dirname(j)

        if (!analysis_dirname)
            return(analysis_dirname)
    }
    

    if (!dir.create(analysis_dirname, mode="750"))
        return(err(FALSE,
                   sprintf("Couldn't create dir %s\n",
                           analysis_dirname)))

    setwd(analysis_dirname)
    if (!dir.create("grms", mode="750"))
        return(err(FALSE, "Couldn't create dir grms\n"))

    if (!dir.create("preprocess_data", mode="750"))
        return(err(FALSE, "Couldn't create dir preprocess_data\n"))

    if (!dir.create("postprocess_data", mode="750"))
        return(err(FALSE, "Couldn't create dir postprocess_data\n"))

    if (!dir.create("seq", mode="750"))
        return(err(FALSE, "Couldn't create dir seq\n"))

    if (!dir.create(file.path("seq", "pos"), mode="750"))
        return(err(FALSE, "Couldn't create dir seq\n"))

    if (!dir.create("results", mode="750"))
        return(err(FALSE, "Couldn't create dir results\n"))
    
    return(TRUE)
}


init <- function(phenotype) {
    if (dir.exists(phenotype))
        stop(sprintf("Director for %s already exists\n", phenotype))

    if (!is_valid_dirname(phenotype))
        stop(sprintf("Invalid direcotry name\n"))

    if (!make_analysis_directory(phenotype))
        stop(sprintf("couldn't make project\n"))


    ## copy scripts
    scripts_path <- system.files("scripts", package = "hwas")
    for (w in list.files(scripts_path, "*.R"))
        file.copy(file.path(scripts_path, w), ".",
                  copy.mode=TRUE,
                  overwrite=FALSE)


}
