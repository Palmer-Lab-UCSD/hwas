# Data validation tests
#
# Implements a simple model for data validation.  Each of the 
# is_valid_<data_type> function tests whether the user inputs
# satisfies the requirements of <data_type>.  The return value
# of the is_valid_<data_type> functions is an instance of the
# IsValid S3 class.  This takes boolean values and also a 
# message attribute.  The message to print to users can be
# contstucted from the invalid_data_msg function.


NewIsValid <- function(is_valid, msg=NULL) {
    return(structure(is_valid,
                     class="IsValid",
                     msg = msg))
}


invalid_data_msg <- function(is_valid) {
    return(paste("INVALID DATA", attr(is_valid, "msg")))
}


is_valid_phenotypes <- function(y) {
    if (!is.vector(y, mode="numeric"))
        return(NewIsValid(FALSE,
                     "Phenotype must be a numeric vector"))

    if (any(is.nan(y)))
        return(NewIsValid(FALSE,
                     "Phenotype has NaN values"))

    if (any(is.na(y)))
        return(NewIsValid(FALSE,
                     "Phenotype has missing, NA, values."))
    
    if (any(is.infinite(y)))
        return(NewIsValid(FALSE,
                     "Phenotype has infinite values"))

    if (is.null(names(y)))
        return(NewIsValid(FALSE,
                      "Phenotype elements are not named"))

    return(NewIsValid(TRUE))
}


# @brief Determine whether input grm is valid grm
# @description A valid grm must satisfy the normalization requirement
#   of Kang et al. 2010 Nature Genetics (EMMAX procedure) and have
#   rows / col names labeled in order of the provided sample names.
# @param grm is an instance of the grm class
# @param samp_names is a character vector of sample names
# @param tol is a tolerance for the grm normalization test
# @return bool stored as an IsValid S3 class
is_valid_grm <- function(grm, samp_names, tol=1e-8) {

    if (!identical(snames, rownames(grm)))
        return(NewIsValid(FALSE, 
                          paste("GRM rownames are not aligned with",
                                "sample names.")))

    if (!identical(snames, colnames(grm)))
        return(NewIsValid(FALSE, 
                          paste("GRM colnames are not aligned with",
                                "sample names.")))

    if (any(is.nan(grm)))
        return(NewIsValid(FALSE,
                          "GRM has NaN values"))

    if (any(is.na(grm)))
        return(NewIsValid(FALSE,
                          "GRM has missing, NA, values"))

    if (!is_grm_normalized(grm, tol))
        return(NewIsValid(FALSE,
                          paste("GRM is not correctly normalized,",
                            "use hwas::gsim_to_grm to normalize")))

    return(NewIsValid(TRUE))
}


is_valid_covars <- function(covars, samp_names) {
    if (is.null(covars))
        return(NewIsValid(TRUE))

    if (!is.vector(covars, mode="numeric") && !is.matrix(covars))
        return(NewIsValid(FALSE,
                          paste("Invalid covariate data structure, it",
                                "must be one of the following: ",
                                "NULL, numeric vector, or matrix")))

    if (is.vector(covars) && !identical(names(covars), samp_names))
        return(NewIsValid(FALSE,
                          paste("The labels of covariate vector",
                                "do not match input sample names")))

    if (is.matrix(covars) & !identical(rownames(covars), samp_names))
        return(NewIsValid(FALSE,
                          paste("The rownames of covariate matrix",
                                "do not match input sample names")))

    return(NewIsValid(TRUE))
}


is_valid_bconn <- function(bconn, samp_names) {
    if (!is_open(bconn))
        return(NewIsValid(FALSE,
                          "Bcf file is not open for reading"))

    if (!identical(samp_names, sample_names(bconn)))
        return(NewIsValid(FALSE,
                          paste("Sample names configured for bcf I/O",
                                "do not match input sample names.")))

    return(NewIsValid(TRUE))
}
