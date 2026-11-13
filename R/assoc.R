

# @brief Run association analysis 
# @description The association analysis tests the hypothesis that
#
#   H_0: \beta_1 = \beta_2 = ... = \beta_k 
#   H_A: \beta_1 or \beta_2 or ... or \beta_k not equal 0
#
#   by computing the F-statistic at each genomic position.  Here k
#   is the dimension of the genetic data.  For example, if one uses
#   GP for the genotype posterior probability of a biallelic variant
#   k could be equal to 3.  Another example is the Haplotype Dose
#   (HD) where k is then the number of haplotypes.  The genomic
#   positions associations are calculated are specified in the bconn 
#   Bcf file connection object.
# @param bconn is a Bcf file connection object 
# @param format_id is the format key/id of the genetic record 
#   used in the association analysis.
# @param phenotypes is an n sample numeric vector of trait measurements
# @param grm the genetic relationship matrix rows and column 
#   names should be sample names
# @param covariates NULL or an n sample by p covariate matrix where row
#   names are sample id.  If NULL, intercept must be TRUE.
# @param intercept is a boolean on whether an intercept should be modeled
# @return List(contig_name_1 = List(fstats = numeric vector,
#                                   pos = numeric vector of genomic 
#                                         coordinates, vcf 1-based),
#              contig_name_2 = List(fstats = numeric vector,
#                                   pos = numeric vector of genomic 
#                                         coordinates, vcf 1-based))
stream_ctg_assoc <- function(bconn, 
                             format_id, 
                             phenotypes, 
                             grm, 
                             covariates = NULL,
                             intercept = TRUE) {

    # TODO: add check for format_id in set of formats stored
    # in bcf.  To do this I need to add the format key enumeration
    # tool in the Rbcfio.h file

    # Validate user input
    if (!(is_val <- is_valid_phenotypes(phenotypes)))
        stop(invalid_data_msg(is_val))

    samp_names <- names(phenotypes)
    if (!(is_val <- is_valid_bconn(bconn, samp_names)))
        stop(invalid_data_msg(is_val))

    if (!(is_val <- is_valid_grm(grm, samp_names)))
        stop(invalid_data_msg(is_val))

    if (!(is_val <- is_valid_covars(covariates, samp_names)))
        stop(invalid_data_msg(is_val))

    # I prefer the code does what the user expects, as such I don't
    # want to set intercept to TRUE when covariates is NULL in case
    # the user was not expecting this.
    if (is.null(covariates) && !intercept)
        stop("When covariates are NULL, the intercept must be TRUE")

    # construct reduced model design matrix, this will be some
    # combination of an intercept term and input covariates
    if (intercept) {
        design_matrix <- matrix(1., nrow = nsamps, ncol = 1)

        if (!is.null(covariates))
            design_matrix <- cbind(design_matrix, covariates)
    } else if (is.vector(covariates))
        design_matrix <- matrix(covariates, 
                                nrow=length(covariates),
                                ncol=1)
    else
        design_matrix <- covariates

    vpars <- reml_vars_(phenotypes, 
                        grm,
                        design_matrix)

    if (is.null(covariates)) {
        return(stream_ctg_assoc_no_covs_(bconn,
                                         format_id,
                                         phenotypes,
                                         grm,
                                         vpars))
    }

    return(stream_ctg_assoc_(bconn,
                             format_id,
                             phenotypes,
                             covariates,
                             grm,
                             vpars))
}


# @brief Compute association statistics of a single locus
assoc <- function(gdata,
                  phenotypes,
                  grm,
                  covariates = NULL,
                  intercept = TRUE) {
    stop("Needs work")

    # Validate user input
    if (!(is_val <- is_valid_phenotypes(phenotypes)))
        stop(invalid_data_msg(is_val))

    samp_names <- names(phenotypes)
    if (!(is_val <- is_valid_grm(grm, samp_names)))
        stop(invalid_data_msg(is_val))

    if (!(is_val <- is_valid_covars(covariates, samp_names)))
        stop(invalid_data_msg(is_val))

    # I prefer the code does what the user expects, as such I don't
    # want to set intercept to TRUE when covariates is NULL in case
    # the user was not expecting this.
    if (is.null(covariates) && !intercept)
        stop("When covariates are NULL, the intercept must be TRUE")

    # construct reduced model design matrix, this will be some
    # combination of an intercept term and input covariates
    if (intercept) {
        design_matrix <- matrix(1., nrow = nsamps, ncol = 1)

        if (!is.null(covariates))
            design_matrix <- cbind(design_matrix, covariates)
    } else if (is.vector(covariates))
        design_matrix <- matrix(covariates, 
                                nrow=length(covariates),
                                ncol=1)
    else
        design_matrix <- covariates

    vpars <- reml_vars_(phenotypes, 
                        grm,
                        design_matrix)



    # if (!is.matrix(gdata))
    #     stop("genetic data is not a matrix")

    # samp_names <- names(phenotypes)

    # if (!all(samp_names == rownames(gdata)))
    #     stop("Phenotype sample names not aligned with genetic data matrix")

    # if (!all(samp_names == rownames(grm)))
    #     stop("Phenotype sample names not aligned with loco_grm")

    # if (!all(samp_names == colnames(grm)))
    #     stop("Phenotype sample names not aligned with loco_grm")


    # # note we get sample vector length checks with these validation
    # if (is.null(covariates)) {
    #     vpars <- est_reml_vars(phenotypes,grm)
    #     return(assoc_no_covars_(gdata, phenotypes, grm, vpars))
    # }

    # stop("not implemented yet")
    #$ if (!all(samp_names == rownames(covariates)))
    #$     stop("Phenotype sample names not aligned with covariates")

    #$ vpars <- est_reml_vars(phenotypes, 
    #$                        grm,
    #$                        covariates)
    #$ return(assoc_(gdata, phenotypes, grm, covariates, vpars))
}
