

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
# @param phenotypes Is an n sample numeric vector of trait 
#   measurements
# @param grm the genetic relationship matrix rows and column 
#   names should be sample names
# @param covariates an n sample by p covariate matrix where row
#   names are sample id
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
                             covariates = NULL) {

    if (!is_grm(grm))
        stop("Not a valid grm")
    if (!is.vector(phenotypes))
        stop("Phenotypes must be a vector.")

    samp_names <- names(phenotypes)
    if (is.null(samp_names))
        stop("Phenotype elements must be labeled with sample name")

    if (!identical(samp_names, sample_names(bconn)))
        stop("Phenotype sample names not aligned with bcf file.")

    if (!identical(samp_names, rownames(grm)))
        stop("Phenotype sample names not aligned with grm rownames.")

    if (!identical(samp_names, colnames(grm)))
        stop("Phenotype sample names not aligned with grm colnames.")


    # note we get sample vector length checks with these validation
    if (!is.null(covariates) && !identical(samp_names, rownames(covariates)))
        stop("Phenotype sample names not aligned with covariates")

    vpars <- est_reml_vars(phenotypes, 
                           grm,
                           covariates)

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
                             vpars));
}


assoc <- function(gdata,
                  phenotypes,
                  grm,
                  covariates = NULL) {

    if (!is_grm(grm)) 
        stop("Not a valid grm")

    if (!is.matrix(gdata))
        stop("genetic data is not a matrix")

    samp_names <- names(phenotypes)

    if (!all(samp_names == rownames(gdata)))
        stop("Phenotype sample names not aligned with genetic data matrix")

    if (!all(samp_names == rownames(grm)))
        stop("Phenotype sample names not aligned with loco_grm")

    if (!all(samp_names == colnames(grm)))
        stop("Phenotype sample names not aligned with loco_grm")


    # note we get sample vector length checks with these validation
    if (is.null(covariates)) {
        vpars <- est_reml_vars(phenotypes,grm)
        return(assoc_no_covars_(gdata, phenotypes, grm, vpars))
    }

    stop("not implemented yet")
    #$ if (!all(samp_names == rownames(covariates)))
    #$     stop("Phenotype sample names not aligned with covariates")

    #$ vpars <- est_reml_vars(phenotypes, 
    #$                        grm,
    #$                        covariates)
    #$ return(assoc_(gdata, phenotypes, grm, covariates, vpars))
}
