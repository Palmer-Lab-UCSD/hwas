



# Search interval data structure
NewIntervalResult <- function(lower, upper=NULL) {
    if (is.null(upper))
        return(structure(lower,
                         class="interval",
                         interval=FALSE))

    return(structure(c(lower, upper), 
                     class="interval",
                     interval=TRUE))
}

is_interval <- function(x) {
    if (!is_class(x, "interval"))
        return(FALSE)

    return(attr(x, "interval"))
}

all_reml_delta_normal_descending <- function(ysq, eigvals, df, delta_grid) {

    for (delta in delta_grid) {
        deriv <- reml_delta_normal_eqn(delta, ysq, eigvals, df)
        if (deriv >= 0)
            return(FALSE)
    }

    return(TRUE)
}


# @brief use grid based search to define interval for root finding.
# @description Our REML estimate of the \delta parameter is performed
#   by finding the root of its normal equation, an equation implemented 
#   by the reml_delta_normal_eqn function.  The challenge is that
#   delta exists on the inteval of [0,\infinity) and the uniroot
#   function in R requires a search interval to be specified.  This
#   function defines a grid of delta values across several orders of 
#   magnitude.  For each value of delta, the normal equation is 
#   computed and a search interval defined by the first interval pair
#   in which delta_i-1 > 0 and delta_i < 0.
#
#   The case where no values of delta meet this criteria is possible.
#   For example, if the optimal value of delta is greater than any
#   delta in the grid.  To see how to procede in such a circumstance
#   recall the definition of 
#
#   delta = variance(environment) / variance(genetics)
#   
#   which diverges as variance(genetics) -> 0.  If the phenotypic 
#   variance attributed to genetics goes to zero, so does the 
#   heritability and the concept of delta meaningless.  In such
#   a case we expect that the derivative of the REML likelihood 
#   to be positive for all values of delta.
# @param ysq is a vector of length n samples of squared phenotypes.
# @param eigenvals are the n eigenvalues of the grm
# @param df is the degress of freedom defined as 
#   n samples - dim(colspace(X covariates))
# @return the value delta, an interval, NULL. Cases are enumerated below,
#   for simplicity let d be the derivative of REML likelihood implemented
#   in reml_delta_normal_eqn.
#   0   represents heritability = 1, no enviroment variance
#       d < 0 for all values of delta in grid
#       d == 0 for delta = 0 and d < 0 for remaining delta
#   Inf represents heritability = 0, no genetic variance
#       d > 0 for all values of delta in grid
#   c(delta lower bound, delta upper bound) d > 0 for all delta <= 
#       delta lower bound and d < 0 for all delta >= delta upper
#       bound
#   NULL for complete failure
find_interval <- function(ysq, eigenvals, df) {
    ngrid_points <- 100
    delta_grid = 10^seq(-8, 8, length.out=ngrid_points)

    delta_lower_bound <- 0
    deriv_lower_bound <- reml_delta_normal_eqn(delta_lower_bound,
                                               ysq,
                                               eigenvals,
                                               df)
    if (is.null(deriv_lower_bound) 
        || is.na(deriv_lower_bound) 
        || is.nan(deriv_lower_bound))
        return(NULL)

    if (deriv_lower_bound <= 0) {
        if (all_reml_delta_normal_descending(ysq, eigenvals, df, delta_grid)) {
            return(NewIntervalResult(0))
        } else {
            return(NULL)
        }
    }

    for (delta_upper_bound in delta_grid) {

        deriv_upper_bound <- reml_delta_normal_eqn(delta_upper_bound,
                                                   ysq,
                                                   eigenvals,
                                                   df)
        # something went wrong, the returned value should be a number
        if (is.null(deriv_lower_bound) 
            || is.na(deriv_lower_bound) 
            || is.nan(deriv_lower_bound))
            return(NULL)
    
        if (deriv_upper_bound < 0)
            return(NewIntervalResult(delta_lower_bound, 
                                     delta_upper_bound))


        # note, if the deriv_upper_bound is zero then we want to still
        # check that the derivative at the next value of delta is negative
        # otherwise the returned inteval would have an extremum at the
        # interval boundary, which we can't have.
        if (deriv_upper_bound != 0) {
            delta_lower_bound <- delta_upper_bound
            deriv_lower_bound <- deriv_upper_bound
        }
    }

    # case 1: derivative pattern +[0+]*0
    if (deriv_upper_bound == 0) {
        delta_upper_bound <- delta_upper_bound + 100
        deriv_upper_bound <- reml_delta_normal_eqn(delta_upper_bound,
                                                   ysq,
                                                   eigenvals,
                                                   df)
        if (deriv_upper_bound < 0)
            return(NewIntervalResult(delta_lower_bound, delta_upper_bound))
    }

    # case 2
    # if deriv_upper_bound >= 0 for all grid points then assume delta->Inf,
    # heritability is zero
    return(NewIntervalResult(Inf))
}


# @brief REML log likelihood
reml_delta_normal_eqn <- function(delta, ysq, eigvals, df) {
    s <- delta + eigvals
    denom <- sum(ysq / s)
    numerator <- sum(ysq / (s * s))
    inv_sum <- sum(1 / s)
    return(df * numerator / denom - inv_sum)
}



# @brief Internal function call to estimate variance component pars
# @description Estimate REML variance component from data assumed
#   to be correctly formated.  This should only be called by internal
#   functions.
reml_vars_computation <- function(phenotypes,
                                  grm,
                                  design_matrix) {

    # Use QR decomposition for the matrix that projects the response
    # variable Y to the residual vector space.  The data are defined
    # by a n dim sample vector space, V. The the columns of data 
    # covariate matrix spans a t dimensional subspace of V, I define
    # the residual vector space (Res) as the subspace of V that complements
    # the column space of X.  Let S be the projection matrix from 
    # vector space V to Res
    #
    # S = I - X(X^T X)^{-1} X^T
    #
    # Under the complete QR decomposition of X the n x n orthonormal 
    # matrix Q can be partitioned into [Q_1, Q_2] with Q_1 being an
    # orthonormal basis of the t dimensional column space of X and
    # Q_2 being n-t orthonormal basis vectors of the complementary
    # subspace.
    #
    # S = Q_2 Q_2^T
    #
    # no inverse or squaring necessary!
    nsamps <- dim(grm)[1]
    qr_decomp <- qr(design_matrix)
    qmat <- qr.Q(qr_decomp, complete=TRUE)[, (qr_decomp$rank+1):nsamps]

    # Add identity matrix to grm for eigenvector computation.  This will
    # help eigenvector solution in the case that the grm has rank < N.
    # I remove the identity be subtracting 1 from the eigen values
    # below.
    eig_reml <- eigen(t(qmat) %*% (grm + diag(1, nsamps)) %*% qmat, 
                      TRUE)
    eig_reml$values <- eig_reml$values - 1
    reml_transformation = qmat %*% eig_reml$vectors

    transformed_y <- as.vector(t(reml_transformation) %*% phenotypes)
    transformed_y_sq <- transformed_y * transformed_y

    df <- nsamps - qr_decomp$rank

    # use grid to find region of root
    out <- find_interval(transformed_y_sq, eig_reml$values, df)
    if (is.null(out))
        return(NULL)

    if (is_interval(out)) {
        fout <- uniroot(reml_delta_normal_eqn, 
                        out,
                        transformed_y_sq,
                        eig_reml$values,
                        df)

        ## check for convergence
        delta = fout$root
        var_g <- sum(transformed_y_sq / (eig_reml$values + delta))
        var_e <- var_g * delta
        herit <- 1 / (1 + delta)

    } else if (is.infinite(out)) {
        # Under this case var_g is zero and consequently covariance matrix
        # is independent of the GRM.  Under the model the remaining source of
        # variation is the environment cov(Y) = var_e I.  I use the n x n-p
        # Q matrix from QR decomposition of the n x p covariates matrix X, 
        # intercept included, as REML P in Patterson and Thompson, Biometrika
        # 1971.  We see that Q satisfies the requirements of Patterson and 
        # Thompson
        #   1) Q^T Q = I by the orthonormal properties of the QR decomp
        #   2) Q Q^T = S, I haven't figure out how to prove this.  I do know 
        #       by (1) that Q Q^T is idempotent as is required by S.  We do
        #       know that rank(Q Q^T) = n-k, which is also the rank of S.  
        #       Lastly, I know the eigen values of Q Q^T and S are equal as
        #       both are idempotent matrices.  Idempotent matrics have eigenvalues 
        #       {0,1} and consequently their sum is the rank.  Given that both
        #       matrices, QQ^T and S, have equal rank they must have equal 
        #       eigenvalues.
        yreml <- as.vector(t(qmat) %*% phenotypes)
        var_e <- sum(yreml * yreml) / df
        var_g <- 0
        herit <- 0
        delta <- Inf
        fout <- NULL
    } else if (out == 0) {
        ## Check the rank of t(qmat) %*% grm %*% qmat to determine if any
        ## eigenvalues should be interpreted as zero
        ## TODO: double check that I am accounting for rank deficient GRM 
        ## correctly
        qr_tmp <- qr(eig_reml$vectors)
        if (qr_tmp$rank == df) {
            var_g <- sum(transformed_y_sq / eig_reml$values) / df
        } else {
            idx <- seq(1, qr_tmp$rank)
            var_g <- sum(transformed_y_sq[idx] / eig_reml$values[idx]) / qr_tmp$rank
        }

        delta <- 0
        var_e <- 0
        herit <- 1
        fout <- NULL
    }

    return(structure(list(herit = herit,
                          vare = var_e,
                          varg = var_g, 
                          delta = delta,
                          fout = fout),
        class = "var_components"))
}


# @brief estimate the variance componets with REML
# @description 
# @param phenotypes vector of phenotype measurements, each element
#   must be labeled with sample name
# @param grm is a symmetric positive semi-definite matrix that satisfies
#   the normalization requirement Kang et al. Nature Genetics 2010
# @param covariates NULL or an n sample by p covariate matrix where row
#   names are sample id.  If NULL, the interecept must be TRUE.
# @param intercept is a boolean on whether an intercept should be modeled
# @return S3 class of var_components, data are in list with attributes:
#   herit: estimated narrow sense heritability
#   vare: the estimated environment variance
#   varg: the estimated variance of polygenic effects
#   delta: value optimized by the EMMA algorithm for REML estimation of
#       variance components, Kang et al. 2008 Genetics.
#   fout: NULL or output of R's uniroot function
est_reml_vars <- function(phenotypes,
                          grm,
                          covariates = NULL,
                          intercept = TRUE) {

    # Validate user input
    if (!is.vector(phenotypes))
        stop("Phenotypes must be a vector.")

    snames <- names(phenotypes)
    if (is.null(snames))
        stop("Phenotype elements must be labeled with sample name")

    if (!is_grm(grm))
        stop("Input is not a valid the grm class")

    # note we get sample vector length checks with these validation
    if (!identical(snames, rownames(grm)))
        stop("Phenotype sample names not aligned with grm rownames")

    if (!identical(snames, colnames(grm)))
        stop("Phenotype sample names not aligned with grm colnames")

    # I prefer the code does what the user expects, as such I don't
    # want to set intercept to TRUE when covariates is NULL in case
    # the user was not expecting this.
    if (is.null(covariates) && !intercept)
        stop("If no covariates provided intercept must be TRUE")

    if (!is.null(covariates) && !identical(snames, rownames(covariates)))
        stop("Phenotype sample names not aligned with covariates")

    nsamps <- length(phenotypes)

    # add the intercept
    if (intercept) {
        design_matrix <- matrix(1., nrow = nsamps, ncol = 1)

        if (is.matrix(covariates) || is.vector(covariates))
            design_matrix <- cbind(design_matrix, covariates)
    } else
        design_matrix <- covariates

    return(reml_vars_computation(phenotypes, grm, design_matrix))
}
