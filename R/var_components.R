# REML estimation of variance components in a linear mixed model
#
# NOTE: This is a refactor of my original implementation using
# recommendations of Claude Opus 5.5.  
#
# The aim is to implement the variance component REML estimates
# described by Kang et al. 2008 Genetics, the so-called EMMA 
# method.  The entry point for variance component estimation for
# users is the "est_reml_vars" function as it validates user
# input and then calls the internal implementation function
# reml_vars_.
# 
# The reml_vars_ function estimates the REML variance components
# of the model
#
# y = X b + g + e
# g ~ N(0, var_g A)
# e ~ N(0, var_e I)
# 
# define V = var_g A + vbar_e I
#
# by strategy of Kang et al..  This strategy is an implementation
# of that defined by Patterson and Thompson 1971 Biometrika that
# separates the log likelihood into two components that span
# orthogonal subspaces of the data space Omega: 
#
#   1. the subspace Omega_x in Omega that spans colspace(X),
#   2. the subspace Omega_v in Omega that is the complement 
#       of Omega_x.
# 
# The REML component is the log-likelihood over Omega_v and
# is determined by the transformation P' y, where P' is the
# transpose of the (N x (N-p)) matrix P.  Kang et al. show
# that the matrix of N-p eigenvectors of the transformed 
# covariance matrix SVS, with S = I - X(X' X)^{-1} X', 
# satisfies the requirements for P under Patterson and Thompson. 
# 
# Maximizing the REML component for estimating var_g and var_e
# is implemented using the following definitions
#
# 1. eta = P' y and eta_sq as the square of each element.
# 2. df = N - p where p is the rank of fixed effects 
#       matrix X.
# 3. delta = var_e / var_g
# 4. lambda to be the df eigenvalues of P' V P
# 5. w to be df length vector where w_i = eta_i / (lamba_i + delta)
# 6. u to be df length vecot where u_i = eta_i / sqrt(lambda_i + delta)
#
# Under these definitions the REML log likelihood (Eq 7 of Kang et al.)
#
#   reml_loglik = df/2 [ 
#       log(df / (2 pi)) - 1 
#       - log(u' u) 
#       - (1/df) sum(log(lambda + delta))
#   ]
#
# and its derivative w.r.t. delta (Eq 9 of Kang et al.) 
#
#   dreml_loglik = df/2 (w' w) / (u' u) 
#           - 1/2 sum(1 / (lambda + delta).
#
# A maximum is found by:
#
# 1. Define grid of delta = [0, 
# 10^seq(-8, 8, length.out=ngrid_points), Inf]
# 2. At each interior point determine whether a root exists
# 3. Compute reml_loglik at boundaries and roots.
# 4. 


# The REML log-likelihood as a function of delta.  The genetic
# variance parameter has been profiled out.
reml_loglik_ <- function(delta, eta_sq, lambda, df) {
    const_vals <- -0.5 * df * log(2*pi) - 0.5 * df

    ## delta is infinite then there is no genetic component,
    ## the reml_loglikelihood is then in terms of environment 
    ## variance, var_e I.
    if (is.infinite(delta))
        return(const_vals - 0.5 * df * log(sum(eta_sq) / df))

    l <- lambda + delta
    if (any(l <= 0))
        return(-Inf)

    const_vals <- const_vals + 0.5 * df * log(df)
    return(const_vals - 0.5 * df * log(sum(eta_sq / l)) - 0.5*sum(log(l)))
}


# caller is responsible for checking that delta + lambda > 0
dreml_loglik_ <- function(delta, eta_sq, lambda, df) {
    l <- lambda + delta
    w <- sum(eta_sq / (l * l))
    u <- sum(eta_sq / l)
    return(0.5 * df * w / u - 0.5 * sum(1 / l))
}


# NOTE: This function written by Claude Opus 5.5
# @brief validate and clean eigenvalues of Q2^T K Q2
# @description Eigenvalues within tol of zero are set to exactly zero.  An
#   eigenvalue below -tol means the GRM is not positive semi-definite beyond
#   floating point error of the eigensolver; the REML likelihood is undefined
#   on part of delta's domain and an error is raised.
# @param lambda eigenvalues
# @param tol NULL or non-negative tolerance.  NULL uses the conventional
#   rank tolerance length(lambda) * eps * max|lambda|.
clean_eigenvalues <- function(lambda, tol = NULL) {
    if (is.null(tol))
        tol <- length(lambda) * .Machine$double.eps * max(abs(lambda))

    if (min(lambda) < -tol)
        stop(sprintf(paste0("GRM (projected onto the residual space) is not ",
                            "positive semi-definite: min eigenvalue %.3e < -tol ",
                            "(%.3e). Check the precision used to accumulate ",
                            "the GRM."),
                     min(lambda), tol),
             call. = FALSE)

    lambda[lambda < tol] <- 0
    return(lambda)
}



# NOTE: This function written by Claude Opus 5.5
# @brief maximize the REML log likelihood over delta in [0, Inf]
# @description EMMA strategy: evaluate dl/d delta on a grid in log(delta),
#   refine every + to - sign change (local maxima) with uniroot on the
#   log scale, then select the global maximum among the interior local
#   maxima and the two boundaries delta = 0 (h^2 = 1) and delta = Inf
#   (h^2 = 0).
# @return list(delta, loglik, fout) with fout the uniroot result for the
#   chosen interior root, or NULL for a boundary solution.
optimize_reml_delta <- function(eta_sq,
                                lambda,
                                df,
                                log_delta_range = c(-10, 10),
                                ngrid = 100,
                                tol = 1e-10) {

    log_grid <- seq(log_delta_range[1], log_delta_range[2],
                    length.out = ngrid + 1)
    dl <- vapply(exp(log_grid), dreml_loglik_, numeric(1),
                 eta_sq = eta_sq, lambda = lambda, df = df)
    if (anyNA(dl))
        stop("Non-finite REML derivative on the search grid.", call. = FALSE)

    # boundary candidates
    cand_delta <- c(0, Inf)
    cand_ll <- c(reml_loglik_(0, eta_sq, lambda, df),
                 reml_loglik_(Inf, eta_sq, lambda, df))
    cand_fout <- list(NULL, NULL)

    # interior local maxima: derivative changes sign from + to -
    f <- function(log_delta)
        dreml_loglik_(exp(log_delta), eta_sq, lambda, df)

    idx <- which(dl[-length(dl)] > 0 & dl[-1] <= 0)
    for (i in idx) {
        fout <- uniroot(f,
                        lower = log_grid[i],
                        upper = log_grid[i + 1],
                        f.lower = dl[i],
                        f.upper = dl[i + 1],
                        tol = tol)
        delta <- exp(fout$root)
        cand_delta <- c(cand_delta, delta)
        cand_ll <- c(cand_ll, reml_loglik_(delta, eta_sq, lambda, df))
        cand_fout <- c(cand_fout, list(fout))
    }

    best <- which.max(cand_ll)
    return(list(delta = cand_delta[best],
                loglik = cand_ll[best],
                fout = cand_fout[[best]]))
}

# @brief Internal function call to estimate variance component pars
# @description Estimate REML variance component from data assumed
#   to be correctly formated.  This should only be called by internal
#   functions.
# @param phenotypes is a numeric vector of N valid phenotypes
# @param grm is a numeric N x N normalized genetic similarity matrix
# @param design_matrix is an N x p matrix of fixed effects, this may 
#   include an intercept term.
reml_vars_ <- function(phenotypes, grm, design_matrix) {

    nsamps <- dim(grm)[1]
    qr_decomp <- qr(design_matrix)

    df <- nsamps - qr_decomp$rank
    if (df < 2)
        stop("Fewer than 2 residual degrees of freedom.")

    qmat <- qr.Q(qr_decomp, complete=TRUE)[, (qr_decomp$rank+1):nsamps]

    eig_reml <- eigen(t(qmat) %*% grm %*% qmat, TRUE)
    lambda <- clean_eigenvalues(eig_reml$values)

    reml_transformation = qmat %*% eig_reml$vectors

    transformed_y <- as.vector(t(reml_transformation) %*% phenotypes)
    transformed_y_sq <- transformed_y * transformed_y

    if (sum(transformed_y_sq) == 0)
        stop("Phenotypes lie in the column space of the design matrix.",
             call. = FALSE)

    # use grid to find region of root
    opt <- optimize_reml_delta(transformed_y_sq, lambda, df)
    if (is.infinite(opt$delta)) {
        var_g <- 0
        var_e <- sum(transformed_y_sq) / df
        herit <- 0
    } else {
        # Kang et al. 2008: s2g = sum(eta^2 / (lambda + delta)) / df.
        # delta == 0 is only selected when every lambda > 0 (otherwise
        # loglik(0) == -Inf), so the division is safe.
        var_g <- sum(transformed_y_sq / (lambda + opt$delta)) / df
        var_e <- opt$delta * var_g
        herit <- 1 / (1 + opt$delta)
    }

    return(structure(list(herit = herit,
                          vare = var_e,
                          varg = var_g, 
                          delta = opt$delta,
                          loglik = opt$loglik,
                          fout = opt$fout),
            class = "var_components"))
}


# @brief Estimate the variance componets by REML
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

    nsamps <- length(phenotypes)

    # add the intercept
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

    return(reml_vars_(phenotypes, grm, design_matrix))
}


#   Use QR decomposition for the matrix that projects the response
#   variable Y to the residual vector space.  The data are defined
#   by a n dim sample vector space, V. The the columns of data 
#   covariate matrix spans a t dimensional subspace of V, I define
#   the residual vector space (Res) as the subspace of V that complements
#   the column space of X.  Let S be the projection matrix from 
#   vector space V to Res
#
#   S = I - X(X^T X)^{-1} X^T
#
#   Under the complete QR decomposition of X the n x n orthonormal 
#   matrix Q can be partitioned into [Q_1, Q_2] with Q_1 being an
#   orthonormal basis of the t dimensional column space of X and
#   Q_2 being n-t orthonormal basis vectors of the complementary
#   subspace.
#
#   S = Q_2 Q_2^T
#
#   no inverse or squaring necessary!
