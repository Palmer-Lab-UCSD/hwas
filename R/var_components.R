



# @brief REML log likelihood
reml_delta_normal_eqn <- function(delta, ysq, eigvals, df) {
    s <- delta + eigvals
    denom <- sum(ysq / s)
    numerator <- sum(ysq / (s * s))
    inv_sum <- sum(1/s)
    return(df * numerator / denom - inv_sum)
}


find_interval <- function(ysq, eigenvals, df) {
    delta_vals <- c(0, 10^seq(-8, 8, length.out=100))

    idx <- -1
    num_maxima <- 0
    num_minima <- 0

    prev_deriv <- reml_delta_normal_eqn(delta_vals[1],
                                        ysq,
                                        eigenvals,
                                        df)
    for (i in seq(2, length(delta_vals))) {
        deriv <- reml_delta_normal_eqn(delta_vals[i],
                                       ysq,
                                       eigenvals,
                                       df)

        if (is.na(deriv)) {
            print("reml_delta_normal_eqn returned NA for the following value")
            print(delta_vals[i])
            break
        }

        if (deriv < 0 && prev_deriv > 0) {
            idx <- i
            num_maxima <- num_maxima + 1
        } else if (deriv > 0 && prev_deriv < 0) {
            num_minima <- num_minima + 1
        }

        prev_deriv <- deriv
    }

    if (num_maxima == 0)
        return(NULL)

    return(list(interval = c(delta_vals[idx-1], delta_vals[idx]),
                num_maxima = num_maxima,
                num_minima = num_minima))
}



# @brief estimate the variance componets with REML
# @description Kang et al. Genetics 2008, the EMMA paper
estimate_var_comps <- function(phenotypes,
                               grm,
                               covariates = NULL) {

    if (!is_class(grm, "grm"))
        stop("Input is not of the grm class")


    nsamps <- length(phenotypes)

    # add the intercept
    if (is.matrix(covariates))
        p_covars <- dim(covariates)[2] + 1
    else if (is.vector(covariates))
        p_covars <- 2
    else
        p_covars <- 1

    design_matrix <- matrix(1., nsamps, p_covars)
    if (p_covars > 1)
        design_matrix[,2:p_covars] <- covariates

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
    # matrix Q can be decomposed into [Q_1, Q_2] with Q_1 being the 
    # orthonormal basis of the t dimensional column space of X and
    # Q_2 being n-t orthonormal basis vectors that span Res.  Consequently,
    #
    # S = Q_2 (Q_2^T Q_2)^{-1} Q_2^T
    #   = Q_2 Q_2^T
    #
    # no inverse or squaring necessary!
    #

    qr_decomp <- qr(design_matrix)
    qmat <- qr.Q(qr_decomp, complete=TRUE)[, (qr_decomp$rank+1):nsamps]

    eig_reml <- eigen(t(qmat) %*% grm %*% qmat, TRUE)
    reml_transformation = qmat %*% eig_reml$vectors

    transformed_y <- as.vector(t(reml_transformation) %*% phenotypes)
    transformed_y_sq <- transformed_y * transformed_y

    df <- nsamps - qr_decomp$rank

    # use grid to find region of root
    out <- find_interval(transformed_y_sq, eig_reml$values, df)

    if (is.null(out))
        return(NULL)

    fout <- uniroot(reml_delta_normal_eqn, 
                    out$interval,
                    transformed_y_sq,
                    eig_reml$values,
                    df)
    
    var_g <- sum(transformed_y_sq / (eig_reml$values + fout$root))
    var_e <- var_g * fout$root

    return(structure(list(herit = 1 / (1 + fout$root),
                          varg = var_g, 
                          vare = var_e,
                          fout = fout),
        class = "var_components"))
}
