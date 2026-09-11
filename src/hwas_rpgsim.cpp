#include <Rpgsim.h>


const char* rpgsim::status_msg(rpgsim::Status status) {
    switch (status) {
    case rpgsim::Status::Success:
        return "Success";
    case rpgsim::Status::ErrNotSymmetricMatrix:
        return "Not a symmetric matrix";
    case rpgsim::Status::ErrNotSquareMatrix:
        return "Not a square matrix";
    case rpgsim::Status::ErrHeritabilityOutOfRange:
        return "Heritability must be on interval (0,1).";
    default:
        break;
    }

    return "Unexpected status, please contact maintainers.";
}


bool is_symm_(const Rcpp::NumericMatrix& A) {
    double tol = 1e-10;

    int nrow = A.nrow();
    if (nrow != A.ncol())
        return false;

    double upper_bound = 0;
    double lower_bound = 0;

    for (int i = 0; i < nrow; i++) {
        for (int j = i+1; j < nrow; j++) {

            lower_bound = A(j, i) - tol;
            upper_bound = A(j, i) + tol;

            if (A(i, j) < lower_bound || A(i, j) > upper_bound)
                return false;
        }
    }

    return true;
}

rpgsim::Status matrix_sum(const Rcpp::NumericMatrix& A, double& val) {
    double s = 0;
    for (int i = 0; i < A.nrow(); i++) {
        for (int j = 0; j < A.ncol(); j++)
            s += A(i, j);
    }

    val = s;

    return rpgsim::Status::Success;
}


rpgsim::Status matrix_trace(const Rcpp::NumericMatrix& A, double& val) {
    if (A.nrow() != A.ncol())
        return rpgsim::Status::ErrNotSquareMatrix;

    double mtrace = 0;
    for (int i = 0; i < A.nrow(); i++)
        mtrace += A(i, i);

    val = mtrace;
    return rpgsim::Status::Success;
}


rpgsim::Status calc_pgsim_pars_(const Rcpp::NumericMatrix& grmatrix,
        PgSimParams& params) {

    // n represents number of samples
    int n = grmatrix.nrow();
    double ndouble = static_cast<double>(n);
    
    double grm_trace = 0;
    rpgsim::Status status = matrix_trace(grmatrix, grm_trace);
    if (status != rpgsim::Status::Success)
        return status;

    double grm_sum = 0;
    status = matrix_sum(grmatrix, grm_sum);
    if (status != rpgsim::Status::Success)
        return status;

    double numerator = params.heritability * static_cast<double>(n * (n - 1));
    double denom = static_cast<double>(1) - params.heritability;
    denom = denom * (ndouble * grm_trace - grm_sum);

    params.var_e = 1;
    params.var_u = numerator / denom;
    return rpgsim::Status::Success;
}


// Not public facing and consquently assume that bcf_filename
// is a bcf file, heritability is valid
// int estimate_pars(bcfio::Bcf* input_bid,
//         float heritability, 
//         PolygenicParameters* pars) {
// 
//     // bconn_t local_bid = bopen(input_bid->fname_.c_str(), "r");
//     
//     printf("NOT IMPLEMENTED\n");
//     return 1;
// }


// [[Rcpp::export]]
Rcpp::NumericVector pg_sim(const Rcpp::NumericMatrix& grmatrix,
        const float heritability) {

    if (heritability <= 0 || heritability >=1)
        Rcpp::stop("Heritability, h, must be 0 < h < 1");

    if (!is_symm_(grmatrix))
        Rcpp::stop(rpgsim::status_msg(rpgsim::Status::ErrNotSymmetricMatrix));

    PgSimParams params {};
    params.heritability = heritability;

    rpgsim::Status status = calc_pgsim_pars_(grmatrix, params);
    if (status != rpgsim::Status::Success)
        Rcpp::stop(rpgsim::status_msg(status));

    // Constructing covariance matrix from grm and model parameters
    // Recall that:
    //
    // covmatrix = GRM * var_u + I * var_e
    int nsamples = grmatrix.nrow();
    Eigen::MatrixXd covmat = Eigen::MatrixXd::Zero(nsamples, nsamples);
    for (int i = 0; i < nsamples; i++) {
        for (int j = 0; j < nsamples; j++) {
            // remember that environment variance is only added to diagonal
            if (i == j)
                covmat(i, i) = params.var_e + params.var_u * grmatrix(i, i);
            else
                covmat(i, j) = params.var_u * grmatrix(i, j);
        }
    }

    // LU decomposition
    // Recall that by construction the covmatrix, additional of diagonal
    // environment variance, it is guaranteed to be symmetric and positive
    // positive definite, therefore we can use Cholesky without issue
    
    Eigen::LLT<Eigen::MatrixXd> choleskyL(covmat);

    // Sampling
    //
    // We need to sample from the multivariate normal.  This is simple
    // to do using the Cholesky decompostion of the covariance matrix.
    // The recipe is to:
    //  1) independently sample each element of the N x = sample phenotype 
    //      column vector P from the standard normal, N(0,1).
    //  2) generate the correlated phenotype values by L * P, with L 
    //      being the cholesky decomposition of the covariance matrix.
    //
    // Note on random seed, using R::rnorm under the default 
    // [[Rcpp::export]] statement we inherit the random state of the R program.  
    // Therefore, the seed used for random numbers can be set in R with 
    // set.seed(<number>) and used here.
    //
    // Note: VectorXd is a column vector
    Eigen::VectorXd phenotypes(nsamples);
    for (int i = 0; i < nsamples; i++)
        phenotypes[i] = R::rnorm(0, 1);


    // there probably is a better way to convert the type Eigen::VectorXd to
    // Rcpp::NumericVector, but I can't figure it out, explicity should be
    // sufficient for now.
    phenotypes = choleskyL.matrixL() * phenotypes;
    Rcpp::NumericVector out(nsamples);
    for (int i = 0; i < nsamples; i++)
        out[i] = phenotypes[i];

    out.attr("var_u") = params.var_u;
    out.attr("var_e") = params.var_e;

    // decomposition grm 
    return out;
}


// [[Rcpp::export]]
Rcpp::NumericVector pg_sim_qtl(const bconn_t bconn,
        const Rcpp::NumericMatrix& grmatrix,
        const float qtl_freq,
        const float qtl_effect_size,
        const float heritability,
        const char* id,
        const int seed) {

    Rcpp::stop("not yet implemented");
    // if (!is_open(bconn))
    //     Rcpp::stop(bcfio::status_msg(bcfio::Status::ErrBcfNotOpen));

    // if (heritability <= 0 || heritability >= 1)
    //     Rcpp::stop(rpgsim::status_msg(rpgsim::Status::ErrHeritabilityOutOfRange));

    // uint32_t nsamps = num_samples(bconn);
    // bcfio::bid_t bid = bcfio::replicate(bconn.get());
    // if (!bid)
    //     Rcpp::stop("Internal data structure error");

    // PgSimParams pars { 1, 1, qtl_freq, qtl_effect_size, heritability }; 
    // rpgsim::Status status = calc_pgsim_pars_(grmatrix, pars);
    // if (status == rpgsim::Status::Success)
    //     Rcpp::stop(rpgsim::status_msg(status));
}
