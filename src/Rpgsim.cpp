#include <Rpgsim.h>


Status calc_pgsim_pars_(const Rcpp::NumericMatrix& grmatrix,
        PgSimParams& params) {

    // n represents number of samples
    int n = grmatrix.nrow();
    double ndouble = static_cast<double>(n);
    
    double grm_trace = 0;
    Status status = matutils::matrix_trace(grmatrix, grm_trace);
    if (status != Status::Success)
        return status;

    double grm_sum = 0;
    status = matutils::matrix_sum(grmatrix, grm_sum);
    if (status != Status::Success)
        return status;

    double numerator = params.heritability * static_cast<double>(n * (n - 1));
    double denom = static_cast<double>(1) - params.heritability;
    denom = denom * (ndouble * grm_trace - grm_sum);

    params.var_e = 1;
    params.var_u = numerator / denom;
    return Status::Success;
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

    if (!matutils::is_symm_(grmatrix))
        Rcpp::stop(status_msg(Status::ErrNotSymmetricMatrix));

    PgSimParams params {};
    params.heritability = heritability;

    Status status = calc_pgsim_pars_(grmatrix, params);
    if (status != Status::Success)
        Rcpp::stop(status_msg(status));

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
    //     Rcpp::stop(status_msg(bcfio::Status::ErrBcfNotOpen));

    // if (heritability <= 0 || heritability >= 1)
    //     Rcpp::stop(status_msg(Status::ErrHeritabilityOutOfRange));

    // uint32_t nsamps = num_samples(bconn);
    // bcfio::bid_t bid = bcfio::replicate(bconn.get());
    // if (!bid)
    //     Rcpp::stop("Internal data structure error");

    // PgSimParams pars { 1, 1, qtl_freq, qtl_effect_size, heritability }; 
    // Status status = calc_pgsim_pars_(grmatrix, pars);
    // if (status == Status::Success)
    //     Rcpp::stop(status_msg(status));
}
