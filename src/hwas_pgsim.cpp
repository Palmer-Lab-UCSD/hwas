
#include <Rpgsim.h>


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
        return rpgsim::Status::ErrNotSymmetricMatrix;

    double mtrace = 0;
    for (int i = 0; i < A.nrow(); i++)
        mtrace += A(i, i);

    val = mtrace;
    return rpgsim::Status::Success;
}


rpgsim::Status calc_pgsim_pars(const double heritability, 
        const Rcpp::NumericMatrix& grm,
        PgSimParams& params) {

    // n represents number of samples
    int n = grm.nrow();
    double ndouble = static_cast<double>(n);
    
    double numerator = heritability * static_cast<double>(n * (n - 1));
    double denom = (static_cast<double>(1) - heritability);

    double grm_trace = 0;
    rpgsim::Status status = matrix_trace(grm, grm_trace);
    if (status != rpgsim::Status::Success)
        return status;

    double grm_sum = 0;
    status = matrix_sum(grm, grm_sum);
    if (status != rpgsim::Status::Success)
        return status;

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
Rcpp::NumericVector pg_sim_qtl(const bconn_t bconn,
        const Rcpp::NumericMatrix& grm,
        const float qtl_freq,
        const float qtl_effect_size,
        const float heritability,
        const char* id) {

    if (!is_open(bconn))
        Rcpp::stop(bcfio::status_msg(bcfio::Status::ErrBcfNotOpen));

    if (heritability <= 0 || heritability >= 1)
        Rcpp::stop(rpgsim::status_msg(rpgsim::Status::ErrHeritabilityOutOfRange));

    uint32_t nsamps = num_samples(bconn);
    bcfio::bid_t bid = bcfio::replicate(bconn.get());
    if (!bid)
        Rcpp::stop("Internal data structure error");

    PgSimParams pars { 1, 1 }; 
    rpgsim::Status status = calc_pgsim_pars(heritability, grm, pars);
    if (status == rpgsim::Status::Success)
        Rcpp::stop(rpgsim::status_msg(status));

    Rcpp::NumericVector phenotypes(nsamps);

    return phenotypes;
}
