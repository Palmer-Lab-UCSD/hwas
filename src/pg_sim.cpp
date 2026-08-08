
#include <hwas.h>

struct PolygenicParameters {
    double marker_var = 0;
    double env_var = 0;
};

// Not public facing and consquently assume that bcf_filename
// is a bcf file, heritability is valid
int estimate_pars(bcfio::Bcf* input_bid,
        float heritability, 
        PolygenicParameters* pars) {

    // bcf_conn_t local_bid = bopen(input_bid->fname_.c_str(), "r");
    
    printf("NOT IMPLEMENTED\n");
    return 1;
}

// [[Rcpp::export]]
Rcpp::NumericVector pg_sim(const char* bcf_filename, 
        float heritability) {

    if (heritability < 0 || heritability > 1)
        return Rcpp::NumericVector();

    bcfio::Bcf* bid = bcfio::bopen(bcf_filename, "r");

    if (!bid)
        return Rcpp::NumericVector();

    PolygenicParameters pars {};
    int success = estimate_pars(bid, heritability, &pars);

    uint64_t nsamps = num_samples(bid);
    Rcpp::NumericVector phenotypes(nsamps);

    return phenotypes;
}
