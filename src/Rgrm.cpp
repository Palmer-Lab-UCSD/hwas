
#include <Rgrm.h>


// [[Rcpp::export]]
Rcpp::NumericMatrix calc_gsim(const bconn_t bconn,
        const char* id) {
    // reset signal for subsequent usage.
    if (!is_open(bconn) || !id)
        Rcpp::stop(status_msg(Status::ErrInvalidInput));

    bcfio::bid_t bid = bcfio::replicate(bconn.get());

    bcfio::BcfHdrAttr hattr {};
    Status status = bcfio::decode_hts_idinfo(bid->hdr,
            id,
            BCF_HL_FMT,
            &hattr);
    if (status != Status::Success)
        Rcpp::stop(status_msg(status));

    switch (hattr.type) {
    case BCF_HT_REAL:
        return calc_gsim_<float>(bid.get(), id);
    case BCF_HT_INT:
        return calc_gsim_<int>(bid.get(), id);
    }

    Rcpp::stop("Unsupported Bcf Type"); 
}



// [[Rcpp::export]]
void gsim_to_grm(Rcpp::NumericMatrix& gsim) {
    double mat_sum = 0;
    Status status = matutils::matrix_sum(gsim, mat_sum);
    if (status != Status::Success)
        Rcpp::stop(status_msg(status));

    double mat_tr = 0;
    status = matutils::matrix_trace(gsim, mat_tr);
    if (status != Status::Success)
        Rcpp::stop(status_msg(status));

    int n = gsim.nrow();
    double norm_const = (mat_tr - mat_sum / n) / (n - 1);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            gsim(i, j) = gsim(i, j) / norm_const;
    }

    gsim.attr("class") = "grm";
}
