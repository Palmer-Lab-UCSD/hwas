
#include <Rgrm.h>



// [[Rcpp::export]]
Rcpp::NumericMatrix calc_unnormalized_grm(const bconn_t bconn,
        const char* id) {
    // reset signal for subsequent usage.
    if (!is_open(bconn) || !id)
        Rcpp::stop(bcfio::status_msg(bcfio::Status::ErrInvalidInput));

    bcfio::bid_t bid = bcfio::replicate(bconn.get());

    bcfio::BcfHdrAttr hattr {};
    bcfio::Status status = bcfio::decode_hts_idinfo(bid->hdr,
            id,
            BCF_HL_FMT,
            &hattr);
    if (status != bcfio::Status::Success)
        Rcpp::stop(bcfio::status_msg(status));

    Rcpp::NumericMatrix data;
    switch (hattr.type) {
    case BCF_HT_REAL:
        return calc_unnormalized_grm_<float>(bid.get(), id);
    case BCF_HT_INT:
        return calc_unnormalized_grm_<int>(bid.get(), id);
    }

    Rcpp::stop("Unsupported Bcf Type"); 
}
