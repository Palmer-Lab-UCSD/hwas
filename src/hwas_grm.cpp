
#include <hwas.h>


// For signal handling see
// https://en.cppreference.com/cpp/utility/program/signal
// I need to review this more closely.  The use of atomic
// indicates that this is using mutex free concurrency.
// volatile std::sig_atomic_t signal_received = 0;

// void signal_handler(int signal) {
//     signal_received = 1;
// }


// [[Rcpp::export]]
Rcpp::NumericMatrix calc_unnormalized_grm(const char* filename, 
        const char* id) {
    // reset signal for subsequent usage.
    if (!filename || !id)
        Rcpp::stop(bcfio::status_msg(bcfio::Status::ErrInvalidInput));

    bcfio::bid_t bid = bcfio::bread(filename);
    if (!bcfio::is_open(bid))
        Rcpp::stop(bcfio::status_msg(bcfio::Status::ErrBcfOpenFailure));

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
        data = calc_unnormalized_grm<float>(bid, id);
        break;
    case BCF_HT_INT:
        data = calc_unnormalized_grm<int>(bid, id);
        break;
    default:
        Rcpp::stop("Unsupported Bcf Type"); 
    }

    return data;
}
