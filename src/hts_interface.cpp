
#include <hwas_types.h>


// [[Rcpp::export]]
Rcpp::XPtr<bcfio::Bcf> bopen(const char* filename, const char* mode) {
    htslib::htsFile* fid = htslib::hts_open(filename, mode);
    if (!fid)
        Rcpp::stop("File was not found");

    return Rcpp::XPtr<bcfio::Bcf>(new bcfio::Bcf(filename, fid), true);
}

// [[Rcpp::export]]
void bclose(Rcpp::XPtr<bcfio::Bcf> bid) {
    bid->close();
}

// [[Rcpp::export]]
bool isopen(Rcpp::XPtr<bcfio::Bcf> bid) {
    return bid->isopen();
}

// [[Rcpp::export]]
double k_fmt(Rcpp::XPtr<bcfio::Bcf> bid, const char* id) {
    return static_cast<double>(bid->hdr_.k_fmt(id));
}

// [[Rcpp::export]]
uint32_t num_samples(Rcpp::XPtr<bcfio::Bcf> bid) {
    return bid->hdr_.n_samples();
}

// [[Rcpp::export]]
Rcpp::RObject query_next(Rcpp::XPtr<bcfio::Bcf> bid, const char* id) {
    bcfio::BcfRecord<float> rec {};

    int status = bcfio::next_record(bid.get(), &rec, id);
    if (status != 0)
        return R_NilValue;

    std::optional<float> datum = std::nullopt;
    Rcpp::NumericMatrix expect_hap_counts(rec.nrows(), rec.ncols());

    for (uint64_t i = 0; i < rec.nrows(); i++)
        for (uint64_t j = 0; j < rec.ncols(); j++) {
            if (!(datum = rec.get(i, j)))
                Rcpp::stop("Data retrieval error.");
            expect_hap_counts(i, j) = datum.value();
        }

    return expect_hap_counts;
}

// [[Rcpp::export]]
int subset_samples(Rcpp::XPtr<bcfio::Bcf> bid, const char* filename) {
    return htslib::bcf_hdr_set_samples(bid->hdr_.hts_hdr_, filename, 1);
}

// [[Rcpp::export]]
int set_threads(Rcpp::XPtr<bcfio::Bcf> bid, int n) {
    return htslib::hts_set_threads(bid->fid_, n);
}
