
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
int64_t num_positions(Rcpp::XPtr<bcfio::Bcf> bid) {
    return bcfio::num_records(bid.get());
}


// [[Rcpp::export]]
Rcpp::RObject sample_names(Rcpp::XPtr<bcfio::Bcf> bid) {
    if (!bid.get())
        return R_NilValue;

    uint32_t nsamples = bid->hdr_.n_samples();
    Rcpp::CharacterVector samp_names(nsamples);

    for (uint32_t i = 0; i < nsamples; i++)
        samp_names[i] = std::string(bid->hdr_.hts_hdr_->samples[i]);

    return samp_names;
}

// [[Rcpp::export]]
int subset_samples(Rcpp::XPtr<bcfio::Bcf> bid, const char* filename) {
    return htslib::bcf_hdr_set_samples(bid->hdr_.hts_hdr_, filename, 1);
}

// [[Rcpp::export]]
int set_threads(Rcpp::XPtr<bcfio::Bcf> bid, int n) {
    return htslib::hts_set_threads(bid->fid_, n);
}

// [[Rcpp::export]]
Rcpp::RObject query_next(Rcpp::XPtr<bcfio::Bcf> bid, const char* id) {
    bcfio::BcfRecord<float> rec {};

    int status = bcfio::next_record(bid.checked_get(), &rec, id);
    if (status != 0)
        return R_NilValue;

    std::optional<float> datum = std::nullopt;
    Rcpp::NumericMatrix data(rec.nrows(), rec.ncols());

    for (uint64_t i = 0; i < rec.nrows(); i++)
        for (uint64_t j = 0; j < rec.ncols(); j++) {
            if (!(datum = rec.get(i, j)))
                Rcpp::stop("Data retrieval error.");
            data(i, j) = datum.value();
        }

    data.attr("chrom") = std::string(rec.chrom(bid->hdr_.hts_hdr_));
    data.attr("pos") = rec.pos();
    //data.attr("qual") = rec.qual();
    //data.attr(

    return data;
}



