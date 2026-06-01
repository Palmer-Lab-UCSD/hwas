
#include <hwas_types.h>


/////////////////////////////////////////////////////////////////////
// bcf connection
/////////////////////////////////////////////////////////////////////


// Open a file connection to bcf, vcf, or vcf.gz file
//
// @param filename A character vector of the name, and the path if necessary, of
// the bcf file that a connection will be opened.
// @param mode A character vector with the mode, read (r), TODO
// @return pointer to open bcf file connection
// 
// [X] Documentation: man/hts_conn.Rd
//
// [[Rcpp::export]]
Rcpp::Nullable<Rcpp::XPtr<bcfio::Bcf>> bopen(const char* filename, const char* mode) {
    htslib::htsFile* fid = htslib::hts_open(filename, mode);
    if (!fid)
        return R_NilValue;

    return Rcpp::XPtr<bcfio::Bcf>(new bcfio::Bcf(filename, fid), true);
}

// [X] Documentation: man/hts_conn.Rd
//
// [[Rcpp::export]]
int bclose(Rcpp::XPtr<bcfio::Bcf> bid) {
    if (!bid || !bid.get()) return -1;

    bid->close();
    return 0;
}

// [X] Documentation: man/hts_conn.Rd
//
// [[Rcpp::export]]
bool is_open(Rcpp::XPtr<bcfio::Bcf> bid) {
    return bid->is_open();
}


// [X] Documentation: man/hts_conn.Rd
//
// [[Rcpp::export]]
bool is_bcf(const char* filename) {
    return bcfio::is_bcf(filename);
}


/////////////////////////////////////////////////////////////////////
// bcf header
/////////////////////////////////////////////////////////////////////

// [X] Documentation: man/hts_header.Rd
//
// [[Rcpp::export]]
double k_fmt(Rcpp::XPtr<bcfio::Bcf> bid, const char* format_id) {
    return static_cast<double>(bid->hdr_.k_fmt(format_id));
}

// [X] Documentation: man/hts_header.Rd
//
// [[Rcpp::export]]
uint32_t num_samples(Rcpp::XPtr<bcfio::Bcf> bid) {
    return bid->hdr_.n_samples();
}

// [X] Documentation: man/hts_header.Rd
//
// [[Rcpp::export]]
int64_t num_positions(Rcpp::XPtr<bcfio::Bcf> bid) {
    return bcfio::num_records(bid.get());
}

// [X] Documentation: man/hts_header.Rd
//
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

/////////////////////////////////////////////////////////////////////
// bcf config
/////////////////////////////////////////////////////////////////////

// [] Documentation: man/hts_config.Rd
//
// [[Rcpp::export]]
int subset_samples(Rcpp::XPtr<bcfio::Bcf> bid, const char* filename) {
    return htslib::bcf_hdr_set_samples(bid->hdr_.hts_hdr_, filename, 1);
}

// [] Documentation: man/hts_config.Rd
//
// [[Rcpp::export]]
int set_threads(Rcpp::XPtr<bcfio::Bcf> bid, int n) {
    return htslib::hts_set_threads(bid->fid_, n);
}

/////////////////////////////////////////////////////////////////////
// bcf query records
/////////////////////////////////////////////////////////////////////

// [] Documentation: man/hts_records.Rd
//
// [[Rcpp::export]]
Rcpp::RObject next_record(Rcpp::XPtr<bcfio::Bcf> bid, const char* id) {
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
