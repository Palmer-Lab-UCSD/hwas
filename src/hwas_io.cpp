
#include <hwas.h>


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
bcf_conn_t bopen(const char* filename, const char* mode) {
    bcfio::Bcf* bid = bcfio::bopen(filename, mode);
    if (bid == nullptr)
       return bcf_conn_t(nullptr, true); 
    return bcf_conn_t(bid, true);
}

// [X] Documentation: man/hts_conn.Rd
//
// [[Rcpp::export]]
int bclose(bcf_conn_t bid) {
    if (!bid || !bid.get()) return -1;

    bid->close();
    return 0;
}

// [X] Documentation: man/hts_conn.Rd
//
// [[Rcpp::export]]
bool is_open(bcf_conn_t bid) {
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
double k_fmt(bcf_conn_t bid, const char* id) {
    return static_cast<double>(bid->hdr_.k_fmt(id));
}

// [X] Documentation: man/hts_header.Rd
//
// [[Rcpp::export]]
uint32_t num_samples(bcf_conn_t bid) {
    return bcfio::num_samples(bid.get());
}

// [X] Documentation: man/hts_header.Rd
//
// [[Rcpp::export]]
int64_t num_positions(bcf_conn_t bid) {
    return bcfio::num_pos(bid.get());
}

// [X] Documentation: man/hts_header.Rd
//
// [[Rcpp::export]]
Rcpp::RObject sample_names(bcf_conn_t bid) {
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
int subset_samples(bcf_conn_t bid, 
        Rcpp::CharacterVector samples) {

    bool is_na_val = samples[0] == R_NaString;
    if (samples.size() == 1 && is_na_val)
        return bcfio::subset_samples(bid.get(), nullptr);

    if (is_na_val)
        return -1;

    std::string samp_str = Rcpp::as<std::string>(samples[0]);

    // check for missing values
    for (int i = 1; i < samples.size(); i++) {
        if (samples[i] == R_NaString)
            return -1;

        samp_str += ',' + Rcpp::as<std::string>(samples[i]);
    }

    return bcfio::subset_samples(bid.get(), samp_str.c_str());
}


int subset_samples_from_file(bcf_conn_t bid,
        const char* samples_filename) {
    if (!bid)
        return -1;
    return bcfio::subset_samples_from_file(bid.get(), 
            samples_filename);
}

// [] Documentation: man/hts_config.Rd
//
// [[Rcpp::export]]
int set_threads(bcf_conn_t bid, int n) {
    return htslib::hts_set_threads(bid->fid_, n);
}

/////////////////////////////////////////////////////////////////////
// bcf query records
/////////////////////////////////////////////////////////////////////

// [] Documentation: man/hts_records.Rd
//
// [[Rcpp::export]]
Rcpp::RObject next_record(bcf_conn_t bid, const char* id) {
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
