
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
//
// [[Rcpp::export]]
bcf_conn_t bopen(const char* filename, const char* mode) {
    if (strcmp(mode, "r") != 0) {
        fprintf(stderr, "ERROR: Mode %s is not yet supported.\n");
        return R_NilValue;
    }

    std::unique_ptr<bcfio::Bcf> bid = bcfio::bopen(filename, mode);
    if (bid == nullptr)
       return bcf_conn_t(nullptr, true); 
    return bcf_conn_t(bid.release(), true);
}

// [[Rcpp::export]]
int bclose(bcf_conn_t bconn) {
    if (!bconn || !bconn.get()) return -1;

    bconn->close();
    return 0;
}

// [[Rcpp::export]]
bool is_open(bcf_conn_t bconn) {
    return bcfio::is_open(bconn.get());
}


// [[Rcpp::export]]
bool is_bcf(const char* filename) {
    return bcfio::is_bcf(filename);
}


/////////////////////////////////////////////////////////////////////
// bcf header
/////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
uint16_t k_fmt(bcf_conn_t bconn, const char* id) {
    uint16_t k = 0;
    if (bcfio::k_fmt(bconn.get(), id, &k) != 0)
        return R_NilValue;
    return k;
}

// [[Rcpp::export]]
uint32_t num_samples(bcf_conn_t bconn) {
    uint32_t nsamps = 0;
    if (bcfio::num_samples(bconn.get(), &nsamps) != 0)
        return R_NilValue;
    return nsamps;
}

// [[Rcpp::export]]
int64_t num_positions(bcf_conn_t bconn) {
    int64_t npos = 0;
    if (bcfio::num_pos(bconn.get(), &npos) != 0)
        return R_NilValue;
    return npos;
}

// [[Rcpp::export]]
Rcpp::RObject sample_names(bcf_conn_t bconn) {
    if (!bconn || !bconn.get())
        return R_NilValue;

    uint32_t nsamples = 0;
    if (bcfio::num_samples(bconn.get(), &nsamples) != 0)
        return R_NilValue;

    Rcpp::CharacterVector samp_names(nsamples);

    for (uint32_t i = 0; i < nsamples; i++)
        samp_names[i] = std::string(bconn->hdr->samples[i]);

    return samp_names;
}

/////////////////////////////////////////////////////////////////////
// bcf config
/////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int subset_samples(bcf_conn_t bconn, 
        Rcpp::CharacterVector samples) {

    if (samples.size() == 1 && samples[0] == R_NaString)
        return bcfio::subset_samples(bconn.get(), nullptr);

    int i = 0;
    for (; i < s.size() && s[i] == R_NaString; i++)
        ;
    
    if (i == s.size())
        return -1;

    std::string samp_str = Rcpp::as<std::string>(samples[i]);
    i++;
    for (; i < samples.size(); i++) {
        if (samples[i] == R_NaString)
            continue;

        samp_str += ',' + Rcpp::as<std::string>(samples[i]);
    }

    return bcfio::subset_samples(bconn.get(), samp_str.c_str());
}


int subset_samples_from_file(bcf_conn_t bconn,
        const char* samples_filename) {
    if (!bconn || !samples_filename)
        return -1;
    return bcfio::subset_samples_from_file(bconn.get(), 
            samples_filename);
}

// [[Rcpp::export]]
int set_threads(bcf_conn_t bconn, int n) {
    return htslib::hts_set_threads(bconn->fid, n);
}

/////////////////////////////////////////////////////////////////////
// bcf query records
/////////////////////////////////////////////////////////////////////



// [] Documentation: man/hts_records.Rd
// TODO: update to new interface
// [[Rcpp::export]]
Rcpp::RObject next_record(bcf_conn_t bconn, const char* id) {

    bcfio::Bcf* bid = bconn->checked_get();

    htslib::bcf_fmt_t* fmt_cfg = htslib::bcf_get_fmt(bid->hdr, id);

    // htslib macro (htslib/vcf.h line 1259
    // TODO int_id, i need to transform id string to integer
    // representation for htslib
    //
    // uint32_t fmt_type = bcf_hdr_id2type(bid->hdr_.hts_hdr_, 
    //         BCF_HL_FMT, 
    //         int_id);
    // int status = -1;


    Rcpp::NumericMatrix data;

    switch (fmt_cfg->type) {
    case BCF_HT_REAL:
        data = get_float_matrix(bid, id, fmt_cfg);
        break;
    case BCF_HT_INT:
        data = get_int_matrix();
        break;
    default:
        data = R_NilValue;
    }
    return data;
}
