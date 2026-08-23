
#include <hwas.h>


// Open a file connection to bcf, vcf, or vcf.gz file
//
// @param filename A character vector of the name, and the path if necessary, of
// the bcf file that a connection will be opened.
// @param mode A character vector with the mode, read (r), TODO
// @return pointer to open bcf file connection
// 
//
// [[Rcpp::export]]
bcf_conn_t bread(const char* filename) {
    bcfio::bid_t bid = bcfio::bread(filename);

    if (bid == nullptr) {
        bcfio::Status status = bcfio::Status::ErrInvalidInput;
        Rcpp::stop(bcfio::status_msg(status));
    }

    return bcf_conn_t(bid.release(), true);
}

// [[Rcpp::export]]
int bclose(bcf_conn_t bconn) {
    if (!bconn) {
        bcfio::Status status = bcfio::Status::ErrInvalidInput;
        Rcpp::stop(bcfio::status_msg(status));
    }

    bconn->close();
    return 0;
}

// [[Rcpp::export]]
bool is_open(bcf_conn_t bconn) {
    if (!bconn) {
        bcfio::Status status = bcfio::Status::ErrInvalidInput;
        Rcpp::stop(bcfio::status_msg(status));
    }
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
    bcfio::Status status = bcfio::Status::ErrInvalidInput;

    if (!bconn)
        Rcpp::stop(bcfio::status_msg(status));

    uint16_t k = 0;
    status = bcfio::k_fmt(bconn.get(), id, &k);
    if (status != bcfio::Status::Success)
        Rcpp::stop(bcfio::status_msg(status));

    return k;
}

// [[Rcpp::export]]
uint32_t num_samples(bcf_conn_t bconn) {
    bcfio::Status status = bcfio::Status::ErrInvalidInput;
    if (!bconn)
        Rcpp::stop(bcfio::status_msg(status));

    uint32_t nsamps = 0;
    status = bcfio::num_samples(bconn.get(), &nsamps);
    if (status != bcfio::Status::Success)
        Rcpp::stop(bcfio::status_msg(status));

    return nsamps;
}

// [[Rcpp::export]]
int64_t num_positions(bcf_conn_t bconn) {
    bcfio::Status status = bcfio::Status::ErrInvalidInput;
    if (!bconn)
        Rcpp::stop(bcfio::status_msg(status));

    int64_t npos = 0;
    status = bcfio::num_pos(bconn.get(), &npos);
    if (status != bcfio::Status::Success)
        Rcpp::stop(bcfio::status_msg(status));

    return npos;
}

// [[Rcpp::export]]
Rcpp::RObject sample_names(bcf_conn_t bconn) {
    bcfio::Status status = bcfio::Status::ErrInvalidInput;
    if (!bconn)
        Rcpp::stop(bcfio::status_msg(status));

    uint32_t nsamples = 0;
    status = bcfio::num_samples(bconn.get(), &nsamples);
    if (status != bcfio::Status::Success)
        Rcpp::stop(bcfio::status_msg(status));

    Rcpp::CharacterVector samp_names(nsamples);

    for (uint32_t i = 0; i < nsamples; i++)
        samp_names[i] = bconn->hdr->samples[i];

    return samp_names;
}

/////////////////////////////////////////////////////////////////////
// bcf config
/////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int subset_samples(bcf_conn_t bconn, 
        Rcpp::CharacterVector samples) {
    bcfio::Status status = bcfio::Status::ErrInvalidInput;
    if (!bconn)
        Rcpp::stop(bcfio::status_msg(status));

    if (samples.size() == 1 && samples[0] == R_NaString) {
        status = bcfio::subset_samples(bconn.get(), nullptr);

        if (status != bcfio::Status::Success)
            Rcpp::stop(bcfio::status_msg(status));
    }

    // I need to figure out if all the elements of the CharacterVecor
    // are NA values
    int i = 0;
    for (; i < samples.size() && samples[i] == R_NaString; i++)
        ;
    // I expect to enumerate through the entire CharacterVector if
    // each value is R_NaString.  If this happens then i == s.size() 
    if (i == samples.size())
        return -1;

    // get the number of characters required for comma delimited list
    size_t ntotal = 0;
    
    for (i = 0; i < samples.size(); i++) {
        if (samples[i] == R_NaString)
            continue;

        // remember that strlen does not count the null terminator, for
        // our list, all but the last sample string null terminator
        // will be replaced by a comma, the last element will retain 
        // null character meaning that we need to add 1 to strlen return
        // value.
        ntotal += std::strlen(Rcpp::String(samples[i]).get_cstring()) + 1;
    }

    std::unique_ptr<char[]> sample_list = std::make_unique<char[]>(ntotal);
    char* sample_list_elem = sample_list.get();
    char* end_point = sample_list_elem + ntotal - 1;
    const char* sample_str = nullptr;
    size_t idx = 0;
    size_t str_n = 0;
    for (i = 0; i < samples.size(); i++ ) {
        if (samples[i] == R_NaString)
            continue;

        sample_str = Rcpp::String(samples[i]).get_cstring();
        str_n = std::strlen(sample_str);

        if (sample_list_elem + str_n > end_point)
            Rcpp::stop("yikes!");

        std::strncpy(sample_list_elem, sample_str, str_n);
        sample_list_elem[str_n] = ',';
        sample_list_elem += str_n + 1;
    }
    sample_list[ntotal-1] = '\0';

    status = bcfio::subset_samples(bconn.get(), sample_list.get());
    if (status != bcfio::Status::Success)
        Rcpp::stop(bcfio::status_msg(status));

    return 0;
}


int subset_samples_from_file(bcf_conn_t bconn,
        const char* samples_filename) {
    bcfio::Status status = bcfio::Status::ErrInvalidInput;
    if (!bconn)
        Rcpp::stop(bcfio::status_msg(status));

    status = bcfio::subset_samples_from_file(bconn.get(), 
            samples_filename);

    if (status != bcfio::Status::Success)
        Rcpp::stop(bcfio::status_msg(status));

    return 0;
}

// [[Rcpp::export]]
int set_threads(bcf_conn_t bconn, int n) {
    bcfio::Status status = bcfio::Status::ErrInvalidInput;
    if (!bconn)
        Rcpp::stop(bcfio::status_msg(status));

    if (htslib::hts_set_threads(bconn->fid, n) != 0)
        Rcpp::stop(bcfio::status_msg(bcfio::Status::ErrHtslib));

    return 0;
}


// [] Documentation: man/hts_records.Rd
// TODO: update to new interface
// [[Rcpp::export]]
Rcpp::Nullable<Rcpp::NumericMatrix> next_record(bcf_conn_t bconn, 
        const char* id) {
    bcfio::Status status = bcfio::Status::ErrInvalidInput;
    if (!bconn)
        Rcpp::stop(bcfio::status_msg(status));

    if (!bcfio::is_open(bconn.get()))
        Rcpp::stop(bcfio::status_msg(bcfio::Status::ErrBcfNotOpen));

    bcfio::BcfHdrAttr hattr {};
    status = bcfio::decode_hts_idinfo(bconn->hdr,
            id,
            BCF_HL_FMT,
            &hattr);
    if (status != bcfio::Status::Success)
        Rcpp::stop(bcfio::status_msg(status));

    if (hattr.type == BCF_HT_REAL)
        return get_matrix<float>(bconn.get(), id);

    if (hattr.type == BCF_HT_INT)
        return get_matrix<int>(bconn.get(), id);

    Rcpp::stop("Unsupported Bcf Type"); 
}

