//
// By: Robert Vogel
// Affiliation: Palmer Lab at UCSD
// Date: 2025-01-09
//
// Input argument
//    filename: vcf with haplotpye
//
//
//

#include <hwas_types.h>

constexpr uint8_t POSIT_BUFF_SIZE = 100

///////////////////////////////////////////////////////////////////
// BcfHeader
///////////////////////////////////////////////////////////////////
//
int bcfio::BcfHeader::decode_hts_idinfo_(const char *name, 
        const int bcf_dt_type, 
        bcfio::BcfHdrAttr *ptr) const {

    // BCF_DT_ID is the C macro for the ID dictionary index defined 
    // by htslib see htslib/vcf.h line 86
    int idx = htslib::bcf_hdr_id2int(hts_hdr_, BCF_DT_ID, name);

    if (idx < 0)
        return idx;

    uint64_t val = hts_hdr_->id[BCF_DT_ID][idx].val->info[bcf_dt_type];

    ptr->number = val >> 12 & 0xfffff;
    ptr->vl_type = val >> 8 & 0xf;
    ptr->type = val >> 4 & 0xf;
    ptr->coltype = val & 0xf;

    return 0;
}

int bcfio::BcfHeader::get_format_attr(const char *id, BcfHdrAttr *ptr) const {
    return decode_hts_idinfo_(id, BCF_HL_FMT, ptr);
}

int bcfio::BcfHeader::get_info_attr(const char *id, BcfHdrAttr *ptr) const {
    return decode_hts_idinfo_(id, BCF_HL_INFO, ptr);
}

int bcfio::BcfHeader::get_filter_attr(const char *id, BcfHdrAttr *ptr) const {
    return decode_hts_idinfo_(id, BCF_HL_FLT, ptr);
}

int32_t bcfio::BcfHeader::k_fmt(const char *id) const {
    if (!id)
        return -1;

    BcfHdrAttr fmt {};

    int32_t status { -1 };

    if ((status = get_format_attr(id, &fmt)) < 0)
        return status;

    return static_cast<int32_t>(fmt.number);
}


// const std::unique_ptr<std::string[]> bcfio::BcfHeader::sample_names() const {
// 
//     std::unique_ptr<std::string[]> samp_names = 
//         std::make_unique<std::string[]>(n_samples()); 
// 
//     for (size_t i = 0; i < n_samples(); i++)
//         samp_names[i] = std::string(*(hdr_->samples + i));
// 
//     return samp_names;
// }
// 
// ///////////////////////////////////////////////////////////////////
// // BcfFloatRecord
// ///////////////////////////////////////////////////////////////////
// 
bcfio::BcfFloatRecord::~BcfFloatRecord() {
    if (rec_) htslib::bcf_destroy(rec_);
    if (dst_) free(dst_);
    rec_ = nullptr;
    dst_ = nullptr;
}

const float* bcfio::BcfFloatRecord::array() const {
    return dst_; 
}

std::optional<float> bcfio::BcfFloatRecord::get(const uint64_t row_idx,
        const uint64_t col_idx) const {
    size_t idx = row_idx * col_num_ + col_idx;
    if (idx >= size()) return std::nullopt;

    return *(dst_ + idx);
}

int bcfio::BcfFloatRecord::load_data_(bcfio::BcfHeader *hdr, const char *id) {
    int status { 0 };
    col_num_ = row_num_ = 0;

    status = htslib::bcf_get_format_values(hdr->hts_hdr_, 
            rec_, 
            id, 
            (void**)(&dst_),
            &ndst_, 
            BCF_HT_REAL);

    if (status < 0)
        return status;

    int32_t k { 0 };
    if ((k = hdr->k_fmt(id)) < 0) {
        col_num_ = row_num_ = 0;
        return k;
    }

    col_num_ = static_cast<uint64_t>(k);
    row_num_ = hdr->n_samples();

    return 0;
}

////////////////////////////////////////////////////////////////////
// Postions: load positions file
////////////////////////////////////////////////////////////////////
 
bcfio::Positions::Positions(FILE* fid)
    : fid_(fid) {
    if (fid_) {
        buffer = new char[buffsize+1];
        std::memset(buffer, '\0', sizeof(char) * (buffsize+1));
    }
}

bcfio::Positions::~Positions() { 
    if (fid) fclose(fid); 
    if (buffer) delete[] buffer;
    if (pos) delete[] pos;
};


//@return   0       success
//          -1      error: file stream ended on value != EOF
//          -2      error: failed to return to file handle origin
int bcfio::count_lines(bcfio::Positions* posits) {
    uint64_t counter = 0;
    int prev_char = '\0';
    int c = '\0';
    FILE* fid = posits->fid_;

    // don't count empty lines
    while ((c = getchar(fid)) != EOF) {
        if (c == '\n' && prev_char != '\n') counter++;
        prev_char = c;
    }

    // EOF is serving as a end of line character
    // that is why we need to increment
    if (c == EOF && prev_char != '\n')
        counter++;

    // return file handle to the beginning of file for parsing
    if (std::fseek(fid, 0, SEEK_SET) != 0) return -2;
    if (c != EOF) return -1;

    posits->size = counter;
    return 0;
}


// @title parse the positions file
// @description The postion file is a tab delimited text file
//      with two columns.  The first column is the chromosome/
//      contig string, the second column is the locus of interest,
//      i.e. genomic position.
// @param posits 
// @return  0       success
//          -1      error: counting lines
//          -2      error: column mismatch
int bcfio::parse(bcfio::Positions* posits) {

    if (bcfio::count_lines(posits) < 0)
        return -1;

    uint64_t num_pos = posits->size;
    bcfio::Position* posit_array = new Position[posits->size];

    CharBuffer cbuf { POSIT_BUFF_SIZE };

    std::string chromosome {};
    bool has_read_chrom_col = false;
    uint64_t lc = 0;// line counter
    int c = '\0';
    while ((c = getchar(fid)) != EOF && lc < num_pos) {

        // this means that we have more than two columns
        if (c == '\t' && has_read_chrom_col) {
            Rprintf("[Warning]\tIs there more than two columns at"
                    "%s\n?", cbuf.buf_);
            continue;
        }

        if (c == '\t') {
            cbuf.append('\0');
            chromosome = std::string(cbuf.buf_);
            cbuf.reset();
            has_read_chrom_col = true;
        }

        if (c == '\n' && !has_read_chrom_col) {
            Rprintf("[Warning]\tMissing field with buffer %s.\n",
                    cbuf.buf_);
            continue;
        }

        if (c == '\n') {
            cbuf.append('\0');
            posit_array[] = Position{ chromosome, 
                                    std::atoll(cbuf.buf_)};
            cbuf.reset();
            has_read_chrom_col = false;
        }

        if (cbuf.append(c) < 0) 
            Rcpp::stop("[ERROR]\tbuffer overflow\n");
    }

    return 0;
}

// [[Rcpp::export]]
Rcpp::XPtr<bcfio::Positions> popen(const char* filename, 
                                const char* mode) {
    FILE* fid = fopen(filename, mode);
    if (!fid)
        return R_NilValue;

    bcfio::Positions* posits = new bcfio::Positions(fid);
    if (parse(posits) < 0) {
        delete posits;
        return R_NilValue;
    }

    return Rcpp::XPtr<bcfio::Position>(pos);
}
 
// ///////////////////////////////////////////////////////////////////
// // Bcf
// ///////////////////////////////////////////////////////////////////
// ///
// 
//
bcfio::Bcf::Bcf()
    : fname_(""), fid_(nullptr), hdr_() {};
 
bcfio::Bcf::Bcf(const char *filename, htslib::htsFile* fid)
    : fname_(filename),
    fid_(fid),
    hdr_(fid_) {};

bcfio::Bcf::~Bcf() {
    if (fid_) htslib::hts_close(fid_);
}

bool bcfio::Bcf::isopen() const {
    return fid_ != nullptr;
}
// 
// // TODO: subset samples by those in sample_fname file
// int bcfio::ReadBcf::set_samples(const char *sample_fname) {
// 
//     // Subset samples with those found in the file sample_fname 
//     if (!sample_fname || *sample_fname == '\0') {
//         fprintf(stdout, "No file with sample names detected, retreiving"
//                 " records for all samples.\n");
//         return -1;
//     }
// 
//     return hdr_.subset_samples(sample_fname);
// };
// 
// @title load next record
// @description The algorithm assumes that pos is a subset of possible
//      positions.
//
// return values:
//      0 for success
//      -1 for end of file
//      -2 error: nullptr argument
//      -3 error: record retrieval
//      -4 error: include position is out of order
//      -5 error: bcf positions out of order
//      -6 error: unpack error
int bcfio::next_record(bcfio::Bcf* bid, 
        bcfio::BcfFloatRecord *ptr, 
        const char *id) {

    if (!bid || !ptr || !id)
        return -2;

    int status = -3;
    if (!ptr->dst_)
        int64_t prev_rec_pos_ = -1;
    else
        int64_t prev_rec_pos_ = ptr->rec_->pos;

    int64_t prev_include_pos_ = -1;

    status = htslib::bcf_read(bid->fid_, 
                bid->hdr_.hts_hdr_, 
                ptr->cur_rec());

    bool requested_position { false };

    // if user doesn't specify a position file, then bid->pos 
    // is nullptr and a record at each position is returned in
    // turn.
    if (!bid->pos)
        requested_position = true;
    else if (bid->pos == ptr->pos) {
        requested_position = true;
        include_pos_count_++;
        ++(bid->pos_idx_);
    }

    while (!requested_position && status == 0) {
        if (bid->pos < ptr->pos)
            return bid->pos_idx_ < bid->num_pos_ ? -2 : ;

        if (bid->pos[pos_idx_] == ptr->pos) {
            status = htslib::bcf_read(bid->fid_, 
                                    bid->hdr_.hts_hdr_, 
                                    ptr->cur_rec());
            bid->pos_idx_++;
        }
    }

    if (status != 0)
        return status;

    // Unpacking options defined in htslib/vcf.h line 419
    if (htslib::bcf_unpack(ptr->cur_rec(), BCF_UN_ALL) < 0)
        return -6;

    return ptr->load_data_(&(bid->hdr_), id);
}


// [[Rcpp::export]]
Rcpp::XPtr<bcfio::Bcf> bopen(const char* filename, const char* mode) {
    htslib::htsFile* fid = htslib::hts_open(filename, mode);
    if (!fid)
        Rcpp::stop("File was not found");

    return Rcpp::XPtr<bcfio::Bcf>(new bcfio::Bcf(filename, fid), true);
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
    bcfio::BcfFloatRecord rec {};

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

// [[Rcpp::export]]
int subset_posits(Rcpp::XPtr<bcfio::Bcf> bid, 
                    Rcpp::XPtr<bcfio::Positions> posits) {
    bid.pos_ = posits.get();
    return 0;
}
