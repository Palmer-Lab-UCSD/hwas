//
// By: Robert Vogel
// Affiliation: Palmer Lab at UCSD
// Date: 2025-01-09
//
//

#include <hwas_types.h>



////////////////////////////////////////////////////////////////////
// Postions: load positions file
////////////////////////////////////////////////////////////////////
 
constexpr uint8_t POSIT_BUFF_SIZE = 100;

bcfio::Positions::Positions(uint64_t n) 
    : size(n), pos(size > 0 ? new Position[size] : nullptr) {};

bcfio::Positions::~Positions() { 
    if (pos) delete[] pos;
};

// @title parse the positions file
// @description The postion file is a tab delimited text file
//      with two columns.  The first column is the chromosome/
//      contig string, the second column is the locus of interest,
//      i.e. genomic position.
// @param posits 
// @return  bcfio::Positions*   success
//          nullptr             error
std::unique_ptr<bcfio::Positions> bcfio::parse_posits(textio::FileIO* fio) {

    // lc = line count
    uint64_t lc = 0;
    if (textio::count_lines(fio, &lc) < 0)
        return nullptr;

    std::unique_ptr<bcfio::Positions> posits = std::make_unique<bcfio::Positions>(lc);
    textio::Buffer<char> cbuf { bcfio::MAX_POSIT_BUFF };
    textio::Buffer<uint64_t> ibuf { bcfio::MAX_COL_NUM };
    
    int status = -1;
    while (textio::get_line(fid, cbuf) != 0) {
        status = parse_line(cbuf, ibuf, delimiter = '\t');

        if (status < 0)
            return nullptr;

        if (ibuf->buf_idx_ != 2)
            return nullptr 

        posits[i] = bcfio::Position(std::string(cbuf->buf_ + ibuf->buf_[0]),
                                    std::atoll(cbuf->buf_ + ibuf->buf_[1]));
    }

    return std::move(posits);
}

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
int bcfio::next_record(bcfio::Bcf* bid, bcfio::BcfFloatRecord *ptr) {

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
        return -1;

    return ptr->load_data_(&(bid->hdr_), id);
}


//      -1 error: file 
//      -2 error: unpack error
int bcfio::unpack(bcfio::Bcf* bid, bcfio::BcfFloatRecord* ptr, const char* id) {
    if (!bid || !prt || !id)
        return -1;

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

    int status = bcfio::next_record(bid.get(), &rec);
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
uint64_t num_samples(Rcpp::XPtr<bcfio::Bcf>) {
    if (bid->hdr_->pos_)
        return bid->hdr_->pos_->size;
    
    uint64_t n = 0; 
    while ();
}

// [[Rcpp::export]]
int subset_posits(Rcpp::XPtr<bcfio::Bcf> bid, const char* filename) {
    FileIO fio = textio::fiopen(filename, "r");
    bid->hdr_->pos_ = parse_posits(fio);
    return 0;
}
