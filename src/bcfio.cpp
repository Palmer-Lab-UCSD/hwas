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

#include <bcfio.h>


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

void bcfio::BcfHeader::close() noexcept {
    if (!isnull()) {
        htslib::bcf_hdr_destroy(hts_hdr_);
        hts_hdr_ = nullptr;
    }
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
// // template <typename T> BcfRecord
// ///////////////////////////////////////////////////////////////////
// 
template <typename T>
bcfio::BcfRecord<T>::~BcfRecord() {
    if (rec_) htslib::bcf_destroy(rec_);
    if (dst_) free(dst_);
    rec_ = nullptr;
    dst_ = nullptr;
}

template <typename T>
const T* bcfio::BcfRecord<T>::array() const {
    return dst_; 
}

template <typename T>
std::optional<T> bcfio::BcfRecord<T>::get(const uint64_t row_idx,
        const uint64_t col_idx) const {
    size_t idx = row_idx * col_num_ + col_idx;
    if (idx >= size()) return std::nullopt;

    return *(dst_ + idx);
}

template <typename T>
int bcfio::BcfRecord<T>::load_data_(bcfio::BcfHeader *hdr, 
        const char *id) {
    int status = -1;
    col_num_ = row_num_ = 0;

    status = htslib::bcf_get_format_values(hdr->hts_hdr_, 
            rec_, 
            id, 
            (void**)(&dst_),
            &ndst_, 
            bcf_record_type);

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

 
/////////////////////////////////////////////////////////////////////
//// Bcf
/////////////////////////////////////////////////////////////////////
// 
//
bcfio::Bcf::Bcf()
    : fname_(""), fid_(nullptr), hdr_() {};
 
bcfio::Bcf::Bcf(const char *filename, htslib::htsFile* fid)
    : fname_(filename),
    fid_(fid),
    hdr_(fid_) {};

// for the bcf instance to be open, then both the
// header and file connection resources must be in 
// valid states
bool bcfio::Bcf::is_open() {
    if (!fid_) {
        hdr_.close();
        return false;
    }

    // if, for some reason the hdr_ is closed but the file handle
    // isn't, then close the file handle before return false
    if (hdr_.isnull()) {
        htslib::hts_close(fid_);
        fid_ = nullptr;
        return false;
    }

    return true;
}

void bcfio::Bcf::close() noexcept {
    if (is_open()) {
        htslib::hts_close(fid_);
        fid_ = nullptr;
        hdr_.close();
    }
}


// title: load next record
template <typename T>
int bcfio::next_record(bcfio::Bcf* bid,
        bcfio::BcfRecord<T>* ptr, 
        const char* id) {

    if (!bid->is_open())
        return -1;

    int status = htslib::bcf_read(bid->fid_, 
            bid->hdr_.hts_hdr_,
            ptr->cur_rec());
    if (status != 0)
        return status;

    // Unpacking options defined in htslib/vcf.h line 429
    // BCF_UN_STR:      unpack up to ALT, inclusive
    // BCF_UN_FLT:      unpack up to FILTER 
    // BCF_UN_INFO:     unpack up to INFO
    // BCF_UN_FMT:      unpaack FORMAT for each sample
    //
    // BCF_UN_SHR ==> (BCF_UN_STR | BCF_UN_FLT | BCF_UN_INFO)
    // BCF_UN_ALL ==> (BCF_UN_SHR | BCF_UN_FMT)
    if (htslib::bcf_unpack(ptr->cur_rec(), BCF_UN_ALL) < 0)
        return -1;

    return ptr->load_data_(&(bid->hdr_), id);
}


template struct bcfio::BcfRecord<float>;
template int bcfio::next_record<float> (bcfio::Bcf* bid,
        bcfio::BcfRecord<float>* ptr,
        const char* id);

int64_t bcfio::num_records(bcfio::Bcf* bid) {
    if (!bid)
        return -1;
    // open a new file handle, then I can iterate without affecting
    // the current position of bid
    htslib::htsFile* fid = htslib::hts_open(bid->fname_.c_str(), "r");
    if (!fid)
        return -2;

    htslib::bcf_hdr_t* hdr = htslib::bcf_hdr_read(fid);
    if (!hdr) {
        htslib::hts_close(fid);
        return -2;
    }

    // dummy record
    htslib::bcf1_t* brec = htslib::bcf_init();
    if (!brec) {
        htslib::bcf_hdr_destroy(hdr);
        htslib::hts_close(fid);
        return -2;
    }

    uint64_t n = 0;
    int status;
    for (n = 0; status == 0; n++)
        status = htslib::bcf_read(fid, hdr, brec);

    htslib::bcf_destroy(brec);
    htslib::bcf_hdr_destroy(hdr);
    htslib::hts_close(fid);
    
    // report error occured while parsing
    if (status == -2) return status;

    return n-1;
}


bool bcfio::is_bcf(const char* filename) {
    htslib::hFILE* fh = htslib::hopen(filename, "r");
    if (!fh)
        return false;

    htslib::htsFormat fmt {};
    if (htslib::hts_detect_format(fh, &fmt) != 0) {
        htslib::hclose(fh);
        return false;
    }

    if (fmt.format == htslib::bcf) {
        htslib::hclose(fh);
        return true;
    }
        
    htslib::hclose(fh);
    return false;
}
