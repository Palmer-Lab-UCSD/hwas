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
bcfio::Bcf::Bcf(): fid_(nullptr), hdr_() {};
 
bcfio::Bcf::Bcf(htslib::htsFile* fid)
    : fid_(fid), hdr_(fid) {};

// bcfio::Bcf::Bcf(const bcfio::Bcf& bid) {
// }

std::string bcfio::Bcf::filename() const {
    if (!fid_)
        return std::string();
    return std::string(fid_->fn);
}


// for the bcf instance to be open, then both the header and file
// resources must be in valid states.
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


std::unique_ptr<bcfio::Bcf> bcfio::bopen(const char* filename, 
        const char* mode) {
    if (!bcfio::is_bcf(filename))
        return nullptr;

    htslib::htsFile* fh = htslib::hts_open(filename, mode);
    if (!fh)
        return nullptr;

    return std::make_unique<bcfio::Bcf>(fh);
}


// title: load next record
template <typename T>
int bcfio::next_record(bcfio::Bcf* bid,
        bcfio::BcfRecord<T>* ptr, 
        const char* id) {

    if (!bid || !bid->is_open())
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

// TODO I don't remember why this is necessary, should this go
// in header?
template struct bcfio::BcfRecord<float>;
template int bcfio::next_record<float> (bcfio::Bcf* bid,
        bcfio::BcfRecord<float>* ptr,
        const char* id);

// htslib accepts a file name with samples to include / exclude or
// a list of comma delimited sample names
int bcfio::subset_samples_from_file(bcfio::Bcf* bid, const char* samples_filename){
    if (!bid || !samples_filename)
        return -1;

    // Recall that 1 indicates that samples are enumerated in file
    return htslib::bcf_hdr_set_samples(bid->hdr_.hts_hdr_,
            samples_filename, 
            1);
}


int bcfio::subset_samples(bcfio::Bcf* bid, const char* samples){
    if (!bid)
        return -1;

    if (!samples)
        samples = NULL;

    // Recall that 1 indicates that samples are enumerated in file
    return htslib::bcf_hdr_set_samples(bid->hdr_.hts_hdr_,
            samples, 
            0);
}


int64_t bcfio::num_samples(bcfio::Bcf* bid) {
    if (!bid)
        return -1;
    return static_cast<int64_t>(bid->hdr_.n_samples());
}

int64_t bcfio::num_pos(bcfio::Bcf* bid) {
    if (!bid)
        return -1;
    // open a new file handle, then I can iterate without affecting
    // the current position of bid
    std::unique_ptr<bcfio::Bcf> fid;
    fid = bcfio::bopen(bid->filename().c_str(), "r");

    if (!fid)
        return -2;
    // bcf_hdr_set_samples
    int status = bcfio::subset_samples(fid.get(), nullptr);
    if (status != 0)
        return -1;

    // dummy record
    htslib::bcf1_t* brec = htslib::bcf_init();
    if (!brec)
        return -2;

    int64_t n = 0;
    for (n = 0; status == 0; n++)
        status = htslib::bcf_read(fid->fid_, fid->hdr_.hts_hdr_, brec);

    htslib::bcf_destroy(brec);
    
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
        static_cast<void>(htslib::hclose(fh));
        return false;
    }

    if (fmt.format == htslib::bcf ||
            fmt.format == htslib::vcf) {
        static_cast<void>(htslib::hclose(fh));
        return true;
    }
        
    static_cast<void>(htslib::hclose(fh));
    return false;
}

