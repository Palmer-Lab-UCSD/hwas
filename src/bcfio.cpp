//
//

#include <bcfio.h>


bcfio::HFileReadConn::~HFileReadConn() {
    if (fid_) 
        static_cast<void>(htslib::hclose(fid_)); 
}

bool bcfio::HFileReadConn::is_bcf() const {
    htslib::htsFormat fmt {};
    if (htslib::hts_detect_format(fid_, &fmt) != 0)
        return false;

    if (fmt.format == htslib::bcf || fmt.format == htslib::vcf)
        return true;

    return false;
}

bcfio::hfile_conn_t bcfio::hread(const char* filename) {
    htslib::hFILE* fh = htslib::hopen(filename, "r");
    if (fh == nullptr)
        return nullptr;

    bcfio::HFileReadConn* hfile = 
        new(std::nothrow) bcfio::HFileReadConn(fh);
    if (hfile == nullptr) {
        static_cast<void>(htslib::hclose(fh));
        return nullptr;
    }

    return bcfio::hfile_conn_t(hfile);
}

bcfio::Status bcfio::decode_hts_idinfo(const htslib::bcf_hdr_t* hdr,
        const char* id, 
        const int bcf_dt_type, 
        bcfio::BcfHdrAttr* ptr) {

    // BCF_DT_ID is the C macro for the ID dictionary index defined 
    // by htslib see htslib/vcf.h line 86
    int idx = htslib::bcf_hdr_id2int(hdr, BCF_DT_ID, id);
    if (idx == -1)
        return bcfio::Status::ErrInvalidId;

    uint64_t val = hdr->id[BCF_DT_ID][idx].val->info[bcf_dt_type];

    ptr->number = val >> 12 & 0xfffff;
    ptr->vl_type = val >> 8 & 0xf;
    ptr->type = val >> 4 & 0xf;

    // col type is the BCF_HL_* value (line 1252 in htslib/vcf.h)
    ptr->coltype = val & 0xf;

    return bcfio::Status::Success;
}

// // const std::unique_ptr<std::string[]> bcfio::BcfHeader::sample_names() const {
// // 
// //     std::unique_ptr<std::string[]> samp_names = 
// //         std::make_unique<std::string[]>(n_samples()); 
// // 
// //     for (size_t i = 0; i < n_samples(); i++)
// //         samp_names[i] = std::string(*(hdr_->samples + i));
// // 
// //     return samp_names;
// // }
// // 
// ///////////////////////////////////////////////////////////////////
// // template <typename T> BcfRecord
// ///////////////////////////////////////////////////////////////////
// 



// template <typename T>
// std::optional<T> bcfio::BcfRecord<T>::get(const uint64_t row_idx,
//         const uint64_t col_idx) const {
//     size_t idx = row_idx * col_num_ + col_idx;
//     if (idx >= size()) return std::nullopt;
// 
//     return *(dst_ + idx);
// }

 
/////////////////////////////////////////////////////////////////////
//// Bcf
/////////////////////////////////////////////////////////////////////
// 
//
void bcfio::Bcf::close() noexcept {
    if (fid != nullptr) {
        static_cast<void>(htslib::hts_close(fid));
        fid = nullptr;
    }

    if (hdr != nullptr) {
        htslib::bcf_hdr_destroy(hdr);
        hdr = nullptr;
    }
}


/////////////////////////////////////////////////////////////////////
// API
/////////////////////////////////////////////////////////////////////
///

const char* bcfio::get_filename(const Bcf* bid) {
    if (!bid || !bid->fid)
        return nullptr;
    return bid->fid->fn;
}


bcfio::bid_t bcfio::bread(const char* filename) {
    if (!bcfio::is_bcf(filename))
        return nullptr;

    htslib::htsFile* fh = htslib::hts_open(filename, "r");
    if (!fh)
        return nullptr;

    htslib::bcf_hdr_t* hdr = htslib::bcf_hdr_read(fh);
    if (hdr == nullptr) {
        static_cast<void>(htslib::hts_close(fh));
        return nullptr;
    }

    bcfio::Bcf* bid = new(std::nothrow) bcfio::Bcf(fh, hdr);
    if (bid == nullptr) {
        static_cast<void>(htslib::hts_close(fh));
        htslib::bcf_hdr_destroy(hdr);
        return nullptr;
    }

    return bcfio::bid_t(bid);
}


bool bcfio::is_open(const bcfio::Bcf* bid) {
    return bid != nullptr && bid->fid != nullptr && bid->hdr != nullptr;
}


bcfio::Status bcfio::k_fmt(const bcfio::Bcf* bid, 
        const char *id, 
        uint16_t* k) {

    if (!bcfio::is_open(bid))
        return bcfio::Status::ErrBcfNotOpen;
    if (!id || k == nullptr)
        return bcfio::Status::ErrInvalidInput;

    BcfHdrAttr fmt {};
    bcfio::Status status = bcfio::decode_hts_idinfo(bid->hdr, 
            id, 
            BCF_HL_FMT, 
            &fmt);
    if (status != bcfio::Status::Success)
        return status;

    *k = static_cast<uint16_t>(fmt.number);
    return bcfio::Status::Success;
}

// htslib accepts a file name with samples to include / exclude or
// a list of comma delimited sample names
bcfio::Status bcfio::subset_samples_from_file(bcfio::Bcf* bid, 
        const char* samples_filename){
    if (!bcfio::is_open(bid))
        return bcfio::Status::ErrBcfNotOpen;

    if (samples_filename == nullptr)
        return bcfio::Status::ErrInvalidInput;

    // Recall that 1 indicates that samples are enumerated in file
    int status = htslib::bcf_hdr_set_samples(bid->hdr,
            samples_filename, 
            1);
    if (status < 0)
        return bcfio::Status::ErrHtslib;

    if (status > 0)
        return bcfio::Status::WarnSampleSetMismatch;

    return bcfio::Status::Success;
}


bcfio::Status bcfio::subset_samples(bcfio::Bcf* bid, 
        const char* samples) {

    if (!bcfio::is_open(bid))
        return bcfio::Status::ErrBcfNotOpen;

    if (!samples)
        samples = NULL;

    // Recall that 1 indicates that samples are enumerated 
    int status = htslib::bcf_hdr_set_samples(bid->hdr,
            samples, 
            0);

    if (status < 0)
        return bcfio::Status::ErrHtslib;

    if (status > 0)
        return bcfio::Status::WarnSampleSetMismatch;

    return bcfio::Status::Success;
}


bcfio::Status bcfio::num_samples(const bcfio::Bcf* bid, uint32_t* n) {
    if (!bcfio::is_open(bid))
        return bcfio::Status::ErrBcfNotOpen;
       
    if ( n == nullptr)
        return bcfio::Status::ErrInvalidInput;

    *n = static_cast<uint32_t>(bid->hdr->n[BCF_DT_SAMPLE]);
    return bcfio::Status::Success;
}


bcfio::Status bcfio::num_pos(bcfio::Bcf* bid, int64_t* n) {
    if (!bcfio::is_open(bid))
        return bcfio::Status::ErrBcfNotOpen;

    const char* filename = bcfio::get_filename(bid);
    // open a new file handle, then I can iterate without affecting
    // the current position of bid
    bcfio::bid_t fid = bcfio::bread(filename);
    if (fid == nullptr)
        return bcfio::Status::ErrInternal;

    // bcf_hdr_set_samples
    bcfio::Status status = bcfio::subset_samples(fid.get(), nullptr);
    if (status != bcfio::Status::Success)
        return bcfio::Status::ErrInternal;

    // dummy record
    htslib::bcf1_t* brec = htslib::bcf_init();
    if (!brec)
        return bcfio::Status::ErrHtslib;

    int hts_status = htslib::bcf_read(fid->fid, fid->hdr, brec);
    int64_t npos = 1;
    for (; hts_status == 0; npos++)
        hts_status = htslib::bcf_read(fid->fid, fid->hdr, brec);

    htslib::bcf_destroy(brec);
    
    // remember that -1 here is htslib signal for EOF
    if (hts_status != -1) 
        return bcfio::Status::ErrParseBcf;

    *n = npos - 1;
    return bcfio::Status::Success;
}


bool bcfio::is_bcf(const char* filename) {
    bcfio::hfile_conn_t fh = bcfio::hread(filename);
    if (fh == nullptr)
        return false;

    return fh->is_bcf();
}

