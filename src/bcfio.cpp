//
//

#include <bcfio.h>


int bcfio::decode_hts_idinfo(const htslib::bcf_hdr_t* hdr,
        const char* id, 
        const int bcf_dt_type, 
        bcfio::BcfHdrAttr* ptr) {

    // BCF_DT_ID is the C macro for the ID dictionary index defined 
    // by htslib see htslib/vcf.h line 86
    int idx = htslib::bcf_hdr_id2int(hdr, BCF_DT_ID, id);
    if (idx < 0)
        return idx;

    uint64_t val = hdr->id[BCF_DT_ID][idx].val->info[bcf_dt_type];

    ptr->number = val >> 12 & 0xfffff;
    ptr->vl_type = val >> 8 & 0xf;
    ptr->type = val >> 4 & 0xf;

    // col type is the BCF_HL_* value (line 1252 in htslib/vcf.h)
    ptr->coltype = val & 0xf;

    return 0;
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
bcfio::Bcf::Bcf(): fid(nullptr), hdr(nullptr) {};
 
bcfio::Bcf::Bcf(htslib::htsFile* hts_fid)
    : fid(hts_fid),
    hdr(hts_fid ? htslib::bcf_hdr_read(hts_fid) : nullptr) {};

void bcfio::Bcf::close() noexcept {
    if (fid) {
        htslib::hts_close(fid);
        fid = nullptr;
    }

    if (hdr) {
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


// TODO double check return statemnent
bcfio::bid_t bcfio::bopen(const char* filename, const char* mode) {
    if (!bcfio::is_bcf(filename))
        return nullptr;

    htslib::htsFile* fh = htslib::hts_open(filename, mode);
    if (!fh)
        return nullptr;

    bcfio::Bcf* bid = new Bcf(fh);
    if (!bcfio::is_open(bid)) {
        delete bid;
        return nullptr;
    }

    return bcfio::bid_t(bid);
}


bool bcfio::is_open(const bcfio::Bcf* bid) {
    return bid != nullptr && bid->fid != nullptr && bid->hdr != nullptr;
}


int bcfio::k_fmt(const bcfio::Bcf* bid, const char *id, uint16_t* k) {
    if (!bcfio::is_open(bid) || !id || k == nullptr)
        return -1;

    BcfHdrAttr fmt {};
    int status = bcfio::decode_hts_idinfo(bid->hdr, id, BCF_HL_FMT, &fmt);
    if (status < 0)
        return status;

    *k = static_cast<uint16_t>(fmt.number);
    return 0;
}

// title: load next record
// @return -1 end of file, < -1 error, 0 success
// 
// TODO I don't remember why this is necessary, should this go
// in header?
// template struct bcfio::BcfRecord<float>;
// template int bcfio::next_record<float> (bcfio::Bcf* bid,
//         bcfio::BcfRecord<float>* ptr,
//         const char* id);

// htslib accepts a file name with samples to include / exclude or
// a list of comma delimited sample names
int bcfio::subset_samples_from_file(bcfio::Bcf* bid, 
        const char* samples_filename){
    if (!bid || !samples_filename || !bcfio::is_open(bid))
        return -1;

    // Recall that 1 indicates that samples are enumerated in file
    return htslib::bcf_hdr_set_samples(bid->hdr,
            samples_filename, 
            1);
}


int bcfio::subset_samples(bcfio::Bcf* bid, const char* samples) {
    if (!bid || !bcfio::is_open(bid))
        return -1;

    if (!samples)
        samples = NULL;

    // Recall that 1 indicates that samples are enumerated 
    return htslib::bcf_hdr_set_samples(bid->hdr,
            samples, 
            0);
}


int bcfio::num_samples(const bcfio::Bcf* bid, uint32_t* n) {
    if (!bid || n == nullptr)
        return -1;
    *n = static_cast<uint32_t>(bid->hdr->n[BCF_DT_SAMPLE]);
    return 0;
}


int bcfio::num_pos(bcfio::Bcf* bid, int64_t* n) {
    if (!bid)
        return -1;

    const char* filename = bcfio::get_filename(bid);
    // open a new file handle, then I can iterate without affecting
    // the current position of bid
    bcfio::bid_t fid = bcfio::bopen(filename, "r");
    if (!fid)
        return -2;

    // bcf_hdr_set_samples
    int status = bcfio::subset_samples(fid.get(), nullptr);
    if (status != 0)
        return -3;

    // dummy record
    htslib::bcf1_t* brec = htslib::bcf_init();
    if (!brec)
        return -3;

    status = htslib::bcf_read(fid->fid, fid->hdr, brec);
    int64_t npos = 1;
    for (; status == 0; npos++)
        status = htslib::bcf_read(fid->fid, fid->hdr, brec);

    htslib::bcf_destroy(brec);
    
    if (status != -1) return -3;

    *n = npos - 1;
    return 0;
}


bool bcfio::is_bcf(const char* filename) {
    std::unique_ptr<bcfio::HFileReadConn> fh = bcfio::hread(filename);
    if (!fh)
        return false;

    htslib::htsFormat fmt {};
    if (htslib::hts_detect_format(fh->fid, &fmt) != 0)
        return false;

    if (fmt.format == htslib::bcf || fmt.format == htslib::vcf)
        return true;

    return false;
}


std::unique_ptr<bcfio::HFileReadConn> bcfio::hread(const char* filename) {
    htslib::hFILE* fh = htslib::hopen(filename, "r");
    if (!fh)
        return nullptr;

    return std::make_unique<bcfio::HFileReadConn>(fh);
}
