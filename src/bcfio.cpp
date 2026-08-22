//
//

#include <bcfio.h>

bcfio::HFileReadConn::~HFileReadConn() {
    if (fid != nullptr) 
        close_err_ = htslib::hclose(fid_);

    if (fid_ != nullptr && close_err_ == 0) {
        fid_ = nullptr;
        closed_ = true;
    }
}

bcfio::hfile_conn_t bcfio::hread(const char* filename) {
    if (filename == nullptr)
        return nullptr;

    htslib::hFILE* fh = htslib::hopen(filename, "r");
    if (fh == nullptr)
        return nullptr;

    bcfio::HFileReadConn* hconn = 
        new(std::nothrow) bcfio::HFileReadConn(fh);
    if (hconn == nullptr)
        return nullptr;

    return bcfio::hfile_conn_t(hconn);
}


void bcfio::Bcf::close() noexcept {
    if (hfid_ != nullptr)
        close_err = htslib::hts_close(hfid_);

    if (hfid_ != nullptr && close_err == 0) {
        closed = true;
        hfid_ = nullptr;
    }

}


const char* bcfio::Bcf::filename(const Bcf* bid) const {
    if (fid == nullptr || closed)
        return nullptr;

    // remember that n does not include null character
    size_t n = std::strlen(fid->fn) + 1;
    char* fname = new(std::nothrow) char[n];

    if (fname == nullptr)
        return nullptr;

    return std::strncpy(fname, fid->fn, n);
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
    return bid != nullptr && bid->fid != nullptr && !bid->closed; 
}


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

bcfio::bhdr_t bcfio::BcfHeader::init() {
    if (!is_open(bid))
        return nullptr;

    htslib::bcf_hdr_t* hdr = htslib::bcf_hdr_read(bid->hfid_);
    if (hdr == nullptr)
        return nullptr;

    BcfHeader* bhdr = new(std::nothrow) BcfHeader(hdr);
    if (bhdr == nullptr)
        return nullptr;

    return bhdr_t(bhdr);
}

bcfio::BcfHeader::~BcfHeader() {
    if (hdr_ != nullptr) {
        htslib::bcf_hdr_destroy(hdr_);
        hdr_ = nullptr;
    }
}


bcfio::Status bcfio::BcfHeader::k_fmt(const char *id, 
        uint16_t* k) const {
    if (id == nullptr || k == nullptr)
        return Status::ErrInvalidArg;

    bcfio::BcfHdrAttr fmt {};
    int status = bcfio::decode_hts_idinfo(hdr_, id, BCF_HL_FMT, &fmt);
    if (status < 0)
        return Status::ErrHtslib;

    *k = static_cast<uint16_t>(fmt.number);
    return Status::Success;
}


// htslib accepts a file name with samples to include / exclude or
// a list of comma delimited sample names
bcfio::Status bcfio::BcfHeader::subset_samples_from_file(bcfio::Bcf* bid, 
        const char* filename) {
    if (filename == nullptr)
        return Status::ErrInvalidArg;

    // Recall that 1 indicates that samples are enumerated in file
    int status = htslib::bcf_hdr_set_samples(hdr_,
            filename, 
            1);
    if (status < 0)
        return Status::ErrHtslib;

    if (status > 0)
        return Status::WarnSampleListMismatch;

    return Status::Success;
}


bcfio::Status bcfio::BcfHeader::subset_samples(const char* samples) {
    if (samples == nullptr)
        samples = NULL;

    status = htslib::bcf_hdr_set_samples(hdr_, samples, 0);
    if (status < 0)
        return Status::ErrHtslib;

    if (status > 0)
        return Status::WarnSampleListMismatch;

    return Status::Success;
}


int bcfio::BcfHeader::num_samples(uint32_t* n) {
    if (n == nullptr)
        return Status::ErrInvalidArg;
    *n = static_cast<uint32_t>(hdr_->n[BCF_DT_SAMPLE]);
    return Status::Success;
}


bcfio::Status bcfio::num_pos(const bcfio::Bcf* bid, int64_t* n) {
    if (!bcfio::is_open(bid) || n == nullptr)
        return Status::ErrInvalidArg;

    // open a new file handle, then I can iterate without affecting
    // the current position of bid
    bcfio::bid_t fid = bcfio::bopen(bid->filename(), "r");
    if (!fid)
        return Status::ErrBcfNotOpen;

    bcfio::bhdr_t bhdr = bcfio::BcfHeader::init(fid);
    // bcf_hdr_set_samples
    bcfio::Status status = bhdr->subset_samples(nullptr);
    if (status != Status::Success)
        return status;

    // dummy record
    htslib::bcf1_t* brec = htslib::bcf_init();
    if (!brec)
        return Status::ErrHtslib;

    fprintf(stderr,
            "Note: This may take a time."
            "Expect ~ 30s per 250K positions in ~ 7GB file.\n");
    bcfio::SignalInterupt sigint {};

    status = htslib::bcf_read(fid->hfid_, bhdr->hdr_, brec);
    int64_t npos = 1;
    for (; status == 0; npos++)
        status = htslib::bcf_read(fid->hfid_, bhdr->hdr_, brec);

    htslib::bcf_destroy(brec);
    
    // remember that status == -1 is EOF
    if (status != -1) return Status::ErrNotEOF;

    *n = npos - 1;
    return Status::Success;
}


bool bcfio::is_bcf(const char* filename) {
    bcfio::hfile_conn_t hconn = bcfio::hread(filename);
    if (hconn == nullptr)
        return false;

    htslib::htsFormat fmt {};
    if (htslib::hts_detect_format(hconn->fid_, &fmt) != 0)
        return false;

    if (fmt.format == htslib::bcf || fmt.format == htslib::vcf)
        return true;

    return false;
}


