//
//

#include <bcfio.h>

const char* bcfio::status_msg(bcfio::Status status) {
    switch (status) {
    case bcfio::Status::WarnEmptyLine:
        return "Warning: Line is empty";
    case bcfio::Status::WarnSampleSetMismatch:
        return "Warning: Sample list contains names not in bcf";
    case bcfio::Status::Success:
        return "Success";
    case bcfio::Status::EndOfFile:
        return "Reached end of file.";
    case bcfio::Status::ErrNotImplemented:
        return "Not implemented.";
    case bcfio::Status::ErrHtslib:
        return "Likely a problem with htslib interface, please"
            " contact the maintainers.";
    case bcfio::Status::ErrBcfNotOpen:
        return "Bcf file not open for reading";
    case bcfio::Status::ErrBcfRecordInvalid:
        return "Likely invalid Bcf Record.";
    case bcfio::Status::ErrInternal:
        return "Internal error, please contact maintainers.";
    case bcfio::Status::ErrInvalidInput:
        return "Invalid input value";
    case bcfio::Status::ErrParseBcf:
        return "Error parsing Bcf file, please check whether"
            " the file is correctly formatted.  If formatted"
            " correctly please contact maintainers.";
    case bcfio::Status::ErrInvalidId:
        return "Invalid id for the bcf query";
    case bcfio::Status::ErrBcfOpenFailure:
        return "Failed trying to open file, please check that"
            " the specified file is a valid vcf, vcf.gz, or bcf"
            " formatted file.";
    case bcfio::Status::ErrDuplicatePositions:
        return "Positions file has dupliate positions.";
    case bcfio::Status::ErrParsePositionsFileInvalidCoord:
        return "Invalid contig:pos detected in positions file.";
    case bcfio::Status::ErrParsePositionsFileCoordStrTooLong:
        return "contig:pos string too long";
    case bcfio::Status::ErrCouldNotReadFile:
        return "Could not open file for reading.";
    case bcfio::Status::ErrCouldNotInsertCoordInPosSet:
        return "Could not add coordinate to position set, may be duplicate.";
    case bcfio::Status::ErrParseUnrecoverable:
        return "File egregiously violages expected contents, exit.";
    }

    return "Unexpected status, please contact maintainers.";
}

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


bool bcfio::GenomicCoord::operator<(const bcfio::GenomicCoord& other) const {
    if (ctg < other.ctg)
        return true;

    if (ctg == other.ctg && pos < other.pos)
        return true;

    return false;
}

bool bcfio::GenomicCoord::operator==(const bcfio::GenomicCoord& other) const {
    return ctg == other.ctg && pos == other.pos;
}


bool bcfio::GenomicCoord::operator!=(const bcfio::GenomicCoord& other) const {
    return !(*this == other);
}

bool bcfio::GenomicCoord::operator>(const bcfio::GenomicCoord& other) const {
    if (*this < other)
        return false;
    if (*this == other)
        return false;
    return true;
}



bcfio::pos_file_t bcfio::PositionsFile::read(const char* filename) {
    if (filename == nullptr)
        return nullptr;

    FILE* fid = std::fopen(filename, "r");
    if (fid == nullptr)
        return nullptr;
    
    PositionsFile* pfid = 
        new(std::nothrow) PositionsFile(fid);
    if (pfid == nullptr) {
        std::fclose(fid);
        return nullptr;
    }

    std::memset(pfid->buf_, '\0', pfid->buf_size_ * sizeof(char));

    return bcfio::pos_file_t(pfid);
}


void bcfio::PositionsFile::close() {
    if (fid_ != nullptr) {
        std::fclose(fid_);
        fid_ = nullptr;
    }
}

bcfio::PositionsFile::~PositionsFile() {
    close();
}


bcfio::Status bcfio::PositionsFile::getline() {
    buf_len_ = 0;
    buf_[0] = '\0';
    int c;
    for (;buf_len_ < buf_size_; buf_len_++) {
        c = std::fgetc(fid_); 
        
        if (c == '\n' || c == '\r' || c == EOF) {
            buf_[buf_len_] = '\0';
            break;
        }

        buf_[buf_len_] = c;
    }

    // note, the remaining characters in the line are still in read buffer
    // we need to flush these values before returning;
    if (buf_len_ == buf_size_) {

        buf_[0] = '\0';
        buf_len_ = 0;

        const int max_itr = 10000;
        int j = 0;
        for (; j < max_itr; j++) {
            c = std::fgetc(fid_);
            if (c == '\n' || c == '\r' || c == EOF)
                break;
        }
        if (j == max_itr)
            return bcfio::Status::ErrParseUnrecoverable;

        return bcfio::Status::ErrParsePositionsFileCoordStrTooLong;
    }

    if (c == EOF && buf_len_ == 0)
        return bcfio::Status::EndOfFile;

    if (buf_len_ == 0)
        return bcfio::Status::WarnEmptyLine;

    return bcfio::Status::Success;
}


bcfio::Status bcfio::PositionsFile::next_record(bcfio::GenomicCoord* gc) {

    gc->ctg = std::string("");
    gc->pos = 0;

    bcfio::Status status = getline();
    if (status != bcfio::Status::Success)
        return status;

    int i = 0;

    // parse contig, remember that contig is terminated by ':'
    // I don't need to check the last non-null element of the as
    // if the character is or is not ':' there cannot be a position
    // string present => Invalid coordinate
    for (i = 0; i < buf_len_ && buf_[i] != ':'; i++)
        ;

    // Cases:
    //  ':' is first character => no contig specified
    //  ':' is last character => no position specified
    if (i == 0 || i == buf_len_)
        return bcfio::Status::ErrParsePositionsFileInvalidCoord;


    // remember that atoll returns zero on failure.  Moreover, I think
    // the coordinate numbers in the vcf/bcf are 1-based, so zero is
    // not a valid contig position.
    int64_t pos = std::atoll(buf_ + i+1);
    if (pos == 0)
        return bcfio::Status::ErrParsePositionsFileInvalidCoord;

    buf_[i] = '\0';

    gc->ctg = std::string(buf_);
    gc->pos = pos;

    return bcfio::Status::Success;
}

const char* bcfio::PositionsFile::buf() const {
    return buf_;
}

bool bcfio::is_open(const PositionsFile* pfid) {
    if (pfid == nullptr)
        return false;
    return pfid->fid_ != nullptr;
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

bcfio::Status bcfio::num_samples(const bcfio::Bcf* bid, uint32_t* n) {
    if (!bcfio::is_open(bid))
        return bcfio::Status::ErrBcfNotOpen;
       
    if ( n == nullptr)
        return bcfio::Status::ErrInvalidInput;

    *n = static_cast<uint32_t>(bid->hdr->n[BCF_DT_SAMPLE]);
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
    htslib::bcf1_t* rec = htslib::bcf_init();
    if (!rec)
        return bcfio::Status::ErrHtslib;

    int hts_status = 0;
    int64_t npos = 0;
    if (bid->pos.empty()) {
        while (true) {
            hts_status = htslib::bcf_read(fid->fid, fid->hdr, rec);
            if (hts_status != 0)
                break;
            npos++;
        }
    } else {
        const char* ctg = nullptr;
        while (true) {
            hts_status = htslib::bcf_read(fid->fid, fid->hdr, rec);
            if (hts_status != 0)
                break;

            ctg = htslib::bcf_hdr_id2name(fid->hdr, rec->rid);
            bcfio::GenomicCoord gc {std::string(ctg), rec->pos + 1};

            if (bid->pos.count(gc) == 1)
                npos++;
        }
    }
    
    htslib::bcf_destroy(rec);

    // remember that -1 here is htslib signal for EOF
    if (hts_status != -1) 
        return bcfio::Status::ErrParseBcf;

    *n = npos;
    return bcfio::Status::Success;
}


bcfio::Status bcfio::set_pos_from_file(bcfio::Bcf* bid, const char* filename) {
    if (!bcfio::is_open(bid))
        return bcfio::Status::ErrInvalidInput;

    // empty all contents, if any, in the pos set
    bid->pos.clear();

    bcfio::pos_file_t pfid = bcfio::PositionsFile::read(filename);
    if (pfid == nullptr)
        return bcfio::Status::ErrCouldNotReadFile;

    GenomicCoord gc {};

    int64_t idx = 1;
    bcfio::Status status = bcfio::Status::ErrInternal;
    while ((status = pfid->next_record(&gc)) != bcfio::Status::EndOfFile) {

        if (status == bcfio::Status::ErrParseUnrecoverable)
            return status;

        if (status != bcfio::Status::Success) {
            fprintf(stderr, "Line %ld, record %s are excluded due"
                    " to error:\n%s\n", 
                    idx,
                    pfid->buf(),
                    bcfio::status_msg(status));
            continue;
        }
        
        // Assume that each insert makes a copy of gc
        auto res = bid->pos.insert(gc);
        if (!res.second) {
            if (bid->pos.count(gc) == 0)
                return bcfio::Status::ErrCouldNotInsertCoordInPosSet;
            else {
                fprintf(stderr, "Line %ld, record %s excluded due"
                        " to error:\n%s\n", 
                        idx,
                        pfid->buf(),
                        bcfio::status_msg(bcfio::Status::ErrDuplicatePositions));
            }
        }
    }

    return bcfio::Status::Success;
}


bool bcfio::is_bcf(const char* filename) {
    bcfio::hfile_conn_t fh = bcfio::hread(filename);
    if (fh == nullptr)
        return false;

    return fh->is_bcf();
}

