// manage bcf file input / output
//
// *NOTE*  
// I intend to add data write capabilities but need to table for
// now as reading files is the priority.
//
// This library manages the lifetime of htslib data structures
// with wrapper classes and RAII.  Any user of the library should
// only need to refer to the API, see API secition below.
//
//
#ifndef BCFIO_H
#define BCFIO_H

#include <cctype>
#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <cstring>

#include <memory>
#include <new>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

namespace htslib {
extern "C" {
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
}
}

// TODO: I NEED TO ADD A FEATURE FOR LISTIING FORMAT ID's

namespace bcfio {

// TODO: migrate int return values for status codes to the 
// enum struct Status element.  long term slowly integrate
enum struct Status : int {
    WarnSampleSetMismatch   = 4,
    Success                 = 0,
    EndOfFile               = -1,
    ErrHtslib               = -4,
    ErrBcfNotOpen           = -5,
    ErrBcfRecordInvalid     = -6,
    ErrInternal             = -7,
    ErrInvalidInput         = -8,
    ErrParseBcf             = -9,
    ErrInvalidId            = -10,
    ErrBcfOpenFailure       = -11
};

const char* status_msg(Status status);

// @title The meta data on a BCF attribute
// @description BCF, VCF, and VCF.GZ files hold metadata in the
//  header that specify the type and format of data in records.
//  I call each unique piece of data in a record a record attribute, 
//  e.g. an INFO column or FORMAT column of a record are attributes
//  of that record.  HTSLIB encodes attribute information in an
//  unsigned 64 bit integer, and to access any value one needs to 
//  correctly implement bit shifting and masking.  This struct contains 
//  bit-fields representing each value stored in the uint64_t.
// @bitfield number: the number of distinct values required to
//  specify a sample 
//  record at loci i.  For example, a SNP genotype is specified by a
//  single
//  string, e.g. 0/1, while the posterior genotype (0/0, 0/1, 1/1) 
//  probabilities requires three numbers.
// @bitfield vl_type: Specifies whether a variable is fixed length 
//  (BCF_VL_FIXED, in htslib/vcf.h line 68), variable length, etc.
// @bitfield type: the type of variable: binary flag (BCF_HT_FLAG), 
//  integer, real number, string, and 64 bit integers.  Note that HT 
//  is header type.
// @bitfield coltype:
struct BcfHdrAttr { 
    uint64_t number : 20,
    vl_type : 4, 
    type : 4, 
    coltype : 4;
};


class HFileReadConn {
public:
    HFileReadConn()                                 = delete;
    HFileReadConn(const HFileReadConn&)             = delete;
    HFileReadConn& operator=(const HFileReadConn&)  = delete;
    HFileReadConn(HFileReadConn&&)                  = delete;
    HFileReadConn& operator=(HFileReadConn&&)       = delete;

    ~HFileReadConn();

    bool is_bcf() const;
private:
    HFileReadConn(htslib::hFILE* hfid): fid_(hfid) {};

    friend std::unique_ptr<HFileReadConn> hread(const char* filename);
    htslib::hFILE* fid_;
};

typedef std::unique_ptr<HFileReadConn> hfile_conn_t;
hfile_conn_t hread(const char* filename);


// @title: Interface and manage htslib bcf1_t
// @description: The htslib bcf1_t data structure requires manual memory
//  management, knowledge of several bit-packed values, knowledge of
//  several functions for querying data.  This class simplifies 
//  memory management using C++ RAII idiom and provides a simplified,
//  albeit non-comprehensive, interface for loading and querying data
//  stored in the bcf1_t struct.
template <typename T>
class BcfRecord {
public:
    BcfRecord()                                 = delete;
    BcfRecord(const BcfRecord&)                 = delete;
    BcfRecord& operator=(const BcfRecord&)      = delete;
    BcfRecord(BcfRecord&&)                      = delete;
    BcfRecord& operator=(BcfRecord&&)           = delete;

    static std::unique_ptr<BcfRecord<T>> init(); 
    ~BcfRecord();

    htslib::bcf1_t* rec;

    // These attributes below store htslib access points to rec data
    // The number of values in memory
    int data_cap = 0;

    // an array of length ndata with values from the loaded record,
    // data are stored in row major form, with a row being the data
    // for a single sample.
    T* data = nullptr;
    
    uint32_t nrow = 0;  // n samples
    uint16_t ncol = 0;  // k_fmt

private:
    BcfRecord(htslib::bcf1_t* hrec): rec(hrec) {};
};

template <typename T>
using brec_t = std::unique_ptr<BcfRecord<T>>;

template <typename T>
brec_t<T> BcfRecord<T>::init() {
    htslib::bcf1_t* hrec = htslib::bcf_init();
    if (hrec == nullptr)
        return nullptr;

    BcfRecord<T>* brec = new(std::nothrow) BcfRecord<T>(hrec);
    if (brec == nullptr) {
        htslib::bcf_destroy(hrec);
        return nullptr;
    }

    return brec_t<T>(brec);
}

template <typename T>
BcfRecord<T>::~BcfRecord() {
    if (rec != nullptr) 
        htslib::bcf_destroy(rec);
    if (data != nullptr) 
        free(data);
    rec = nullptr;
    data = nullptr;
    data_cap = 0; 
    ncol = 0;
    nrow = 0;
}

template <typename T>
bool is_valid_brec(const BcfRecord<T>* brec) {
    bool truth_val = brec != nullptr && brec->rec != nullptr;

    if (truth_val && brec->data_cap == 0) {
        truth_val = truth_val && brec->data == nullptr;
        truth_val = truth_val && brec->ncol == 0;
        truth_val = truth_val && brec->nrow == 0;
    } else if (truth_val && brec->data_cap > 0)
        truth_val = truth_val && brec->data != nullptr;

    return truth_val;
}



//     float qual() const { return rec_->qual; };


// @title Interface with htslib bcf
// @description ReadBCF manages the lifetime of an open htslib file
//      and organizes the bcf file header and any one record for easy
//      and memory safe parsing.
// @param bcfname: the path and filename to the bcf file to be read.
// @param sample_fname: the path and filename of the text file listing
//  samples id's of records to be retreived.  If this is not included
//  all sample records are retrieved.
struct Bcf
{
    Bcf()                               = delete;
    Bcf(const Bcf&)                     = delete;
    Bcf& operator=(const Bcf&)          = delete;
    Bcf(Bcf&&)                          = delete;
    Bcf& operator=(Bcf&&)               = delete;

    ~Bcf() { close(); };

    void close() noexcept;

    htslib::htsFile* fid;
    htslib::bcf_hdr_t* hdr;
private:
    Bcf(htslib::htsFile* hfid, htslib::bcf_hdr_t* hhdr): 
        fid(hfid), hdr(hhdr) {};
    friend std::unique_ptr<Bcf> bread(const char* filename);
};

typedef std::unique_ptr<Bcf> bid_t;


// @param bcf_dt_type: values are defined in htslib/vcf.h, they
//  specify the type of structure data one wishes to query.  For
//  example, structured meta data, or header lines, can be the
//  following key and encoding:
//
//  VCF key string          key encoding
//  ------------------------------------
//  INFO                    BCF_HL_INFO
//  FORMAT                  BCF_HL_FMT
//  FILTER                  BCF_HL_FLT
//
// @param id: query name for a specific bcf meta data record
// @param ptr: where to store the retrieved meta data
// @param status code of operation
Status decode_hts_idinfo(const htslib::bcf_hdr_t* hdr,
        const char* id,
        const int bcf_dt_type,
        BcfHdrAttr* ptr);

//////////////////////////////////////////////////////////////////
// API
//////////////////////////////////////////////////////////////////


// Check whether file is vcf or bcf
bool is_bcf(const char* filename);

// @title open a connection to an hts genotype file
// @param filename: name of the vcf, vcf.gz, or bcf file to open
// @return nullptr on error otherwise unique_ptr to opened file
//  managed by the Bcf class
bid_t bread(const char* filename);
bool is_open(const Bcf* bid);

// @title Filename from Bcf class
// @param bid: pointer to open file connection
// @return C-string of Bcf filename
const char* get_filename(const Bcf* bid);

// @title Retrieve the number of records for the specified format
// @param bid: pointer to open file connection
// @param id: pointer to character string of a valid format id, e.g.
//  "GT" would recover genotypes, "HD" haplotype dose, etc.
// @param k: the allocated memory for which the result is stored
// @return error code < 0, and 0 upon success
Status k_fmt(const Bcf* bid, const char* id, uint16_t* k);

// @title Retrieve the number of samples 
// @description The number of samples is not the number of samples for
//  which the bcf file holds records, but instead the number of samples
//  that bid is configured to retrieve data.
// @param bid: pointer to open file connection
// @param n: memory to write results
// @return error code < 0 and 0 upon success
Status num_samples(const Bcf* bid, uint32_t* n);


// @title Total number of genomic positions in bcf
// @param bid: pointer to open file connection
// @param n: pointer to integer that stores result
// @return Status
Status num_pos(Bcf* bid, int64_t* n);


// @title Subset samples to samples enumerated in file
// @param bid: pointer to Bcf class
// @param sample_filename: filename with sample names enumerated one
//      per line
// @return Status
Status subset_samples_from_file(Bcf* bid, const char* sample_filename);


// @title Subset samples to samples in comma delimited string
// @param bid: pointer to Bcf class
// @param samples: samples names to included in a comma delimited 
//  string, if the string starts with ^ then an exclusion list.  If
//  sample(s) in the list are not found in the bcf file, these names
//  are ignored.
// @return Status
Status subset_samples(Bcf* bid, const char* samples);


// int32_t exclude_samples(Bcf* bid, const char* sample_filename);

// @title Query the next records at the next position
// @param bid: the pointer to open htslib file
// @param rec: the pointer to the location in memory that data 
//  are loaded
// @param the format id for data to be loaded into record.
// @return: status codes
template <typename T>
Status next_record(Bcf* bid,
        BcfRecord<T>* brec, 
        const char* id) {

    brec->ncol = 0;
    brec->nrow = 0;

    if (!is_open(bid))
        return Status::ErrBcfNotOpen;
        
    if (!is_valid_brec(brec))
        return Status::ErrBcfRecordInvalid;

    // htslib/vcf.h line 413 for bcf_read return values
    int status = htslib::bcf_read(bid->fid, 
            bid->hdr,
            brec->rec);
    if (status == -1)
        return Status::EndOfFile;
    
    if (status < -1)
        return Status::ErrHtslib;

    // Unpacking options defined in htslib/vcf.h line 429
    // BCF_UN_STR:      unpack up to ALT, inclusive
    // BCF_UN_FLT:      unpack up to FILTER 
    // BCF_UN_INFO:     unpack up to INFO
    // BCF_UN_FMT:      unpaack FORMAT for each sample
    //
    // BCF_UN_SHR ==> (BCF_UN_STR | BCF_UN_FLT | BCF_UN_INFO)
    // BCF_UN_ALL ==> (BCF_UN_SHR | BCF_UN_FMT)
    //
    // For simplicity just unpack all values
    if (htslib::bcf_unpack(brec->rec, BCF_UN_ALL) < 0)
        return Status::ErrHtslib;

    // The code below makes format values accssible to the users as
    // C array attribute of the BcfRecord type.
    //
    // remember that `if constexpr ` are evaluated at compile time
    // statements that evaluate to false are not include in code for
    // runtime
    //
    // htslib note: the return value n from bcf_get_format_values is
    // then number of values "written", that is the number of values
    // that the record at this loci contains.  This is distinct from
    // the capacity (brec->data_cap) of the buffer where data are
    // written.  As such, upon success n = num vals per sample * num
    // samples.
    //
    int n = -4;
    if constexpr (std::is_same_v<T, int32_t>)
        n = htslib::bcf_get_format_values(bid->hdr, 
                brec->rec, 
                id, 
                (void**)(&brec->data),
                &brec->data_cap, 
                BCF_HT_INT);

    if constexpr (std::is_same_v<T, float>)
        n = htslib::bcf_get_format_values(bid->hdr, 
                brec->rec, 
                id, 
                (void**)(&brec->data),
                &brec->data_cap, 
                BCF_HT_REAL);


    // TODO: to add string support, I need to first come up with
    //  test examples.
    // if constexpr (std::is_same_v<T, char>)
    //     n = htslib::bcf_get_format_values(bid->hdr, 
    //             brec->rec, 
    //             id, 
    //             (void**)(&brec->data),
    //             &brec->data_cap, 
    //             BCF_HT_STR);

    if (n < 0)
        return Status::ErrHtslib;

    uint32_t n_rec_vals = static_cast<uint32_t>(n);

    uint32_t nsamps = 0;
    if (num_samples(bid, &nsamps) != bcfio::Status::Success)
        return Status::ErrInternal;

    htslib::bcf_fmt_t* fmt_cfg = htslib::bcf_get_fmt(bid->hdr, 
            brec->rec, id);
    if (fmt_cfg == nullptr)
        return Status::ErrHtslib;

    // recall that fmt_cfg->n is the number of values per sample
    // TODO: NEED TO CHECK for OVERFLOW
    if (n_rec_vals != fmt_cfg->n * nsamps)
        return Status::ErrInternal;

    uint16_t k = 0;
    if (k_fmt(bid, id, &k) != bcfio::Status::Success)
        return Status::ErrInternal;

    if (k != fmt_cfg->n)
        return Status::ErrInternal;

    brec->nrow = nsamps;
    brec->ncol = k;

    return Status::Success;
}


// @title: Retrieve the genomic position of a record
// @return < 0 upon error and 0 upon success
template <typename T>
int pos(const BcfRecord<T>* brec, int64_t* p) {
    if (!bcfio::is_valid_brec(brec) || p == nullptr)
        return -1;

    *p = brec->rec->pos + 1;
    return 0;
}


// @return nullptr upon error, cstring upon success
template <typename T>
const char* chrom(const Bcf* bid, const BcfRecord<T>* brec) {
    if (!is_open(bid))
        return nullptr;

    if (!is_valid_brec<T>(brec))
        return nullptr;

    return htslib::bcf_hdr_id2name(bid->hdr, brec->rec->rid);
}

}

#endif
