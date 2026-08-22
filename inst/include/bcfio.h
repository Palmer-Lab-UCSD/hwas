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
#include <csignal>
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


namespace bcfio {

enum struct Status: int {
    Success         = 0,
    ErrBcfNotOpen   = -3,
    ErrInvalidArg   = -4,
    ErrHtslib       = -5,
    ErrNotBcf       = -6,
    ErrNotEOF       = -7
};


struct SignalInterupt {
    SignalInterupt() {
        std::signal(SIGINT, handler)

    void handler(int signal) {
        signal_recieved = 1;
    }

    std::sig_atomic_t signal_recieved = 0;
};

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
    uint64_t number : 20, vl_type : 4, type : 4, coltype : 4; 
};


class HFileReadConn {
public:
    HFileReadConn()                                 = delete;
    HFileReadConn(const HFileReadConn&)             = delete;
    HFileReadConn& operator=(const HFileReadConn&)  = delete;
    HFileReadConn(HFileReadConn&&)                  = delete;
    HFileReadConn& operator=(HFileReadConn&&)       = delete;

    ~HFileReadConn()
private:
    HFileReadConn(htslib::hFILE* hfid): fid_(hfid) {};

    htslib::hFILE* fid_;
    bool closed_ = false;
    int close_err_ = 0;

    friend bool is_bcf(const char* filename);
};

template <typename T>
using hfile_conn_t = std::unique_ptr<HFileReadConn>

hfile_conn_t hread(const char* filename);


// @title Interface with htslib bcf
// @description ReadBCF manages the lifetime of an open htslib file
//      and organizes the bcf file header and any one record for easy
//      and memory safe parsing.
// @param bcfname: the path and filename to the bcf file to be read.
// @param sample_fname: the path and filename of the text file listing
//  samples id's of records to be retreived.  If this is not included
//  all sample records are retrieved.
class Bcf
{
public:
    Bcf()                       = delete;
    Bcf(const Bcf&)             = delete;
    Bcf& operator=(const Bcf&)  = delete;
    Bcf(Bcf&&)                  = delete;
    Bcf& operator=(Bcf&&)       = delete;

    ~Bcf() { close(); };

    void close();

    // @return C-string copy of Bcf filename
    const char* filename() const;

private:
    Bcf(htslib::htsFile* fid) : hfid(fid) {};

    htslib::htsFile* hfid;

    bool closed = false;
    int close_err = 0;

    friend bid_t bopen(const char* filename, const char* mode);
    friend void close(Bcf* bid);
    friend bool is_open(const Bcf* bid);

    friend Status num_pos(const Bcf* bid, int64_t* n);
    friend const char* chrom(const Bcf* bid, 
            const BcfRecord<T>* brec);

    friend Status next_record<T>(Bcf* bid,
            BcfRecord<T>* brec,
            const char* id);
};

typedef std::unique_ptr<Bcf> bid_t;


class BcfHeader {
public:
    static std::unique_ptr<BcfHeader> init(const Bcf* bid);

    BcfHeader()                             = delete;
    BcfHeader(const BcfHeader&)             = delete;
    BcfHeader& operator=(const BcfHeader&)  = delete;
    BcfHeader(BcfHeader&&)                  = delete;
    BcfHeader& operator=(BcfHeader&&)       = delete;

    ~BcfHeader();

    // @title Retrieve the number of records for the specified format
    // @param id: pointer to character string of a valid format id, e.g.
    //  "GT" would recover genotypes, "HD" haplotype dose, etc.
    // @param k: the allocated memory for which the result is stored
    // @return see enum struct Status definition above
    Status k_fmt(const char* id, uint16_t* k) const;

    // @title Subset samples to samples enumerated in file
    // @param filename: name of file with sample names enumerated
    //  one per line
    // @return ErrInvalidArg, ErrHtslib, WarnSampleListMismatch, Success
    Status subset_samples_from_file(const char* filename);

    // @title Subset samples to samples in comma delimited string
    // @param samples: samples names to included in a comma delimited 
    //  string, if the string starts with ^ then an exclusion list.  If
    //  sample(s) in the list are not found in the bcf file, these names
    //  are ignored.
    // @return ErrHtslib, WarnSampleListMismatch, Success
    Status subset_samples(const char* samples);

    // @title Retrieve the number of samples 
    // @description The number of samples is not the number of samples 
    //  for which the bcf file holds records, but instead the number of
    //  samples that hts header is configured to retrieve data.
    // @param n: memory to write results
    // @return ErrInvalidArg, Success
    Status num_samples(uint32_t* n);

private:
    BcfHeader(htslib::bcf_hdr_t* hdr) hdr_(hdr) {};

    htslib::bcf_hdr_t* hdr_;

    friend Status num_pos(const Bcf* bid, int64_t* n);
};

typedef std::unique_ptr<BcfHeader> bhdr_t;


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
    // no public class constructor, use factory function `init` 
    static brec_t<T> init();

    BcfRecord()=delete
    BcfRecord(const BcfRecord&)=delete;
    BcfRecord& operator=(const BcfRecord&)=delete;
    BcfRecord(BcfRecord&&)=delete;
    BcfRecord& operator=(BcfRecord&&)=delete;
    ~BcfRecord() 

private:
    htslib::bcf1_t* rec;

    // These attributes below store htslib access points to rec data
    // The number of values in memory
    int data_cap = 0;

    // an array of length ndata with values copied from the loaded 
    // record, data are stored in row major form, with a row being 
    // the data for a single sample.  
    //
    // As the data are copies from the htslib::bcf1_t struct on
    // the heap, freeing memory should have no effect on the bcf1_t
    // data.  See htslib/vcf.c line 6228 for the copy operation in
    // htslib::bcf_get_format_values function source code.  I did a
    // quick validation test of this to confirm my understanding.
    T* data = nullptr;
    
    uint32_t nrow = 0;  // n samples
    uint16_t ncol = 0;  // k_fmt

    BcfRecord(htslib::bcf1_t* bcf_rec): rec(bcf_rec) {};

    friend const char* chrom(const Bcf* bid, 
            const BcfRecord<T>* brec);
    friend int pos(const BcfRecord<T>* brec, int64_t* p);
    friend int next_record<T>(Bcf* bid,
            BcfRecord<T>* brec,
            const char* id);
l
};

template <typename T>
using brec_t = std::unique_ptr<BcfRecord<T>>;

template <typename T>
brec_t<T> BcfRecord<T>::init() {
    htslib::bcf1_t* rec = htslib::bcf_init();
    if (rec == nullptr)
        return nullptr;

    BcfRecord<T>* brec = new(std::nothrow) BcfRecord<T>(rec);
    if (brec == nullptr) {
        htslib::bcf_destroy(rec);
        return nullptr
    }

    return brec_t<T>(brec);
}

template <typename T>
BcfRecord<T>::~BcfRecord() {
    if (rec != nullptr) htslib::bcf_destroy(rec);
    if (data != nullptr) free(data);
    rec = nullptr;
    data = nullptr;
    data_cap = 0; 
    ncol = 0;
    nrow = 0;
}



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
int decode_hts_idinfo(const htslib::bcf_hdr_t* hdr,
        const char* id,
        const int bcf_dt_type,
        BcfHdrAttr* ptr);

//////////////////////////////////////////////////////////////////
// API
//////////////////////////////////////////////////////////////////


// Check whether file is vcf or bcf
bool is_bcf(const char* filename);
// 
// // @title open a connection to an hts genotype file
// // @param filename: name of the vcf, vcf.gz, or bcf file to open
// // @param mode: mode to open the file, e.g. "r"
// // @return nullptr on error otherwise unique_ptr to opened file
// //  managed by the Bcf class
bid_t bopen(const char* filename, const char* mode);
bool is_open(const Bcf* bid);

// @title Total number of genomic positions in bcf
// @param bid: pointer to open file connection
// @param n: pointer to integer that stores result
// @return Success,
//  ErrInvalidArg,
//  ErrBcfNotOpen,
//  ErrHtslib, or
//  WarnSampleListMismatch
Status num_pos(Bcf* bid, int64_t* n);

// int32_t exclude_samples(Bcf* bid, const char* sample_filename);

// // @title Query the next records at the next position
// // @param bid: the pointer to open htslib file
// // @param rec: the pointer to the location in memory that data 
// //  are loaded
// // @param the format id for data to be loaded into record.
// // @return: 0 for success otherwise non-zero
template <typename T>
int next_record(Bcf* bid,
        BcfRecord<T>* brec, 
        const char* id) {

    if (!brec) {
        fprintf(stderr,
                "Error: brec is nullptr");
        return -1;
    }

    // reset, so that failure modes report no data loaded.
    brec->nrow = 0;
    brec->ncol = 0;
    
    if (!is_open(bid) || !is_brec_open(brec) || !id) {
        fprintf(stderr,
                "Error: invalid input, one or more arguments are nullptr");
        return -2;
    }

    // reset bcf record to prepare reading for new data
    htslib::bcf_clear(brec->rec);

    int status = htslib::bcf_read(bid->fid, 
            bid->hdr,
            brec->rec);
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
    //
    // For simplicity just unpack all values
    if (htslib::bcf_unpack(brec->rec, BCF_UN_ALL) < 0)
        return -3;

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
        return n;

    uint32_t n_rec_vals = static_cast<uint32_t>(n);

    uint32_t nsamps = 0;
    if (num_samples(bid, &nsamps) < 0)
        return -5;

    htslib::bcf_fmt_t* fmt_cfg = htslib::bcf_get_fmt(bid->hdr, 
            brec->rec, id);

    // detect error in reading data
    // recall that fmt_cfg->n is the number of values per sample
    // TODO: NEED TO CHECK for OVERFLOW
    if (n_rec_vals != fmt_cfg->n * nsamps) {
        return -6;
    }

    uint16_t k = 0;
    if (k_fmt(bid, id, &k) < 0) {
        brec->ncol = brec->nrow = 0;
        return k;
    }

    if (k != fmt_cfg->n)
        return -7;

    brec->nrow = nsamps;
    brec->ncol = k;

    return 0;
}


// @title: Retrieve the genomic position of a record
// @return < 0 upon error and 0 upon success
template <typename T>
int pos(const BcfRecord<T>* brec, int64_t* p) {
    if (!brec) {
        fprintf(stderr,
                "Error: BcfRecord is nullptr\n");
        return -1;
    }

    if (!brec->rec) {
        fprintf(stderr,
                "Error: BcfRecord invalid internal state."
                "Reload record.\n");
        return -1;
    }

    if (p == nullptr) {
        fprintf(stderr,
                "Error: pointer to int64_t encoded position is"
                "nullptr, an invalid value.");
        return -2;
    }

    *p = brec->rec->pos + 1;
    return 0;
}


// @return nullptr upon error, cstring upon success
template <typename T>
const char* chrom(const Bcf* bid, const BcfRecord<T>* brec) {
    if (!is_open(bid))
        return nullptr;

    if (!is_valid_brec(brec))
        return nullptr;

    const char* chr = htslib::bcf_hdr_id2name(bid->hdr, brec->rec->rid);
    if (chr == nullptr)
        return nullptr;

    size_t n = strlen(chr);
    char* chr_return_val = new(std::nothrow) char[n + 1];
    if (chr_return_val == nullptr)
        return nullptr
    
    return std::strncpy(chr_return_val, chr, n + 1);
}

}

#endif
