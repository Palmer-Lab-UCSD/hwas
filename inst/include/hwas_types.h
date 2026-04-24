// Parse STITCH vcf file
//
// By: Robert Vogel
// Affiliation: Palmer Lab at UCSD
// Date: 2025-01-09
//
//
// Acknowledgment
//
//
#ifndef HEADER_BCFIO_H
#define HEADER_BCFIO_H

#include <cstdlib>
#include <cstdint>
#include <optional>
#include <memory>
#include <string>
#include <Rcpp.h>

namespace htslib {
extern "C" {
#include <htslib/vcf.h>
#include <htslib/hts.h>
}
}

// samples are separated by white space
// const char HAP_CODE[] { "HD" };

namespace bcfio {


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
struct BcfHdrAttr { uint64_t number : 20, vl_type : 4, type : 4, coltype : 4; };


// @title: Interface and manager of htslib bcf_hdr_t
// @description: The bcf header C-struct requires manual allocation 
//  and release of memory.  This class applies RAII, for easy
//  maintenance.
class BcfHeader {
public:
    BcfHeader() 
        : hdr_(nullptr) {};

    BcfHeader(htslib::htsFile *fid)
        : hdr_(fid ? htslib::bcf_hdr_read(fid) : nullptr) {};

    ~BcfHeader() { if (hdr_) htslib::bcf_hdr_destroy(hdr_); };

    bool isnull() const { return hdr_ == nullptr; };


    // @title: Retreive the set of smaple names
    const std::unique_ptr<std::string[]> sample_names() const;

    // @title: "get_*" member functions for attribute retrieval
    // @description:
    // @param id: the id of the formatted data field to retrieve
    // @param ptr: the pointer to memory for which the BcfHdrAttr 
    //  data will be copied into memory.
    // @return 0 for success < 0 for fail
    int get_format_attr(const char *id, BcfHdrAttr *ptr) const;
    int get_info_attr(const char *id, BcfHdrAttr *ptr) const;
    int get_filter_attr(const char *id, BcfHdrAttr *ptr) const;

    int subset_samples(const char *filename);

    // @title: The number of values stored in format id
    // @description: Each bcf format field is able to hold unique
    //  number of values per sample.  This function provides a simple
    //  interface to the bcf file to retrieve that number.
    // @param id: the format field id
    // @return if an error occured that value returned is < 0,
    //  otherwise the number of values of fmt field id recorded per
    //  sample is returned.
    int32_t k_fmt(const char* id) const;

    uint32_t n_samples() const { return hdr_->n[BCF_DT_SAMPLE]; };

    // TODO: what unit test should I do for this?
    const htslib::bcf_hdr_t *hts_hdr() const { return hdr_; };

private:
    htslib::bcf_hdr_t *hdr_;
    // BcfHdrAttr attr_ {};

    // @title: 
    // @description decoder based upon htslib/vcf.h line 100 in the
    //  typedef struct bcf_idinfo_t. 
    // @param name:
    // @param bcf_dt_type
    // @param ptr
    // @return -1 indicates an error has occured and 0 a success
    int decode_hts_idinfo_(const char *name, 
            const int bcf_dt_type, 
            BcfHdrAttr *ptr) const;
};


// @title: Interface and manage htslib bcf1_t
// @description: The htslib bcf1_t data structure requires manual memory
//  management, knowledge of several bit-packed values, knowledge of
//  several functions for querying data.  This class simplifies 
//  memory management using C++ RAII idiom and provides a simplified,
//  albeit non-comprehensive, interface for loading and querying data
//  stored in the bcf1_t struct.
struct BcfFloatRecord {

    BcfFloatRecord(): rec_(htslib::bcf_init()) {};
    ~BcfFloatRecord();

    // provide check-free fast, but unsafe, access to loaded data
    float operator[](const size_t idx) const { return *(dst_ + idx); };

    // provide index checked access to data.
    std::optional<float> get(const uint64_t row_idx, 
            const uint64_t col_idx) const;


    // the total amount of values per record, n_samples * k_founders
    uint64_t size() const { return static_cast<uint64_t>(ndst_); };
    uint64_t ncols() const { return col_num_; };
    uint64_t nrows() const { return row_num_; };

    htslib::bcf1_t *cur_rec() const { return rec_; }; 

    bool is_snp() const { return htslib::bcf_is_snp(rec_); }

    htslib::bcf1_t *rec_;

    // These attributes store htslib access points to record data
    //
    // ndst_: The number of values in memory, the length of dst
    // *dst_: an array of length ndst_ with float values of the
    //  current record
    int ndst_ = 0;
    float *dst_ = nullptr;

    // data that dst_ point to are stored in row major order, with 
    // columns being k_fmt and rows being n_samples.
    uint64_t col_num_ = 0;
    uint64_t row_num_ = 0;


    // @title: Load sample data at the current locus
    // @description: Sample data at the current locus, is not made
    //  available by reading a locus's record and storing in the
    //  bcf1_t type. Instead, we need to supply a pointer variable
    //  and format id to make that id's sample data available. This 
    //  function helps simplify this process.
    // @param hdr: instance of the bcf header to retreive meta data
    // @param tag: the C-string id representing the data we want to 
    //  query.
    // @return 0 upon success and != 0 for failure
    int load_data_(BcfHeader *hdr, const char *tag);
};


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
    Bcf();
    Bcf(const char* filename, htslib::htsFile* fid);

    Bcf(const Bcf&)=delete;
    Bcf& operator=(const Bcf&)=delete;

    Bcf(Bcf&&)=delete;
    Bcf& operator=(Bcf&&)=delete;

    ~Bcf();

    bool isopen() const;

    const std::string fname_;
    htslib::htsFile *fid_;
    BcfHeader hdr_;
};

// @title Query the next record
// @param the pointer to open htslib file
// @param the pointer to the location in memory that data are loaded
// @param the format id for data to be loaded into record.
// @return: 0 for success otherwise non-zero
int next_record(Bcf* bid, BcfFloatRecord *rec, const char *id);

}

// @title: The number of values stored in format id
// @description: Each bcf format field is able to hold unique
//  number of values per sample.  This function provides a 
//  simple interface to the bcf file to retrieve that number.
// @param Rcpp::SEXP for an Rcpp::XPtr<Bcf> instantiated object
//  for the vcf file to be queried
// @param id: the format field id
// @return if an error occured that value returned is < 0, 
//  otherwise the number of values of fmt field id recorded per
//  sample is returned.
double k_fmt(Rcpp::XPtr<bcfio::Bcf> bid, const char* id);

// See htslib/vcf.h line 649
// Remember that n is the number of entries in the triplet of 
// dictionaries in the VCF.  BCF_DT_SAMPLE, provides the index of n
// that correspondes to the number of samples.
// size_t n_samples() const { return hdr_.n_samples(); };
// 
// int set_samples(const char *filename);
// 
// // TODO: sample_names
// const std::unique_ptr<std::string[]> sample_names() const { 
//     return hdr_.sample_names();
// }
// 


// mode according to htslib: quoting from htslib/hts.h line 608
//
//  @example
//      [rw]b  .. compressed BCF, BAM, FAI
//      [rw]bu .. uncompressed BCF
//      [rw]z  .. compressed VCF
//      [rw]   .. uncompressed VCF
//
// End quote
//
Rcpp::XPtr<bcfio::Bcf> bopen(const char* filename, const char* mode);

Rcpp::RObject query_next(Rcpp::XPtr<bcfio::Bcf> bid, const char* id);

#endif
