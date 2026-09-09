#ifndef RBCFIO_H
#define RBCFIO_H


// STL C DEPENDENCIES
#include <cassert>
#include <cstdio>
#include <cstddef>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <cstdint>
#include <csignal>
#include <cinttypes>

// STL C++ DEPENDENCIES
#include <optional>
#include <memory>
#include <utility>
#include <map>

// RCPP INDEPENDENT HEADERS
#include <bcfio.h>

// RCPP AND RCPP DEPENDENT HEADERS
#include <Rcpp.h>


typedef Rcpp::XPtr<bcfio::Bcf> bconn_t;


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
// 
// See htslib/vcf.h line 649
// Remember that n is the number of entries in the triplet of 
// dictionaries in the VCF.  BCF_DT_SAMPLE, provides the index of n
// that correspondes to the number of samples.
// size_t n_samples() const { return hdr_.n_samples(); };
// 
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
bconn_t bread(const char* filename);
int bclose(bconn_t bconn);
bool is_open(const bconn_t bconn);

int64_t num_positions(bconn_t bconn);
uint32_t num_samples(bconn_t bconn);
uint16_t k_fmt(bconn_t bconn, const char* format_id);


Rcpp::RObject sample_names(bconn_t bconn);

// @brief Subset bconn sample list 
// @param[in] bconn is a connection to an open bcf file
// @param[in] samples is a vector character vector 
// @return number of samples upon success and < 0 upon error
int subset_samples(bconn_t bconn, 
        Rcpp::Nullable<Rcpp::CharacterVector> nullable_samples);

// @brief Subset bconn sample list using a sample inclusion file
int subset_samples_from_file(bconn_t bconn, const char* sample_filename);

// @brief register positions to retrieve data from file stream
// @return the number of positions registered
int subset_pos_from_file(bconn_t bconn, const char* filename);

int set_threads(bconn_t bconn, int n);


// @brief Get position-wide data and sample data matrix at current pos
// @param[in,out] bconn the bcf file handle in which data are read.
//  The connection object is updated to point to the data at the next 
//  genomic position in the bcf.
// @param[in] id is the bcf FORMAT field id.  It specifies the
//  subset of data from the record to load.
// @return R_NilValue upon end of file, or
//  n sample by k format values Rcpp::NumericMatrix and attached
//  attributes with position wide data:
//  * contig: a string, the contig name, usually chromosome name
//  * pos: integral type, the 1 based genomic position on contig
Rcpp::Nullable<Rcpp::NumericMatrix> next_record(bconn_t bconn, 
        const char* id);


// @brief gets a bcf record and stores data in NumericMatrix 
// @description this function is used by next_record.  It is needed
//  as records can be stored in distinct types.  To abstract the
//  types from end R user, this function is called by next_record
//  function to load the type specific data into an 
//  Rcpp::NumericMatrix.
// @param[in,out] bid the bcf file handle in which data are read. 
//  The file handle is mutated, it will 
// @param[in] id is the bcf FORMAT field id.  It specifies the
//  subset of data from the record to load.
// @return see next_record documentation 
template <typename T>
Rcpp::Nullable<Rcpp::NumericMatrix> get_record_matrix_(bcfio::Bcf* bid, 
        const char* id) {

    bcfio::brec_t<T> brec = bcfio::BcfRecord<T>::init();
    bcfio::Status status = bcfio::next_record<T>(bid, brec.get(), id);

    if (status == bcfio::Status::EndOfFile)
        return R_NilValue;

    if (status != bcfio::Status::Success)
        Rcpp::stop(bcfio::status_msg(status));

    Rcpp::NumericMatrix data(brec->nrow, brec->ncol);
    uint16_t ncol = brec->ncol;
    uint32_t nrow = brec->nrow;
    T* cur_samp_rec = brec->data;
    for (uint32_t i = 0; i < nrow; i++) {

        for (uint16_t j = 0; j < ncol; j++)
            data(i, j) = cur_samp_rec[j];

        cur_samp_rec += ncol;
    }

    const char* chr = bcfio::chrom(bid, brec.get());
    if (chr == nullptr) 
        Rcpp::stop(bcfio::status_msg(bcfio::Status::ErrInternal));
    data.attr("contig") = Rcpp::String(chr);

    int64_t p = -1;
    status = bcfio::pos(brec.get(), &p);
    if (status != bcfio::Status::Success)
        Rcpp::stop(bcfio::status_msg(status));

    data.attr("pos") = p;

    //TODO: data.attr("qual") = rec.qual();
    //TODO: ref / alt allales

    return data;
}


// TODO: enumerate info keys and be able to query values
// TODO: enumerate format keys and be able to query values


#endif
