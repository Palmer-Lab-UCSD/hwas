//
//
#ifndef HWAS_H
#define HWAS_H

#include <cassert>
#include <cstdio>
#include <cstddef>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <cstdint>
#include <csignal>

#include <optional>
#include <memory>
#include <utility>

#include <Rcpp.h>

#include <bcfio.h>
#include <grm.h>


// TODO: SUBSET BY SAMPLES,
// TODO: SPECIFY POSITION SPEC
// TODO: GET R FRONT END FRO ASSOCATION AND 
//  COMPUTE FOR A SINGLE LOCUS
//


/////////////////////////////////////////////////////////////////////
/// BCFIO
/////////////////////////////////////////////////////////////////////
///
typedef Rcpp::XPtr<bcfio::Bcf> bcf_conn_t;


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
bcf_conn_t bread(const char* filename);
int bclose(bcf_conn_t bid);


int64_t num_positions(bcf_conn_t bid);
uint32_t num_samples(bcf_conn_t bid);
uint16_t k_fmt(bcf_conn_t bid, const char* format_id);


Rcpp::RObject sample_names(bcf_conn_t bid);

int subset_samples(bcf_conn_t bconn, Rcpp::CharacterVector samples);
int subset_samples_from_file(bcf_conn_t bconn, const char* sample_filename);

// @brief register positions to retrieve data from file stream
// @return the number of positions registered
int subset_pos_from_file(bcf_conn_t bconn, const char* filename);

int set_threads(bcf_conn_t bid, int n);


Rcpp::Nullable<Rcpp::NumericMatrix> next_record(bcf_conn_t bid, 
        const char* id);


// return R_NilValue upon end of file
//  matrix for computation
template <typename T>
Rcpp::Nullable<Rcpp::NumericMatrix> get_matrix(bcfio::Bcf* bid, 
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
    data.attr("chrom") = Rcpp::String(chr);

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



/////////////////////////////////////////////////////////////////////
/// GRM
/////////////////////////////////////////////////////////////////////


// @param instance of bcf file handle 
// @param format id for measurment to use for grm
// @return R_NilValue if err otherwise an n sample by n sample
//  Rcpp::NumericMatrix
Rcpp::RObject calc_grm(bcf_conn_t bid, const char* id);



/////////////////////////////////////////////////////////////////////
/// PG_SIM
/////////////////////////////////////////////////////////////////////
///

//Rcpp::NumericVector pg_sim(const char* bcf_filename,
//        float heritability);

#endif
