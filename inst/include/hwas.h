//
//
#ifndef HWAS_H
#define HWAS_H

#include <cstdio>
#include <cstddef>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <cstdint>
#include <csignal>

#include <optional>
#include <memory>
#include <string>
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
bcf_conn_t bopen(const char* filename, const char* mode);
int bclose(bcf_conn_t bid);


int64_t num_positions(bcf_conn_t bid);
uint32_t num_samples(bcf_conn_t bid);
uint16_t k_fmt(bcf_conn_t bid, const char* format_id);


Rcpp::RObject sample_names(bcf_conn_t bid);

int subset_samples(bcf_conn_t bid, Rcpp::CharacterVector samples);
int subset_samples_from_file(bcf_conn_t bid, const char* sample_filename);

int set_threads(bcf_conn_t bid, int n);

Rcpp::RObject next_record(bcf_conn_t bid, const char* id);

template <typename T>
Rcpp::NumericMatrix get_matrix(bcfio::Bcf* bid, 
        const char* id,
        const htslib::bcf_fmt_t* fmt_cfg) {

    std::unique_ptr<bcfio::BcfRecord<T>> brec =
        std::make_unique<bcfio::BcfRecord<T>>();
    uint32_t nrow = 0;
    uint16_t ncol = 0;

    if (bcfio::num_samples(bid, &nrow) != 0) {
        fprintf(stderr,
                "ERROR: parse record errror\n");
        return R_NilValue;
    }

    if (bcfio::k_fmt(bid, id, &ncol) != 0) {
        fprintf(stderr,
                "ERROR: parse record errror\n");
        return R_NilValue;
    }

    if (ncol != fmt_cfg->n) {
        fprintf(stderr, 
                "ERROR: expected values per sample %d got %d\n",
                ncol,
                fmt_cfg->n);
        return R_NilValue;
    }

    if (bcfio::next_record<T>(bid, brec.get(), id) != 0)
        return R_NilValue;

    Rcpp::NumericMatrix data(nrow, ncol);
    T* cur_samp_rec = nullptr;
    for (uint32_t i = 0; i < nrow; i++) {

        cur_samp_rec = brec->data + i;

        for (uint16_t j = 0; j < ncol; j++)
            data(i, j) = cur_samp_rec[j];
    }

    data.attr("chrom") = bcfio::chrom(bid, brec.get());
    int64_t p = -1;
    bcfio::pos(brec.get(), &p);
    data.attr("pos") = p;

    //data.attr("qual") = rec.qual();
    //data.attr(

    return data;
}


// TODO: get ref allele
// TODO: get alt allele
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

Rcpp::NumericVector pg_sim(const char* bcf_filename,
        float heritability);

#endif
