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
Rcpp::Nullable<Rcpp::XPtr<bcfio::Bcf>> bopen(const char* filename, const char* mode);
int bclose(Rcpp::XPtr<bcfio::Bcf> bid);

Rcpp::RObject next_record(Rcpp::XPtr<bcfio::Bcf> bid, const char* id);

uint32_t num_samples(Rcpp::XPtr<bcfio::Bcf> bid);

Rcpp::RObject sample_names(Rcpp::XPtr<bcfio::Bcf>);

int subset_samples(Rcpp::XPtr<bcfio::Bcf> bid, const char* filename);
int set_threads(Rcpp::XPtr<bcfio::Bcf> bid, int n);
double k_fmt(Rcpp::XPtr<bcfio::Bcf> bid, const char* format_id);


/////////////////////////////////////////////////////////////////////
/// GRM
/////////////////////////////////////////////////////////////////////


// @param instance of bcf file handle 
// @param format id for measurment to use for grm
// @return R_NilValue if err otherwise an n sample by n sample
//  Rcpp::NumericMatrix
Rcpp::RObject calc_grm(Rcpp::XPtr<bcfio::Bcf> bid, const char* id);



#endif
