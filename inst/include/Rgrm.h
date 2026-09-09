#ifndef RGRM_H
#define RGRM_H

// TODO clean up include declarations, only include what i need

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
#include <grm.h>

// RCPP AND RCPP DEPENDENT HEADERS
#include <Rcpp.h>
#include <Rbcfio.h>


// @brief compute the internal tool that computes grm
// @description Grm's may have distinct numeric types, and I don't 
//  want to burdern the R user with this fact.  The result is that
//  calc_unnormalized_grm is the user interface to the computation
//  while this function implements according to the type of the 
//  target data.  It is not my intention that this function be
//  called directly. 
// @param instance of bcf file handle 
// @param format id for measurment to use for grm
// @return Rcpp::NumericMatrix or Rcpp::stop error thrown
template <typename T>
Rcpp::NumericMatrix calc_unnormalized_grm_(bcfio::Bcf* bid,
        const char* id) {
    
    uint32_t nsamps = 0;
    bcfio::Status bstatus = bcfio::num_samples(bid, &nsamps);
    if (bstatus != bcfio::Status::Success)
        Rcpp::stop(bcfio::status_msg(bstatus));

    bcfio::brec_t brec = bcfio::BcfRecord<T>::init();
    if (brec == nullptr)
        Rcpp::stop("Internal Error, could not initialize bcf record"
                " for reading bcf data.");

    grm::ugrm_t g = grm::UnnormalizedGrm<T>::init(nsamps);
    if (g == nullptr)
        Rcpp::stop("Internal Error,k could not initialize grm.");
    
    grm::Status gstatus = grm::Status::ErrInternal;

    const uint64_t idx_report_val = 1000;

    bstatus = next_record(bid, brec.get(), id);
    for (uint64_t idx = 1; bstatus == bcfio::Status::Success; idx++) {

        gstatus = g->update(brec.get());

        if (gstatus != grm::Status::Success)
            Rcpp::stop("Internal Error in updating Grm");

        if (idx % idx_report_val == 0) {
            Rcpp::checkUserInterrupt();

            // Update the number of records processed for user
            //
            // recall that \x1b[ is, according to wikipedia ANSI escape
            // code article, is the start of an escape sequence.  The
            // designation of hex encoding \x is specific to how format
            // strings are parsed.  The hex number 0x1b encodes ASCII 
            // <ESC>, in decimal is 27. The control sequence specifies K,
            // where K is the ANSI control sequence to erase in line.
            // K can be preceeded by a number, but no necessary for our
            // purpose.  Together the format string below is read as
            // put the cursor at the beginning of the line (carriage
            // return, \r), erase contents from cursor to end of line
            // (ANSI escape sequence), and then print string starting
            // with a 64 bit unsigned int.  Note: that using the ANSI
            // escape sequence in this fashion was first recommended by
            // Google AI overview and then manually checked by background
            // reading.
            REprintf("\r\x1B[K" "%" PRIu64 " positions processed.",
                    idx);
        }

        bstatus = next_record(bid, brec.get(), id);
    }
    REprintf("\n");

    if (bstatus != bcfio::Status::EndOfFile)
        Rcpp::stop("Matrix computation did not complete");

    // I know this, strictly, isn't necessary, I justj wanted to make
    // clear that nsamps is an invariant of the grm
    nsamps = g->nsamples();
    Rcpp::NumericMatrix data(nsamps);
    for (uint32_t i = 0; i < nsamps; i++) {
        for (uint32_t j = 0; j < nsamps; j++) {
            // using unsafe retrieval of grm data from g. This is ok
            // nsamps is an invariant shared by data and grm
            data(i, j) = (*g)(i, j);
        }
    }

    return data;
}


// @brief compute the unnormalized grm
// @param bconn An open connection to a bcf file
// @param the "format" field of a bcf record id or key.  This C-string
//  specifies what data from a record to retrieve.
// @return grm as an Rcpp::NumericMatrix or R_NilValue upon success
//  and failure, respectively.
Rcpp::NumericMatrix calc_unnormalized_grm(bconn_t bconn, const char* id);



#endif
