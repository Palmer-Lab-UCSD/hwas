#ifndef RASSOC_H
#define RASSOC_H


// STL headers
#include <cstdint>
#include <cmath>
#include <cstring>

#include <memory>
#include <new>
#include <map>
#include <limits>


// R headers
#include <Rcpp.h>
#include <RcppEigen.h>

// HWAS headers
#include <status.h>
#include <bcfio.h>
#include <Rbcfio.h>



typedef Eigen::Map<Eigen::MatrixXd> mapmat_t;
typedef Eigen::Map<Eigen::VectorXd> mapvec_t;


enum struct StatStatus: int {
    Success                         = 0,
    ErrFullAndReducedModelEq        = -2,
    ErrFullModelSpansSampleSpace    = -3,
};


struct OlsResult {
    OlsResult(): sse(0), df(0) {};
    OlsResult(double sse_, int df_): sse(sse_), df(df_) {};
    double sse; // sum squared errors
    int df;   // degrees of freedom 
};


struct AssocRec {
    int64_t pos;
    double fstat;
    double df_numerator;
    double df_denom;
    StatStatus status;
};

struct CtgAssociations{
public:
    CtgAssociations()                                   = delete;
    CtgAssociations(const CtgAssociations&)             = delete;
    CtgAssociations(CtgAssociations&&)                  = delete;
    CtgAssociations& operator=(const CtgAssociations&)  = delete;
    CtgAssociations& operator=(CtgAssociations&&)       = delete;

    ~CtgAssociations();

    static std::unique_ptr<CtgAssociations> init(const char* contig,
            uint32_t cap);

    Status update(int64_t pos, 
            const OlsResult& ols_reduced,
            const OlsResult& ols_full);

    Rcpp::String contig_r_str() const;
    Rcpp::NumericMatrix to_r() const;

private:
    Status realloc_();

    CtgAssociations(const char* contig,
            AssocRec* recs,
            uint32_t cap) :
        contig_(contig),
        recs_(recs),
        cap_(cap),
        size_(0) {};

    const char* contig_;

    AssocRec* recs_;

    uint32_t cap_;
    uint32_t size_;
};

typedef std::unique_ptr<CtgAssociations> ctg_assoc_t;

struct Associations {
    Associations(uint32_t mpos_): mpos(mpos_) {};

    CtgAssociations* switch_contig(bcfio::Bcf* bid, uint32_t ctg_rid);
    Rcpp::List to_r();

    uint32_t mpos;
    std::map<uint32_t, ctg_assoc_t> results;
};

// @brief Compute the sum of squared errors under OLS
// @description Rencher and Schaalje's book "Linear Models in statistics
//  defines the sum squared errors (SSE) as
//
//  SSE = sum_i (y_i - x_i\hat{beta})^2
//      = y^Ty - \hat{\beta}^T X^T y
//
//  Which under the ordinary least squares (OLS) solution to \hat{\beta}
//
//  SSE = y^T (I - H) y
//
//  where
//
//  H = X (X^T X)^{-1} X^T 
//
//  is the idempotent projection matrix to the column space of X. 
//  I - H is then the idempotent matrix projecting the data to the
//  residual / error space of the data y.  
//
//  This understanding is important as it is the basis for the approach
//  in the ols_sse_ function.  The errors are computed by 
//
OlsResult ols_sse_(const Eigen::MatrixXd& x, 
        const Eigen::VectorXd& y,
        Eigen::VectorXd& yerr,
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd>&);


struct MatrixBlock {
    Eigen::Index ridx;
    Eigen::Index cidx;
    Eigen::Index rsize;
    Eigen::Index csize;
};


//
template <typename T>
void update_genetic_data_matrix_(Eigen::MatrixXd& gdata,
        const bcfio::BcfRecord<T>* brec) {

    T* data = brec->data;
    uint16_t kfounders = brec->ncol;
    uint16_t gdata_ncol = kfounders - 1;
    // uint16_t gdata_ncol = kfounders;

    for (uint32_t i = 0; i < brec->nrow; i++) {
        for (uint16_t j = 0; j < gdata_ncol; j++)
            gdata(i, j) = data[j];

        data += kfounders;
    }
}

// Here we load records, one-at-a-time, for association analysis.
// The resulting F - statistic and position are stored as C++
// STL vectors.  When storing records I rely on C++ std::vector
// deep copy of data.  Moreover, I assume that the contigs in
// the bcf are grouped together.
//
// @param bid is a pointer to the open bcf stream
// @param format_id is th format key for retrieving genetic data from
//  bcf file
// @param y is the transformed phenotypes
// @param x is 
template <typename T>
Rcpp::List run_stream_ctg_assoc_(bcfio::Bcf* bid,
        const char* format_id,
        const Eigen::VectorXd& y,
        const Eigen::MatrixXd& covars,
        const Eigen::MatrixXd& eig_transform) {


    // Prepare variables and run OLS on reduced model.  
    Eigen::Index nsamps = y.size();

    Eigen::VectorXd ytransformed = eig_transform * y;

    // The reduced model is an intercept + covariates.  Moreover,
    // the reduced model needs to be transformed into the diagonlized
    // basis using eig_transform
    Eigen::Index p_reduced = covars.cols() + 1;

    Eigen::MatrixXd x_reduced = Eigen::MatrixXd::Ones(nsamps, p_reduced);
    x_reduced.block(0, 1, nsamps, covars.cols()) = covars;
    x_reduced = eig_transform * x_reduced;

    // Allocate memory for error vector
    Eigen::VectorXd yerr(nsamps);

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_reduced(nsamps, p_reduced);
    OlsResult ols_reduced = ols_sse_(x_reduced, ytransformed, yerr, qr_reduced);


    // Prepare genetic data matrices.  The genetic data consists of two
    // matrices.  The first is gdata, it is the matrix that bcf records
    // are copied.  The second is xfull, which is the design matrix of
    // the full model.  The first set of columns are the intercept and
    // reduced model (transformed) and the second set are for genetic
    // data.  When the genetic data are loaded into gdata, they are
    // transformed by the matrix eig_transform and inserted into the xfull
    // design matrix using Eigen matrix blocks.  Block indices are
    // organized in the MatrixBlock struct.
    uint16_t kfeatures = 0;
    Status status = k_fmt(bid, format_id, &kfeatures);
    if (status != Status::Success)
        Rcpp::stop(status_msg(status));

    Eigen::Index gdata_ncols = static_cast<Eigen::Index>(kfeatures - 1);

    Eigen::Index p_full = p_reduced + gdata_ncols;
    MatrixBlock gblc { 0, p_reduced, nsamps, gdata_ncols };

    Eigen::MatrixXd xfull(nsamps, p_full);
    xfull.block(0, 0, nsamps, p_reduced) = x_reduced;

    Eigen::MatrixXd gdata(nsamps, gdata_ncols);

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_full(nsamps, p_full);

    // Prepare for reading records from bcf file
    bcfio::brec_t<T> brec = bcfio::BcfRecord<T>::init();
    if (brec == nullptr)
        Rcpp::stop(status_msg(Status::ErrBcfRecordInvalid));


    // Instantiate contig and position variables
    uint32_t m_pos = 100000;
    if (!bid->pos.empty())
        m_pos = static_cast<uint32_t>(bid->pos.size());

    Associations results { m_pos };
     
    OlsResult ols_full;
     
    status = bcfio::next_record<T>(bid, brec.get(), format_id);
    if (status != Status::Success)
        Rcpp::stop("Could not read bcf data");

    int32_t curr_ctg_rid = brec->rec->rid;
    int32_t prev_ctg_rid = curr_ctg_rid;
    CtgAssociations* curr_res = results.switch_contig(bid, curr_ctg_rid);
    if (curr_res == nullptr)
        Rcpp::stop("this should not be possible");

    // Run the analysis until we have tested all genomic locations
    // specified by the user.

    int64_t pos = 0;
    const uint64_t idx_report_val = 1000;
    REprintf("start association analysis\n");
    for (uint64_t i = 1;; i++) {

        if (i % idx_report_val == 0) {
            Rcpp::checkUserInterrupt();
            REprintf("Positions processed: %" PRIu64 "\n", i);
        }

        // need to +1 as bcf is zero based positions
        pos = brec->rec->pos + 1;
        curr_ctg_rid = brec->rec->rid;

        // Need to switch contigs
        if (prev_ctg_rid != curr_ctg_rid) {
            curr_res = results.switch_contig(bid, curr_ctg_rid);
            if (curr_res == nullptr)
                Rcpp::stop("this should not be possible");
        }

        // update 
        update_genetic_data_matrix_<T>(gdata, brec.get());
        if (status != Status::Success)
            Rcpp::stop(status_msg(status));

        // transform genetic data and store in xfull design matrix
        xfull.block(gblc.ridx, gblc.cidx, 
                gblc.rsize, gblc.csize) = eig_transform * gdata;

        ols_full = ols_sse_(xfull, ytransformed, yerr, qr_full);
        
        status = curr_res->update(pos, ols_reduced, ols_full);
        if (status != Status::Success)
            Rcpp::stop(status_msg(status));

        prev_ctg_rid = curr_ctg_rid;

        status = bcfio::next_record<T>(bid, brec.get(), format_id);
        if (status != Status::Success)
            break;
    }
    
    if (status != Status::EndOfFile)
        Rcpp::stop(status_msg(Status::ErrParseBcf)); 

    return results.to_r();
}



// if (brec == nullptr)
//     Rcpp::stop(status_msg(Status::ErrInternal));

// // Prepare design matrix and coefficient matrix
// uint16_t k = 0;
// Status status = bcfio::k_fmt(bid, format_id, &k)


// // If we can't read the first record abort
// Status status = bcfio::next_record(bid, brec.get(), format_id);
// if (status != Status::Success)
//     Rcpp::stop(status_msg(status));

// int32_t ctg_rid = brec->rec->rid;

// const char* ctg_name= nullptr;

// int64_t nrecords = 10000
// if (!bid->pos.empty())
//     nrecords = bid->pos.size();

// positions = std::vector<int64_t>(nrecords);
// fstat_vals = std::vector<double>(nrecords);
// fstat_val = double;
// // Resets size to 0 without affecting capacity.  The reason is that
// // for some unknown reason I can't figure out how to instantiate a 
// // vector with a specified capacity and not size.  
// pos.clear();
// fstat.clear();

// int32_t contig_rid = -1;
// for (; bstatus == Status::Success; idx++) {
//     tmp_rid = brec->rec->rid;


//     // Register that the contig has changed and that the
//     // data should be stored and in results map.  The vectors
//     // for holding data are reset
//     if (tmp_rid != contig_rid) {
//         ctg_name = bcf_hdr_id2name(bid->hdr, current_rid);

//         // move assignment
//         // note that argument data are all copied to the new instance
//         // of ContigAssocRecord stored at results[contig_rid]
//         results[contig_rid] = ContigAssocRecord(ctg_name,
//                 fstat_vals,
//                 positions);

//         fstat_vals.clear();
//         positions.clear();

//         // Terminate early as contigs aren't grouped together
//         if (results.count(tmp_rid) == 1)
//             break

//         current_rid = tmp_rid;
//     }

//     status = calc_assoc_at_locus_(transformed_design_matrix, 
//             transformed_phenotypes,
//             c_matrix,
//             &fstat);

//     if (status != Status::Success) {
//     }

//     if (idx % idx_report_val == 0) {
//         Rcpp::checkUserInterrupt();
//         REprintf("\r\x1B[K" "%" PRIu64 " positions processed.",
//                 idx);
//     }

//     bstatus = bcfio::next_record(bid, brec.get(), format_id);
// }


// List 


Rcpp::List stream_ctg_assoc_no_covs_(bconn_t bconn,
        const char* format_id,
        const Rcpp::NumericVector& phenotypes,
        const Rcpp::NumericMatrix& grm,
        const Rcpp::List& var_components);


// @brief Interface between C++ implementation and R function
// @description Users call stream_ctg_assoc function in the file R/assoc.R
//  to run association analysis.  The R function is responsible for basic
//  validation of user input, e.g. checking that the row and column names
//  of data identically enumerate sample names, and estimating variance
//  components.  The variance component estimates are done at the R level
//  as the run-time requirements are not as strict as are those for 
//  association analysis. This function is responsible for converting Rcpp
//  objects to Eigen data types for the actual computation. The
//  computation is run by run_stream_ctg_assoc_ C++ function.
// @param bconn is a Rxpp::XPtr to an open Bcf file
// @param phenotypes vector of sample phenotypes
// @param 
Rcpp::List stream_ctg_assoc_(bconn_t bconn,
        const char* format_id,
        const Rcpp::NumericVector& phenotypes,
        const Rcpp::NumericMatrix& covariates,
        const Rcpp::NumericMatrix& grm,
        const Rcpp::List& var_components);


Rcpp::List assoc_no_covars_(const Rcpp::NumericMatrix&,
        const Rcpp::NumericVector& phenotypes,
        const Rcpp::NumericMatrix& grm,
        const Rcpp::List& var_components);

#endif
