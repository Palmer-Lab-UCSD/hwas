
#include <Rassoc.h>


CtgAssociations::~CtgAssociations() {
    if (fstats_ != nullptr)
        delete[] fstats_;

    if (pvals_ != nullptr)
        delete[] pvals_;

    if (pos_ != nullptr)
        delete[] pos_;

    if (status_ != nullptr)
        delete[] status_;

    if (contig_ != nullptr)
        delete[] contig_;

    fstats_ = nullptr;
    pvals_ = nullptr;
    pos_ = nullptr;
    status_ = nullptr;
    contig_ = nullptr;
}


ctg_assoc_t CtgAssociations::init(const char* contig, uint32_t cap) {
    if (contig == nullptr)
        return nullptr;

    size_t contig_len = std::strlen(contig);
    char* new_contig = new(std::nothrow) char[contig_len + 1];
    if (new_contig == nullptr)
        return nullptr;
    std::strncpy(new_contig, contig, contig_len);

    double* fstats = new(std::nothrow) double[cap];
    if (fstats == nullptr) {
        delete[] new_contig;
        return nullptr;
    }

    double* pvals = new(std::nothrow) double[cap];
    if (pvals == nullptr) {
        delete[] new_contig;
        delete[] fstats;
        return nullptr;
    }

    int64_t* pos = new(std::nothrow) int64_t[cap];
    if (pos == nullptr) {
        delete[] new_contig;
        delete[] fstats;
        delete[] pvals;
        return nullptr;
    }

    StatStatus* status = new(std::nothrow) StatStatus[cap];
    if (status == nullptr) {
        delete[] new_contig;
        delete[] fstats;
        delete[] pvals;
        delete[] pos;
        return nullptr;
    }


    CtgAssociations* ctg = new(std::nothrow) CtgAssociations(new_contig, 
            fstats, pvals, pos, status, cap);
    if (ctg == nullptr) {
        delete[] new_contig;
        delete[] fstats;
        delete[] pvals;
        delete[] pos;
        delete[] status;
        return nullptr;
    }

    return ctg_assoc_t(ctg);
}




Status CtgAssociations::realloc_() {
    // As this is just resizing the data arrays and indexes
    // the contig c-string does not need to be updated
    if (cap_ > std::numeric_limits<uint32_t>::max() / 2)
        return Status::ErrOverflow;

    uint32_t new_cap = 2 * cap_;
    double* new_fstats = new(std::nothrow) double[new_cap];
    if (new_fstats == nullptr)
        return Status::ErrHeapAllocation;

    double* new_pvals = new(std::nothrow) double[new_cap];
    if (new_pvals == nullptr) {
        delete[] new_fstats;
        return Status::ErrHeapAllocation;
    }

    int64_t* new_pos = new(std::nothrow) int64_t[new_cap]; 
    if (new_pos == nullptr) {
        delete[] new_fstats;
        delete[] new_pvals;
        return Status::ErrHeapAllocation;
    }

    StatStatus* new_status = new(std::nothrow) StatStatus[new_cap];
    if (new_status == nullptr) {
        delete[] new_fstats;
        delete[] new_pvals;
        delete[] new_pos;
        return Status::ErrHeapAllocation;
    }

    for (uint32_t i = 0; i < cap_; i++) {
        new_fstats[i] = fstats_[i];
        new_pvals[i] = pvals_[i];
        new_pos[i] = pos_[i];
        new_status[i] = status_[i];
    }

    delete[] fstats_;
    fstats_ = new_fstats;

    delete[] pvals_;
    pvals_ = new_pvals;

    delete[] pos_;
    pos_ = new_pos;

    delete[] status_;
    status_ = new_status;

    size_ = cap_;
    cap_ = new_cap;

    return Status::Success;
}


Rcpp::String CtgAssociations::contig_r_str() const {
    return Rcpp::String(contig_);
}


Status CtgAssociations::update(int64_t pos, 
        const OlsResult& ols_reduced,
        const OlsResult& ols_full) {

    Status status = Status::ErrInvalidInput;
    if (ols_reduced.df < ols_full.df)
        return Status::ErrInvalidInput;

    if (size_ == cap_) {
        status = realloc_();
        if (status != Status::Success)
            return status;
    }

    // Cases:
    //
    // The full model does not contain any information not in covariates.
    // This is detected by the degrees-of-freedom of both models being
    // equivalent.  Claude states that NaN is the appropriate result as
    // we can not attribute the association to either genetics or covariates.
    if (ols_reduced.df == ols_full.df) {
        pvals_[size_] = R_NaN;
        pos_[size_] = pos;
        status_[size_] = StatStatus::ErrFullAndReducedModelEq;
        size_++;
        return Status::Success;
    }


    if (ols_full.df == 0) {
        pvals_[size_] = R_NaN;
        pos_[size_] = pos;
        status_[size_] = StatStatus::ErrFullModelSpansSampleSpace;
        size_++;
        return Status::Success;
    }

    // numerator contribution to f-statistic
    int num_df = ols_reduced.df - ols_full.df;
    double fstat = (ols_reduced.sse - ols_full.sse);
    fstat = fstat / num_df;

    // denominator contribution to f-statistic
    fstat = fstat / ols_full.sse;
    fstat = fstat * ols_full.df;

    fstats_[size_] = fstat;

    // Compute P(F > fstat) under the central F distribution
    pvals_[size_] = R::pf(fstat, num_df, ols_full.df, true, false);
    pos_[size_] = pos;
    status_[size_] = StatStatus::Success;
    size_++;

    return Status::Success;
}


Rcpp::List CtgAssociations::to_r() const {
    Rcpp::NumericVector out_fstats(size_);
    Rcpp::NumericVector out_pvals(size_);

    // annoyingly R's integers are limited to 32 bit meaning that
    // the maximum value is about 2.1 billion.  To be safe, I will
    // just store as double
    Rcpp::NumericVector out_pos(size_);

    Rcpp::IntegerVector out_status(size_);
    Rcpp::CharacterVector out_status_msg(size_);

    for (uint32_t i = 0; i < size_; i++) {
        out_fstats[i] = fstats_[i];
        out_pvals[i] = pvals_[i];
        out_pos[i] = static_cast<double>(pos_[i]);

        out_status[i] = static_cast<int>(status_[i]);
        switch (status_[i]) {
        case StatStatus::Success:
            out_status_msg[i] = "success";
            break;
        case StatStatus::ErrFullAndReducedModelEq:
            out_status_msg[i] = "failed: full and reduce model equal df";
            break;
        case StatStatus::ErrFullModelSpansSampleSpace:
            out_status_msg[i] = "failed: full model df = 0";
            break;
        default:
            out_status_msg[i] = "unknown";
        }
    }

    return Rcpp::List::create(Rcpp::Named("contig") = contig_,
            Rcpp::Named("fstats") = out_fstats,
            Rcpp::Named("pos") = out_pos,
            Rcpp::Named("pvals") = out_pvals,
            Rcpp::Named("status") = out_status,
            Rcpp::Named("status_msg") = out_status_msg);
}


CtgAssociations* Associations::switch_contig(bcfio::Bcf* bid,
        uint32_t ctg_rid) {
    // if contig is not already apart of the map, add it
    if (results.count(ctg_rid) == 0) {
        const char* ctg_name = bcf_hdr_id2name(bid->hdr, 
                ctg_rid);

        ctg_assoc_t res_ = CtgAssociations::init(ctg_name, mpos);
        if (res_ == nullptr)
            return nullptr;

        results[ctg_rid] = std::move(res_);
    }

   
    // This value may be nullptr, need to check
    return results[ctg_rid].get();
}


Rcpp::List Associations::to_r() {
    Rcpp::String ctg_name;
    Rcpp::List out {};
    for (auto& [key, ptr] : results) {
        ctg_name = ptr->contig_r_str();
        out[ctg_name] = ptr->to_r();
    }

    return out;
}

// @brief compute the ordinary least squares sum squared errors
// @param x Is the n sample x p predictor design matrix
// @param y Is the n sample column vector of sample responses
// @param qr Is the pre-allocated n sample x n sample qr matrix object
// @param yerr Is the pre-allocated n sample column vector for 
//  sample errors
// @return the sum squared errors of the linear model
OlsResult ols_sse_(const Eigen::MatrixXd& x,
        const Eigen::VectorXd& y,
        Eigen::VectorXd& yerr,
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& qr) {

    qr.compute(x); 
    int rank = static_cast<int>(qr.rank());
    int nsamps = static_cast<int>(y.size());

    auto qmat = qr.householderQ();  
    yerr = qmat.transpose() * y;
    yerr.head(rank).setZero();

    yerr = qmat * yerr;

    return OlsResult(yerr.squaredNorm(), nsamps - rank);
}


Rcpp::List stream_ctg_assoc_no_covs_(bconn_t bconn, 
        const char* format_id,
        const Rcpp::NumericVector& phenotypes,
        const Rcpp::NumericMatrix& grm,
        const Rcpp::List& var_components) {
    Rcpp::stop("not implemented error");
}


// [[Rcpp::export]]
Rcpp::List stream_ctg_assoc_(bconn_t bconn, 
        const char* format_id,
        const Rcpp::NumericVector& phenotypes,
        const Rcpp::NumericMatrix& covariates,
        const Rcpp::NumericMatrix& grm,
        const Rcpp::List& var_components) {

    if (!bcfio::is_open(bconn.get()))
        Rcpp::stop(status_msg(Status::ErrBcfNotOpen));

    if (format_id == nullptr)
        Rcpp::stop(status_msg(Status::ErrInvalidInput));

    // Eigen decompose the phenotype covariance matrix under the model
    const mapmat_t grm_mat = Rcpp::as<mapmat_t>(grm);

    Eigen::Index nsamps = grm_mat.rows();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_decomp(nsamps);
    eig_decomp.compute(grm_mat);
    if (eig_decomp.info() != Eigen::Success)
        Rcpp::stop("Eigendecomposition of grm failed.");

    double delta = var_components["delta"];

    const Eigen::VectorXd& grm_eig_vals = eig_decomp.eigenvalues();
    Eigen::MatrixXd eigval_inv_matrix = Eigen::MatrixXd::Zero(nsamps, nsamps);

    // TODO: Make more robust for when delta is zero, or Inf.
    for (int i = 0; i < nsamps; i++)
        eigval_inv_matrix(i, i) = 1 / std::sqrt(grm_eig_vals[i] + delta);


    Eigen::MatrixXd eig_transform = eigval_inv_matrix * eig_decomp.eigenvectors().transpose();

    const mapvec_t y = Rcpp::as<mapvec_t>(phenotypes);
    const mapmat_t covars = Rcpp::as<mapmat_t>(covariates);

    // Determine the type data are stored in the bcf stream
    bcfio::BcfHdrAttr hattr {};
    Status status = bcfio::decode_hts_idinfo(bconn->hdr,
            format_id,
            BCF_HL_FMT,
            &hattr);
    if (status != Status::Success)
        Rcpp::stop(status_msg(status));

    switch (hattr.type) {
    case BCF_HT_REAL:
        return run_stream_ctg_assoc_<float>(bconn.get(),
            format_id,
            y, 
            covars, 
            eig_transform);
    case BCF_HT_INT:
        Rcpp::stop("Not implemented error");
    }

    Rcpp::stop("Unsupported Bcf Type"); 
}


 //   Eigen::VectorXd y_transformed = eig_vecs_transpose * y;
//    x_reduced = eig_vecs_transpose * x_reduced;
//     // Prepare the reduced model
//     int p_reduced = covariates.ncol() + 1;
//     Eigen::MatrixXd x_reduced = Eigen::MatrixXd::Ones(nsamps, p_reduced);
//     // Add the covariate values to the design matrix with intercept
//     for (int i = 0; i < nsamps; i++) {
//         for (int j = 1; j < p_reduced; j++)
//             x_reduced(i, j) = covariates(i, j-1);
//     }
// 
