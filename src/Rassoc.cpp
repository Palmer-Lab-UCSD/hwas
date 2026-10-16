
#include <Rassoc.h>


CtgAssociations::~CtgAssociations() {
    if (recs_ != nullptr)
        delete[] recs_;

    if (contig_ != nullptr)
        delete[] contig_;

    recs_ = nullptr;
    contig_ = nullptr;
}


ctg_assoc_t CtgAssociations::init(const char* contig, uint32_t cap) {
    if (contig == nullptr)
        return nullptr;

    size_t contig_len = std::strlen(contig);
    char* new_contig = new(std::nothrow) char[contig_len + 1];
    if (new_contig == nullptr)
        return nullptr;
    std::memset(new_contig, '\0', contig_len+1);
    std::strncpy(new_contig, contig, contig_len);

    AssocRec* new_recs = new(std::nothrow) AssocRec[cap];
    if (new_recs == nullptr) {
        delete[] new_contig;
        return nullptr;
    }

    CtgAssociations* ctg = new(std::nothrow) CtgAssociations(new_contig, 
            new_recs, cap);
    if (ctg == nullptr) {
        delete[] new_contig;
        delete[] new_recs;
        return nullptr;
    }

    return ctg_assoc_t(ctg);
}


Status CtgAssociations::realloc_() {
    // As this is just resizing the data arrays and indexes
    // the contig c-string does not need to be updated

    uint32_t new_cap = cap_;
    while (new_cap <= size_) {

        if (new_cap > std::numeric_limits<uint32_t>::max() / 2)
            return Status::ErrOverflow;

        new_cap = new_cap * 2;
    }

    AssocRec* new_recs = new(std::nothrow) AssocRec[new_cap];
    if (new_recs == nullptr)
        return Status::ErrHeapAllocation;

    printf("Ralloc 2\n");
    for (uint32_t i = 0; i < cap_; i++)
        new_recs[i] = recs_[i];

    printf("Ralloc 3\n");
    delete[] recs_;
    recs_ = new_recs;
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

    if (size_ >= cap_) {
        printf("Reallocating\n");
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
    //
    AssocRec& rec = recs_[size_];
    rec.pos = pos;
    rec.df_numerator = ols_reduced.df - ols_full.df;
    rec.df_denom = ols_full.df;

    if (ols_reduced.df == ols_full.df) {
        rec.status = StatStatus::ErrFullAndReducedModelEq;
        rec.fstat = R_NaN;
        size_++;
        return Status::Success;
    }

    if (ols_full.df == 0) {
        rec.status = StatStatus::ErrFullModelSpansSampleSpace;
        rec.fstat = R_NaN;
        size_++;
        return Status::Success;
    }

    // TODO: over / underflow
    rec.fstat = (ols_reduced.sse - ols_full.sse) / rec.df_numerator;
    rec.fstat = rec.fstat * rec.df_denom / ols_full.sse;
    rec.status = StatStatus::Success;
    size_++;

    return Status::Success;
}


Rcpp::NumericMatrix CtgAssociations::to_r() const {
    Rcpp::NumericMatrix out(size_, 6);
    // TODO: Need to set column names
    
    Rcpp::CharacterVector cnames = Rcpp::CharacterVector::create("pos",
            "fstat",
            "df_numerator",
            "df_denom",
            "status",
            "pval");
    Rcpp::colnames(out) = cnames;

    // TODO: warning for when floating point can't faithfully
    // represent integers, e.g. when they exceed a certain number
    // the gaps between success values will be bigger than 1
    AssocRec* rec;
    for (uint32_t i = 0; i < size_; i++) {
        rec = &(recs_[i]);
        out(i, 0) = rec->pos;
        out(i, 1) = rec->fstat;
        out(i, 2) = rec->df_numerator;
        out(i, 3) = rec->df_denom;
        out(i, 4) = static_cast<int>(rec->status);
        out(i, 5) = R::pf(rec->fstat, 
                rec->df_numerator,
                rec->df_denom,
                false,
                false);
    }

    return out;
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
