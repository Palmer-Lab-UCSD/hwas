// Palmer Lab at UCSD 2026
//
#include <hwas_types.h>


////////////////////////////////////////////////////////////////////
// GRM CLASS
////////////////////////////////////////////////////////////////////
//
// Recall that the GRM is a symmetric matrix, therefore we only need
// to store the upper triagonal and diagonal element values.  
// Consequently, the size of the array storing the data is n*(n +1)/2.
//

grm::Grm::Grm(): n_samples(0), data(nullptr) {};

grm::Grm::Grm(uint64_t n_samps)
    : n_samples(n_samps),
        data(size() != 0 ? std::make_unique<float[]>(size()) : nullptr) {

    if (data)
       std::memset(data.get(), 0, size() * sizeof(float));
}

grm::Grm::Grm(grm::Grm&& other)
    : n_samples(other.n_samples), data(std::move(other.data)) {
    other.n_samples = 0;    
};

grm::Grm& grm::Grm::operator=(grm::Grm&& other) {
    if (this == &other)
        return *this;

    n_samples = other.n_samples;
    other.n_samples = 0;

    data = std::move(other.data);

    return *this;
}

// The number of upper diagonal + diagonal elements of the GRM
uint64_t grm::Grm::size() const { 
    return n_samples * (n_samples + 1) / 2; 
};


int grm::Grm::midx_to_arr(const uint64_t i, const uint64_t j,
        uint64_t* idx) const {

    if (i >= n_samples || j >= n_samples)
        return -1;

    // remember that by symmetry, the matrix is equal to its transpose
    if (i <= j)
        *idx = SYMMATRIX_IDX_TO_ARRAY(i, j, n_samples);
    else 
        *idx = SYMMATRIX_IDX_TO_ARRAY(j, i, n_samples);

    return 0;
}


float grm::Grm::operator()(const uint64_t i, const uint64_t j) const {
    if (i > j)
        return data[SYMMATRIX_IDX_TO_ARRAY(j, i, n_samples)];
    return data[SYMMATRIX_IDX_TO_ARRAY(i, j, n_samples)];
}


float& grm::Grm::operator()(const uint64_t i, const uint64_t j) {
    if (i > j)
        return data[SYMMATRIX_IDX_TO_ARRAY(j, i, n_samples)];
    return data[SYMMATRIX_IDX_TO_ARRAY(i, j, n_samples)];
}


int grm::Grm::get(const uint64_t i, const uint64_t j, float *val) const {
    uint64_t idx = 0;
    int status = -1;
    if ((status = midx_to_arr(i, j, &idx)) != 0) 
        return status;

    *val = data[idx];

    return status;
}


int grm::Grm::set(const uint64_t i, const uint64_t j, const float val) {
    uint64_t idx = 0;
    int status = -1;
    if ((status = midx_to_arr(i, j, &idx)) != 0) 
        return status;

    data[idx] = val;

    return status;
}


////////////////////////////////////////////////////////////////////
// GRM COMPUTATION
////////////////////////////////////////////////////////////////////

// @return -1 on error and 0 or success
static int update(grm::Grm* grmat, bcfio::BcfFloatRecord* rec) {

    // instantiate indexing variables used in for loops
    // use static to prevent construction and destruction of variables
    // between function calls
    uint64_t i = 0;
    uint64_t j = 0;
    uint64_t k_hap = 0;

    float val = 0;
    uint64_t k_haps = rec->ncols();
    uint64_t n_samples = rec->nrows();
    if (n_samples != grmat->n_samples)
        return -1;

    // Remember we are using pointer over the n individual
    // by k haplotype data matrix
    const float* rec_arr = rec->array();

    // get the address of the first position of data
    const float* grm_arr = grmat->data.get();

    float* samp_i = nullptr;
    float* samp_j = nullptr;

    // only iterate over upper triangle
    // remember that record data is an n_sample by k haplotype matrix 
    for (i = 0; i < n_samples; i++) {

        samp_i = rec_arr + MATRIX_IDX_TO_ARRAY(i, 0, k_haps);

        for (j = i; j < n_samples; j++) {

            val = 0;
            samp_j = rec_arr + MATRIX_IDX_TO_ARRAY(j, 0, k_haps);

            for (k_hap = 0; k_hap < k_haps; k_hap++)
                val += samp_i[k_hap] * samp_j[k_hap];

            grm_arr[SYMMATRIX_IDX_TO_ARRAY(i, j, n_samples)] += val;
        }
    }

    return 0;
}


// [[Rcpp::export]]
Rcpp::RObject calc_grm(Rcpp::XPtr<bcfio::Bcf> bid, const char* id) {

    // instantiate matrices to hold calculations
    int32_t k { 0 };
    if ((k = k_fmt(bid, id)) < 0) {
        return R_NilValue;
    }
    uint64_t nsamps = num_samples(bid);

    grm::Grm grmat { nsamps };
    bcfio::BcfFloatRecord rec {};

    size_t rec_count = 0;
    while (bcfio::next_record(bid.get(), &rec, id) == 0) {

        if (update(&grmat, &rec) != 0) {
            return R_NilValue;
        }

        if (++rec_count % 1000 == 0)
            Rprintf("Processed %zu records\n", rec_count);
    }



    Rcpp::NumericMatrix grmatrix(nsamps, nsamps);
    // fill in lower diagonal elements
    for (uint64_t i = 0; i < nsamps; i++)
        for (uint64_t j = 0; j < nsamps; j++)
            grmatrix(i, j) = grmat(j, i);
            
    Rprintf("Done\n");
    return grmatrix;
}
