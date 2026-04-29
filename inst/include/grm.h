#ifndef GRM_H
#define GRM_H

#include <cctype>
#include <cstring>
#include <memory>
#include <chrono>

#include <bcfio.h>
// Grm class manages storage and access of GRM matrix
// 
// The GRM as an n_sample by n_sample symmetric, positive semi-definite
// matrix.  Let Z represent the n_sample by m_marker data genetic data.  
// From these data the GRM is computed as GRM = ZZ^T.
//
// @param n_samples of the GRM.
//

namespace grm {
// Symmetri matrix
// n is the number of columns in the matrix
constexpr uint64_t sym_matrix_idx_to_array(const uint64_t i, 
                        const uint64_t j, 
                        const uint64_t n) { 
    return i*n - i*(i-1)/2 + j - i; 
}

// n is the number of columns in the matrix
constexpr uint64_t matrix_idx_to_array(const uint64_t i, 
                              const uint64_t j,
                              const uint64_t n) { 
    return i*n + j; 
};

struct Grm {
    // 
    Grm();
    Grm(uint64_t n_samps);

    Grm(const Grm&)=delete;                          
    Grm& operator=(const Grm&)=delete;

    Grm(Grm&&);
    Grm& operator=(Grm&&);
                                            
    // Unchecked indexes when setting and getting of matrix values
    float operator()(const uint64_t i, const uint64_t j) const;
    float& operator()(const uint64_t i, const uint64_t j);

    // Checked indexes when setting and getting of matrix values
    int set(const uint64_t i, const uint64_t j, const float val); 
    int get(const uint64_t i, const uint64_t j, float *val) const; 

    uint64_t size() const;

    int midx_to_arr(const uint64_t i, 
            const uint64_t j, uint64_t* idx) const;

    uint64_t n_samples;
    std::unique_ptr<float[]> data;
};

int hap_update_kernel(Grm* grmat, const bcfio::BcfRecord<float>* rec);
}

#endif
