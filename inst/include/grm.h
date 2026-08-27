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

enum struct Status : int {
    Success = 0
};

const char* status_msg(Status status);

// TODO: Double if constexpr is the right way to go.  I was getting
//  duplicate symboles linker errors without it.
// @title Compute index of one-dimension array from matrix indices
// @param i: matrix row index
// @param j: matrix col index, note that j >= i
// @param n: number of cols, and as symmetric number of rows, of the
//  matrix, note that n = max(i) + 1 = max(j) + 1
// @return index
constexpr uint32_t sym_matrix_idx_to_array(const uint32_t i, 
                        const uint32_t j, 
                        const uint32_t n) {
    return i*n - i*(i-1)/2 + j - i; 
};

// n is the number of columns in the matrix
constexpr uint32_t matrix_idx_to_array(const uint32_t i, 
                             const uint32_t j,
                             const uint32_t n) {
    return i*n + j; 
};

template <typename T>
struct Grm {
    Grm(): nsamps(0), capacity(0), data(nullptr) {};
    Grm(uint32_t n_samps) : 
        nsamps(n_samps),
        capacity(n_samps == 0 ? 0 : n_samps * (n_samps + 1) / 2),
        data(n_samps == 0 ? nullptr : new T[capacity]) {

        if (data != nullptr) {
            T default_val {};
            std::memset(data, default_val, capacity * sizeof(T));
        }
    }

    T operator()(const uint32_t i, const uint32_t j) const;

    T& operator()(const uint32_t i, const uint32_t j);

    int midx_to_arr(const uint32_t i, 
            const uint32_t j, 
            uint32_t* idx) const;

    uint32_t nsamps;        
    uint32_t capacity;      // size of allocated memory for data
    T* data;
};


template <typename T>
Grm<T>::Grm(uint32_t n_samps): nsamps(n_samps),
    capacity(n_samps == 0 ? 0 : n_samps * (n_samps + 1) / 2),
    data(nsamps == 0 ? nullptr : new T[capacity]) {

    if (data) {
        T default_val {};
        std::memset(data, default_val, nsamps * sizeof(T));
    }
}

// template <typename T>
// Grm<T>::Grm(Grm&& other)
//     : nsamps(other.nsamps),
//     data(std::move(other.data)) {
// 
//     other.nsamps = 0;    
//     grm::Grm& Grm::operator=(Grm&& other) {
//         if (this == &other)
//             return *this;
//     
//         nsamps = other.nsamps;
//         other.nsamps = 0;
//     
//         data = std::move(other.data);
//     
//         return *this;
//     }
// }

template <typename T>
Grm<T>::~Grm() {
    if (data)
        delete[] data;
    data = nullptr;
    nsamps = 0;
    capacity = 0;
}

template <typename T>
T Grm<T>::operator()(const uint32_t i, const uint32_t j) const {
    if (i > j)
        return data[sym_matrix_idx_to_array(j, i, nsamps)];
    return data[sym_matrix_idx_to_array(i, j, nsamps)];
}

template <typename T>
T& Grm<T>::operator()(const uint32_t i, const uint32_t j) {
    if (i > j)
        return data[sym_matrix_idx_to_array(j, i, nsamps)];

    return data[sym_matrix_idx_to_array(i, j, nsamps)];
}


template <typename T>
int Grm<T>::midx_to_arr(const uint32_t i, 
        const uint32_t j, 
        uint32_t* idx) const {

    if (i >= nsamps || j >= nsamps)
        return -1;

    // remember that by symmetry, the matrix is equal to its transpose
    if (i <= j)
        *idx = sym_matrix_idx_to_array(i, j, nsamps);
    else 
        *idx = sym_matrix_idx_to_array(j, i, nsamps);

    return 0;
}

// int hap_update_kernel(Grm* grmat, const bcfio::BcfRecord<float>* rec);

// compute_sum_grm();
// compute_norm_grm();

}

#endif
