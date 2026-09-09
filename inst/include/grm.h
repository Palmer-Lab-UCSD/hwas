#ifndef GRM_H
#define GRM_H

#include <cctype>
#include <cstring>

#include <memory>
#include <limits>
#include <new>

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
    Success                         = 0,
    ErrDimensionsNotEqual           = -2,
    ErrIndexOutofBounds             = -3,
    ErrInternal                     = -4
};

const char* status_msg(Status status);


template <typename T>
class UnnormalizedGrm {
public:
    UnnormalizedGrm<T>()                                        = delete;
    UnnormalizedGrm<T>(const UnnormalizedGrm<T>&)               = delete;
    UnnormalizedGrm<T>(UnnormalizedGrm<T>&&)                    = delete;
    UnnormalizedGrm<T>& operator=(const UnnormalizedGrm<T>&)    = delete;
    UnnormalizedGrm<T>& operator=(UnnormalizedGrm<T>&&)         = delete;

    ~UnnormalizedGrm<T>();
    // T operator()(const uint32_t i, const uint32_t j) const;
    
    // @brief factory function for new UnnormalizedGrm<T>
    // @param[in] n the number of samples for the grm
    // @return unique_ptr to initialized Grm
    static std::unique_ptr<UnnormalizedGrm<T>> init(uint32_t n);

    // @brief Update the UnnormalizedGrm<T> data with a BcfRecord
    // @description
    // @param[in] brec is a record queried from a Bcf file
    // @return Status:  Success, ErrDimensionsNotEqual, 
    Status update(const bcfio::BcfRecord<T>* brec);

    // @brief unsafe retrieval of values
    T operator()(const uint32_t i, const uint32_t j) const;

    // @brief safe get and set of values in grm
    Status get(const uint32_t i, const uint32_t j, T* val) const;
    Status set(const uint32_t i, const uint32_t j, T val);

    uint32_t nsamples() const { return nsamps_; };
    uint32_t cap() const { return capacity_; };

    bool is_null() const { return data_ == nullptr; };
private:

    // @brief Used by the init factor function
    UnnormalizedGrm<T>(uint32_t n, uint32_t caps, T* data_array)
        : nsamps_(n),
        capacity_(caps),
        data_(data_array) {};

    // @brief unsafe mapping of symmetric matrix to data index
    // @description mapping is unsafe on two accounts.  First, I assume
    //  that i <= j, if this cannot be assumed use midx_to_arr_ instead.
    //  Note, that the iltj in the function name is intended to
    //  communicate this assumption as the acronym iltj represents "i" 
    //  less than "j".
    //  Second, I do not verify that the reulting bounds index is within 
    //  the bounds of the data array.  I assume that the code calling this 
    //  has integrating this logic.  
    // @param[in] i matrix row index
    // @param[in] j matrix col index
    // @return index of data array of the requested data
    uint32_t sym_midx_to_arr_iltj_(const uint32_t i, const uint32_t j) const;

    // @brief unsafe mapping of symmetric matrix to data index
    // @description mapping is unsafe as I do not verify that the reulting 
    //  bounds index is within the bounds of the data array.  I assume that
    //  the code calling this has integrating this logic.
    // @param[in] i matrix row index
    // @param[in] j matrix col index
    // @return index of data array of the requ3ested data
    uint32_t sym_midx_to_arr_(const uint32_t i, const uint32_t j) const;

    uint32_t nsamps_;        
    uint32_t capacity_;      // size of allocated memory for data
    T* data_;
};

template <typename T>
using ugrm_t = std::unique_ptr<UnnormalizedGrm<T>>;


template <typename T>
ugrm_t<T> UnnormalizedGrm<T>::init(uint32_t n) {
    if (n == 0)
        return nullptr;

    uint32_t d = std::numeric_limits<uint32_t>::max() / n;
    
    // check for overflow issue
    if (d < (n + 1) / 2) {
        fprintf(stderr, 
                "Number of samples exceeds max, contact maintainer\n");
        return nullptr;
    }

    uint32_t capacity = n * (n + 1) / 2;

    T* data_array = new(std::nothrow) T[capacity];
    if (data_array == nullptr) 
        return nullptr;

    std::memset(data_array, static_cast<T>(0), sizeof(T) * capacity);

    UnnormalizedGrm<T>* ugrm = new(std::nothrow) UnnormalizedGrm<T>(n,
            capacity, data_array);

    if (ugrm == nullptr) {
        delete[] data_array;
        return nullptr;
    }

    return std::unique_ptr<UnnormalizedGrm<T>>(ugrm);
}


template <typename T>
UnnormalizedGrm<T>::~UnnormalizedGrm() {
    if (data_)
        delete[] data_;
    data_ = nullptr;
    nsamps_ = 0;
    capacity_ = 0;
}

template <typename T>
T UnnormalizedGrm<T>::operator()(const uint32_t i, const uint32_t j) const {
    return data_[sym_midx_to_arr_(i, j)];
}

template <typename T>
Status UnnormalizedGrm<T>::get(const uint32_t i, 
        const uint32_t j,
        T* val) const {
    if (i >= nsamps_ || j >= nsamps_)
        return Status::ErrIndexOutofBounds;

    *val = data_[sym_midx_to_arr_(i, j)];

    return Status::Success;
}


template <typename T>
Status UnnormalizedGrm<T>::set(const uint32_t i, 
        const uint32_t j,
        T val) {
    if (i >= nsamps_ || j >= nsamps_)
        return Status::ErrIndexOutofBounds;

    data_[sym_midx_to_arr_(i, j)] = val;
    return Status::Success;
}


template <typename T>
uint32_t UnnormalizedGrm<T>::sym_midx_to_arr_iltj_(const uint32_t i, 
        const uint32_t j) const {
    return i*nsamps_ - i*(i-1)/2 + j - i; 
};


template <typename T>
uint32_t UnnormalizedGrm<T>::sym_midx_to_arr_(const uint32_t i, 
        const uint32_t j) const {
    // remember that by symmetry, the matrix is equal to its transpose
    if (i <= j)
        return sym_midx_to_arr_iltj_(i, j);

    return sym_midx_to_arr_iltj_(j, i);
}

// int hap_update_kernel(Grm* grmat, const bcfio::BcfRecord<float>* rec);
template <typename T>
Status UnnormalizedGrm<T>::update(const bcfio::BcfRecord<T>* rec) {

    // instantiate indexing variables used in for loops
    // use static to prevent construction and destruction of variables
    // between function calls
    if (nsamps_ != rec->nrow)
        return Status::ErrDimensionsNotEqual;

    uint16_t k_cols = rec->ncol;
    uint16_t k_col = 0;

    const T* samp_i = rec->data;
    const T* samp_j = nullptr;

    T val = static_cast<T>(0);

    // only iterate over upper triangle
    // remember that record data is an n_sample by k haplotype matrix 
    for (uint32_t i = 0; i < nsamps_; i++) {

        samp_j = samp_i;
        for (uint32_t j = i; j < nsamps_; j++) {

            val = 0;

            for (k_col = 0; k_col < k_cols; k_col++)
                val += samp_i[k_col] * samp_j[k_col];

            // note that sym_matrix_idx_to_array_ assumes that index
            // i <= j this is implemented by the nested for loops.
            data_[sym_midx_to_arr_iltj_(i, j)] += val;
            samp_j += k_cols;
        }


        samp_i += k_cols;
    }

    return Status::Success;
}

}

#endif
