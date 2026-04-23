/* -*- C++ -*- */
// Palmer Lab at UCSD
//
// This library provides the data structure of a genetic relationship
// matrix and functions for file I/O.
// 
// GRM BINARY FILE SPECIFICATION
//
// A computed GRM is stored in a custom binary format.  The extension
// ".grm" of these files is mandatory.  The file is divided into two 
// components, a header with meta-data necessary to reproduce the grm
// calculation and the the computed grm values, named the payload.  
//
// The .grm file header is defined by the struct Hdr, and contains,
// at a minimum, the following information:
//
//      * program_version: grm program version number
//      * data_type: alt_count, 
//          expected_alt_count, 
//          expected_haplotype_count, 
//          both expected_alt_count and 
//          expected_haplotype_count.
//      * coords: Genomic coordinates used in the grm calculation.
//      * samples: list of sample id's in order of the grm
//
// the coords and samples are defined by their own structs with field
// pointers to heap allocated memory addresses.  Reading and writing 
// such heap allocated structs make use of runtime polymorphism of
// function "read" and "write"
//
//
// ACKNOWLEDGMENT
//
// Code design and original version completed by Robert Vogel,
// reviewed by Claude Opus 4.6, the AI assistant from Anthropic.
// Some recommendations have been incorporated.
// 
#ifndef HEADER_GRM_H
#define HEADER_GRM_H

#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cctype>
#include <memory>
#include <utility>
#include <string>

#include <bcfio.hpp>
#include <logger.hpp>
#include <io.hpp>
#include <constants.hpp>
#include <utils.hpp>



// The algorithm for getting the array idx from matrix indexes is
// simply
//
// idx = i * n_samples - n_skipped_idxs + j
//
// where i is the matrix row index and j is the matrix column index.
// The interesting term is n_skipped_idxs, this is the number of 
// elements that referencing (i, j) skip when only storing upper
// triangle. For example, suppose that i = 3 and j = 2.  Here three
// complete rows of the matrix has been traversed, therefore the number 
// skipped is
//
// n_skipped_idx >= (i-1) * i / 2
//
// if j < i, then as we are only storing the upper triangle and
// by symmetry swap the values of i an j.  If j > i,  then the number
// skipped in row i is i, making the total skipped,
//
// n_skipped_idx = i * (i - 1) / 2 + i
//
// making the equation above read
//  
//  idx = i * n_samples - i * (i - 1)/2 - i + j
// 
// Then let's manually validate
//
// n = 3
// i    j   idx_true    idx
// 0    0   0           0     
// 0    1   1           1
// 0    2   2           2
// 1    0   --- transpose -> 0, 1 ---
// 1    1   3           1*3 - 1*(0)/2 - 1 + 1 = 3
// 1    2   4           1*3 - 1*(0)/2 - 1 + 2 = 4
// 2    2   5           2*3 - 2*(1)/2 - 2 + 2 = 6 - 1 = 5
//
#define MATRIX_IDX_TO_ARRAY(i, j, n)   ((i)*(n) - (i)*((i)-1)/2 + (j) - (i))


namespace grm {


// Recall that the magic number is simply GRM\0 in hex,
// 0x47 = G, 0x52 = R, and so on.  The magic number provides
// a simple means to determine file type when parsing.
constexpr uint32_t FILE_TYPE_SPEC = 0x47524D00;
constexpr utils::Version FILE_VERSION = { 0, 0, 0 };
constexpr char FILE_SUFFIX[] = ".grm";


enum STATUS { 
    SUCCESS, 
    FAILED, 
    UNKNOWN_FAILURE,
    ERROR_IDX_ARR_BOUNDS,
    ERROR_FOPEN,
    ERROR_EOF_NOT_REACHED,
    ERROR_ON_WRITE,
    ERROR_ON_READ,
    ERROR_FILE_NOT_OPEN,
    ERROR_NULLPTR_ARG,
    ERROR_INVALID_ARG,
    ERROR_NOT_A_GRM_FILE,
    ERROR_BCF_ATTR,
    ERROR_BCF_IDX
};


enum GrmType {
    EHC,        // Expected Haplotype Count
    EAC,        // Expected Alternative Allele Count
    BOTH,       // Both EHC AND EAC
    DS,         // Dosage, i.e. Called Alternative Allele Count
    UNSPECIFIED,
}; 

// @title: Store genomic coordinates used in GRM calculation
struct Coordinates {
    Coordinates(): contig(""), len(0), pos(nullptr) {};
    Coordinates(char* contig_in, uint64_t len_in)
        : contig(contig_in == nullptr ? "" : contig_in),
        len(contig == "" ? 0 : len_in), 
        pos(len == 0 ? nullptr : std::make_unique<uint64_t[]>(len)) {};

    Coordinates(const Coordinates&) = delete;
    Coordinates& operator=(const Coordinates&) = delete;

    Coordinates(Coordinates&& other);
    Coordinates& operator=(Coordinates&& other);

    // Data Fields
    std::string contig;
    uint64_t len;
    std::unique_ptr<uint64_t[]> pos;
};


// Coordinates Storage Layout
//
//  type    number  description
//  -----------------------------------------------------------------
//  size_t  1       number of characters (n) in contig name
//  char    n       characters for contig name without null character
//  size_t  1       number of genomic positions (npos) 
//  size_t  npos    the positions used for computation of the grm
//
STATUS write(io::FileIO* fio, const Coordinates* coords);
STATUS read(io::FileIO* fio, Coordinates* coords);


// Samples stores sample id strings and the number of samples
//
struct Samples {
    Samples(): len(0), names(nullptr) {};
    Samples(uint64_t n_samples): 
        len(n_samples), 
        names(len == 0 ? nullptr : std::make_unique<std::string[]>(len)) 
        {};

    Samples(const Samples&) = delete;
    Samples& operator=(const Samples&) = delete;

    Samples(Samples&& other);
    Samples& operator=(Samples&& other);

    // Data Fields
    uint64_t len;               //number of samples
    std::unique_ptr<std::string[]> names;

};


// Sample Storage Layout
//
//  type    number  description
//  --------------------------------------------------------------------
//  size_t  1       represents number of samples
//  size_t  1       the number of characters of longest sample id
//  size_t  1       number of characters (n_1) in first sample id
//  char    n_1     characters of sample id 1 without terminal null '\0'
//  size_t  1       number of characters (n_2) in second sample id
//  char    n_2     characters of sample id 2 without terminal null '\0'
//  ...
//  size_t  1       number of characters (n_N) in N^{th} sample id
//  char    n_N     characters of sample id N without terminal null '\0'

STATUS write(io::FileIO* fio, const Samples* samples);
STATUS read(io::FileIO* fio, Samples* samples);


// Header
struct Hdr {

    Hdr();
    Hdr(const Hdr&) = delete;
    Hdr& operator=(const Hdr&) = delete;

    Hdr(Hdr&&);
    Hdr& operator=(Hdr&&);

    // Data Fields
    utils::Version prog_version;
    utils::Version file_version;

    GrmType grm_type;
    std::unique_ptr<Coordinates> coords;
    std::unique_ptr<Samples> samples;

};

// Header Storage Layout
//
//  type    number  description
//  --------------------------------------------------------------------
//  size_t  1       number of characters (n) in version string
//  char    n       version string without null terminator
//  GrmType 1       type of grm 
//
//  call write for coordinates
//
//  call write for samples
//
STATUS write(io::FileIO *fio, const Hdr *hdr);
STATUS read(io::FileIO *fio, Hdr *hdr);


// Grm class manages storage and access of GRM matrix
// 
// The GRM as an n_sample by n_sample symmetric, positive semi-definite
// matrix.  Let Z represent the n_sample by m_marker data genetic data.  
// From these data the GRM is computed as GRM = ZZ^T.
//
// @param n_samples of the GRM.
//
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
    STATUS set(const uint64_t i, const uint64_t j, const float val); 
    STATUS get(const uint64_t i, const uint64_t j, float *val) const; 

    uint64_t size() const;

    STATUS midx_to_arr(const uint64_t i, 
            const uint64_t j, uint64_t* idx) const;

    uint64_t n_samples;
    std::unique_ptr<float[]> data;
};

// @title: Write meta-data and computed grm elements to file
// @description: The binary file written contains a header and payload:
//      * Header
//          - an instance of grm::Header
//      * Payload
//          - grm data in row major order
// @param filename: name of file that the data are written
// @param hdr: an instance of grm::Hdr with important meta data
// @return grm::STATUS: 
//
STATUS write(io::FileIO* fio, const Hdr* hdr, const Grm* grmatrix);
STATUS read(io::FileIO* fio, Hdr* hdr, Grm* grmatrix);


//
STATUS calc_grm_ehc(Logger* log, bcfio::ReadBcf* bfid, Grm* grmatrix);

//int compute_eac_and_ehc_matrix();
// int compute_genotype_matrix();
// int compute_eac_matrix();

}

#endif
