#ifndef RMATUTILS_H
#define RMATUTILS_H

#include <Rcpp.h>

#include <status.h>

namespace matutils {

bool is_symm_(const Rcpp::NumericMatrix& A);

// @brief compute the sum over all matrix elements
// @param[in] A the matrix 
// @param sum of elements of A
Status matrix_sum(const Rcpp::NumericMatrix& A, double& val);


// @brief compute the trace of a square matrix
// @param[in] the square matrix for computing trace
// @return < 0 upon error, otherwise positive double equal to the 
//  matrix trace
Status matrix_trace(const Rcpp::NumericMatrix& A, double& val);

}

#endif
