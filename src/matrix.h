// Matrix utilities
//
// Copyright (C) 2020 Karl Browman 
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// ###################################################################
// DISCLAIMER ON COPYRIGHT
// ###################################################################
//
// The above copyright notice was added by Robert Vogel to clearly
// distinguish the author's work, and license for which it is published, 
// from my own.  The source was copied from GitHub repository 
//
// https://github.com/rqtl/qtl2
//
// The DESCRIPTION file states that the author and current (April 2026)
// maintainer is Karl Broman, as such I have attributed the copyright
// to him despite not definitively seeing a statement of copyright.
// Moreover, I specified the year 2020 as, according to GitHub, the file
// was last modified 6 years ago.
//
// ###################################################################
// END DISCLAIMER
// ###################################################################
//
#ifndef MATRIX_H
#define MATRIX_H

#include <RcppEigen.h>

// find columns that exactly match previous columns
// returns numeric vector with -1 indicating no match to an earlier column and
//                             >0 indicating matches that earlier column
//                                (indexes starting at 1)
Rcpp::NumericVector find_matching_cols(const Rcpp::NumericMatrix& mat,
                                       const double tol);

// find set of linearly independent columns in a matrix
// returns a vector of column indices (starting at 1)
Rcpp::IntegerVector find_lin_indep_cols(const Rcpp::NumericMatrix& mat,
                                        const double tol);

// form X matrix with intcovar
// has_intercept = true indicates that addcovar has an intercept
//                 and so probs reduced by one column
//               = false means probs has full set of columns
//
// This is maybe a bit confusing.
//
// In the has_intercept=true case, an intercept is included in the
// addcovar matrix, and probs is missing the first column.
// The matrix formed is [A P (P.I)] where A=addcovar, P=probs, I=intcovar
//
// In the has_intercept=false case, no intercept is included in the
// addcovar matrix, and the probs have all columns (so each row sums
// to 1). The matrix formed is [P A (P*.I)] where in the P*.I bit we
// drop the first column of the probs when getting interactions with
// intercovar.
Rcpp::NumericMatrix formX_intcovar(const Rcpp::NumericVector& probs,
                                   const Rcpp::NumericMatrix& addcovar,
                                   const Rcpp::NumericMatrix& intcovar,
                                   const int position,
                                   const bool has_intercept);

// multiply each column of a matrix by a set of weights
Rcpp::NumericMatrix weighted_matrix(const Rcpp::NumericMatrix& mat,
                                    const Rcpp::NumericVector& weights);

// multiply each element of a vector by the corresponding weight
Rcpp::NumericVector weighted_3darray(const Rcpp::NumericVector& array,
                                     const Rcpp::NumericVector& weights);

// expand genotype probabilities with intcovar
Rcpp::NumericVector expand_genoprobs_intcovar(const Rcpp::NumericVector& probs, // 3d array ind x prob x pos
                                              const Rcpp::NumericMatrix& intcovar);


// matrix multiplication
Rcpp::NumericMatrix matrix_x_matrix(const Rcpp::NumericMatrix& X,
                                    const Rcpp::NumericMatrix& Y);

// multiply matrix by vector
Rcpp::NumericVector matrix_x_vector(const Rcpp::NumericMatrix& X,
                                    const Rcpp::NumericVector& y);

// multiply matrix by array
Rcpp::NumericVector matrix_x_3darray(const Rcpp::NumericMatrix& X,
                                     Rcpp::NumericVector& A);

#endif // MATRIX_H
