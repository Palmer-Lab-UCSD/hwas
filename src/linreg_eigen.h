// linear regression via RcppEigen
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
#ifndef LINREG_EIGEN_H
#define LINREG_EIGEN_H

#include <RcppEigen.h>

// calc X'X
Eigen::MatrixXd calc_XpX(const Eigen::MatrixXd& A);

// least squares by "LLt" Cholesky decomposition
// needs to be full rank
Rcpp::List fit_linreg_eigenchol(const Rcpp::NumericMatrix& X,
                                const Rcpp::NumericVector& y,
                                const bool se);

// least squares by "LLt" Cholesky decomposition
// return just the coefficients
Rcpp::NumericVector calc_coef_linreg_eigenchol(const Rcpp::NumericMatrix& X,
                                               const Rcpp::NumericVector& y);

// least squares by "LLt" Cholesky decomposition
// this returns the coefficients and SEs
Rcpp::List calc_coefSE_linreg_eigenchol(const Rcpp::NumericMatrix& X,
                                        const Rcpp::NumericVector& y,
                                        const double tol);

// least squares by "LLt" Cholesky decomposition
// return just the residual sum of squares
// needs to be full rank
double calc_rss_eigenchol(const Rcpp::NumericMatrix& X,
                          const Rcpp::NumericVector& y);

// least squares by "LLt" Cholesky decomposition
// return just the fitted values
Rcpp::NumericVector calc_fitted_linreg_eigenchol(const Rcpp::NumericMatrix& X,
                                                 const Rcpp::NumericVector& y);

// least squares by QR decomposition with column pivoting
Rcpp::List fit_linreg_eigenqr(const Rcpp::NumericMatrix& X,
                              const Rcpp::NumericVector& y,
                              const bool se,
                              const double tol);

// least squares by QR decomposition with column pivoting
// this just returns the coefficients
Rcpp::NumericVector calc_coef_linreg_eigenqr(const Rcpp::NumericMatrix& X,
                                             const Rcpp::NumericVector& y,
                                             const double tol);

// least squares by QR decomposition with column pivoting
// this returns the coefficients and SEs
Rcpp::List calc_coefSE_linreg_eigenqr(const Rcpp::NumericMatrix& X,
                                      const Rcpp::NumericVector& y,
                                      const double tol);

// least squares by QR decomposition with column pivoting
// return just the residual sum of squares
double calc_rss_eigenqr(const Rcpp::NumericMatrix& X,
                        const Rcpp::NumericVector& y,
                        const double tol);

// least squares by QR decomposition with column pivoting
// return just the fitted values
Rcpp::NumericVector calc_fitted_linreg_eigenqr(const Rcpp::NumericMatrix& X,
                                               const Rcpp::NumericVector& y,
                                               const double tol);

// least squares by "LLt" Cholesky decomposition, with matrix Y
// return vector of RSS
Rcpp::NumericVector calc_mvrss_eigenchol(const Rcpp::NumericMatrix& X,
                                         const Rcpp::NumericMatrix& Y);

// least squares by QR decomposition with column pivoting, with matrix Y
// return vector of RSS
Rcpp::NumericVector calc_mvrss_eigenqr(const Rcpp::NumericMatrix& X,
                                       const Rcpp::NumericMatrix& Y,
                                       const double tol);

// least squares by "LLt" Cholesky decomposition, with matrix Y
// return matrix of residuals
Rcpp::NumericMatrix calc_resid_eigenchol(const Rcpp::NumericMatrix& X,
                                         const Rcpp::NumericMatrix& Y);

// least squares by QR decomposition with column pivoting, with matrix Y
// return matrix of residuals
Rcpp::NumericMatrix calc_resid_eigenqr(const Rcpp::NumericMatrix& X,
                                       const Rcpp::NumericMatrix& Y,
                                       const double tol);

#endif // LINREG_EIGEN_H
