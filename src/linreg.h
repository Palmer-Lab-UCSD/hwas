// generic linear regression functions
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
#ifndef LINREG_H
#define LINREG_H

#include <RcppEigen.h>

// Calculate vector of residual sum of squares (RSS) from linear regression of Y vs X
Rcpp::NumericVector calc_rss_linreg(const Rcpp::NumericMatrix& X,
                                    const Rcpp::NumericMatrix& Y,
                                    const double tol);

// Calculate just the coefficients from linear regression of y on X
Rcpp::NumericVector calc_coef_linreg(const Rcpp::NumericMatrix& X,
                                     const Rcpp::NumericVector& y,
                                     const double tol);

// Calculate the coefficients and SEs from linear regression of y on X
Rcpp::List calc_coefSE_linreg(const Rcpp::NumericMatrix& X,
                              const Rcpp::NumericVector& y,
                              const double tol);

// Calculate matrix of residuals from linear regression of Y on X
Rcpp::NumericMatrix calc_resid_linreg(const Rcpp::NumericMatrix& X,
                                      const Rcpp::NumericMatrix& Y,
                                      const double tol);

// use calc_resid_linreg for a 3-dim array
Rcpp::NumericVector calc_resid_linreg_3d(const Rcpp::NumericMatrix& X,
                                         const Rcpp::NumericVector& P,
                                         const double tol);


// least squares, returning everything
// output is list of (coef, fitted, resid, rss, sigma, rank, df, SE)
//
// argument se indicates whether to calculate standard errors (SE)
Rcpp::List fit_linreg(const Rcpp::NumericMatrix& X,
                      const Rcpp::NumericVector& y,
                      const bool se,
                      const double tol);

#endif // LINREG_H
