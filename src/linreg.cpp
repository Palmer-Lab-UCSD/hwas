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

// [[Rcpp::depends(RcppEigen)]]

#include "linreg.h"
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

#include "linreg_eigen.h"

// Calculate vector of residual sum of squares (RSS) from linear regression of Y on X
// [[Rcpp::export]]
NumericVector calc_rss_linreg(const NumericMatrix& X, const NumericMatrix& Y,
                              const double tol=1e-12)
{
    return calc_mvrss_eigenqr(X, Y, tol);
}

// Calculate just the coefficients from linear regression of y on X
// [[Rcpp::export]]
NumericVector calc_coef_linreg(const NumericMatrix& X, const NumericVector& y,
                               const double tol=1e-12)
{
    return calc_coef_linreg_eigenqr(X, y, tol);
}

// Calculate coefficients and SEs from linear regression of y on X
// [[Rcpp::export]]
List calc_coefSE_linreg(const NumericMatrix& X, const NumericVector& y,
                        const double tol=1e-12)
{
    return calc_coefSE_linreg_eigenqr(X, y, tol);
}

// Calculate matrix of residuals from linear regression of Y on X
// [[Rcpp::export]]
NumericMatrix calc_resid_linreg(const NumericMatrix& X, const NumericMatrix& Y,
                                const double tol=1e-12)
{
    return calc_resid_eigenqr(X, Y, tol);
}

// use calc_resid_linreg for a 3-dim array
// [[Rcpp::export]]
NumericVector calc_resid_linreg_3d(const NumericMatrix& X, const NumericVector& P,
                                   const double tol=1e-12)
{
    const int nrowx = X.rows();
    if(Rf_isNull(P.attr("dim")))
        throw std::invalid_argument("P should be a 3d array but has no dim attribute");
    const Dimension d = P.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("P should be a 3d array");
    if(d[0] != nrowx)
        throw std::range_error("nrow(X) != nrow(P)");

    NumericMatrix pr(nrowx, d[1]*d[2]);
    std::copy(P.begin(), P.end(), pr.begin()); // FIXME I shouldn't need to copy

    NumericMatrix result = calc_resid_eigenqr(X, pr, tol);
    result.attr("dim") = d;

    return result;
}

// least squares, returning everything
// output is list of (coef, fitted, resid, rss, sigma, rank, df, SE)
//
// argument se indicates whether to calculate standard errors (SE)
//
// [[Rcpp::export]]
List fit_linreg(const NumericMatrix& X, const NumericVector& y,
                const bool se=true, const double tol=1e-12)
{
    return fit_linreg_eigenqr(X, y, se, tol);
}
