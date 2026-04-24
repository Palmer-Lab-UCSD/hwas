// fit a single-QTL model at a single position by LMM
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
#ifndef FIT1_PG_H
#define FIT1_PG_H

#include <Rcpp.h>

// fit single-QTL model at a single position
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// eigenvec  = eigenvectors from eigen decomposition of kinship matrix
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list with a bunch of stuff
//
Rcpp::List fit1_pg_addcovar(const Rcpp::NumericMatrix& genoprobs,
                            const Rcpp::NumericVector& pheno,
                            const Rcpp::NumericMatrix& addcovar,
                            const Rcpp::NumericMatrix& eigenvec,
                            const Rcpp::NumericVector& weights,
                            const bool se,
                            const double tol);


// fit single-QTL model at a single position
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// intcovar  = interactive covariates (should also be included in addcovar)
// eigenvec  = eigenvectors from eigen decomposition of kinship matrix
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list of a bunch of stuff
//
Rcpp::List fit1_pg_intcovar(const Rcpp::NumericMatrix& genoprobs,
                            const Rcpp::NumericVector& pheno,
                            const Rcpp::NumericMatrix& addcovar,
                            const Rcpp::NumericMatrix& intcovar,
                            const Rcpp::NumericMatrix& eigenvec,
                            const Rcpp::NumericVector& weights,
                            const bool se,
                            const double tol);

#endif // FIT1_PG_H
