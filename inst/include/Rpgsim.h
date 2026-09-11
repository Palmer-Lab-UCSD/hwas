#ifndef RPGSIM_H
#define RPGSIM_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Rbcfio.h>


namespace rpgsim {
enum struct Status : int {
    Success                     = 0,
    ErrNotSymmetricMatrix       = -2,
    ErrNotSquareMatrix          = -3,
    ErrHeritabilityOutOfRange   = -4
};

const char* status_msg(Status status);
}

bool is_symm_(const Rcpp::NumericMatrix& A);

// @brief compute the sum over all matrix elements
// @param[in] A the matrix 
// @param sum of elements of A
rpgsim::Status matrix_sum(const Rcpp::NumericMatrix& A, double& val);


// @brief compute the trace of a square matrix
// @param[in] the square matrix for computing trace
// @return < 0 upon error, otherwise positive double equal to the 
//  matrix trace
rpgsim::Status matrix_trace(const Rcpp::NumericMatrix& A, double& val);


struct PgSimParams {
    double qtl_freq;
    double qtl_effect_size;
    double var_u;   // variance of random polygenic effects
    double var_e;   // environment variance
    double heritability;
};




// @brief Computer variance component parameters
// @description To simulate phenotype data we need to know the the 
//  the fixed effects and random effects.  This function computers the
//  parameters relevant to the random effects.  The reason it is
//  necessary is that the random component parameters, i.e. the variance 
//  of random effects per marker (var_u) and the environment variance 
//  (var_e), are not easily interpretable.  Instead, the heritability is 
//  the more interpretable quantity, and often is our quantity of interest.  
//  This function computes var_u from the GRM, heritability, and under
//  the assignment var_e = 1.  The derivation is as follows.
//
//  Recall the definition of narrow sense heritability 
//
//  h^2 := var(G) / (var(G) + var(E))
//
//  Yang et al. 2010 established that additive polgenic effects can
//  effectively modeled as random effects in a linear mixed model.
//  Similarly, Broman et al. Genetics 2019 provides a random effects
//  model for mapping population founder haplotypes to phennotypes.  
//  In similar fashion we apply a random effects model to describe
//  the polygenic effects of expected count of founder haplotypes,
//  Z_{ijk}, of sample i, locus j, and haplotype k by defining
//
//  G_i = \sum_k \sum_j Z_{ijk}^T U_{jk},
//  U_{jk}  ~  Normal(0, var_u).
//  E_i ~ Normal(0, var_e)
//  
//  Under the definition of narrow sense heritability and the polygenic
//  model above we are able solve var_u.  To do this we need to compute
//  variances under model assumptions.
//
//  The variance of G under the model is computed as the expected
//  value of the empirical variance over the population of N samples
//
//  E(var(G)) = E( \frac{1}{N-1} \sum_{i} (G_i - \bar{G})^2 )
//
//  which after some arithmetic and plugging in model values and applying
//  definitions
//
//  E(var(G)) = \frac{\sigma_u^2}{N-1} [
//                  Trace(A) - \frac{1}{N} S
//              ].
//  
//  where, 
//
//  A := \sum_k Z_k Z_k^T 
//
//  A is our defined N x N GRM and the matrix notation Z_k is an N sample
//  by M marker matrix of expected counts of founder haplotype k.
//  The quantity 
//
//  S = \sum_i\sum_j A_{ij} is 
//
//  the sum over every element of the GRM.  The variance of the
//  phenotype attributed to the environment is computed simiarly.  Due
//  to the assumption of independence,
//
//  E(var(E)) = \sigma_e^2
//
//  Assessing our goal.  We have two unknown parameters, \sigma_u^2 and
//  \sigma_e^2 and one equation.  We cannot solve for both parameters.
//  However, as we do not care about the exact value of the parameters 
//  themeselves, but instead only the heritability, we can set 
//
//  \sigma_e^2 = 1
//
//  and solve for \sigma_u^2 in terms of the heritability and GRM.  This
//  procedure helps us tune the parameters of the polygenic effects, 
//  however exclude the QTL effects.  
// @param[in] grm N x N sample genetic relationship matrix defined as 
//  matrix A above
// @param[in,out] params are the simulation parameters updated by the
//  this function
// @return rpgsim namespace status code
rpgsim::Status calc_pgsim_pars_(const Rcpp::NumericMatrix& grm,
        PgSimParams& params);


// @brief Simulate sample phenotypes
// @description The simulation of phenotype y_i from sample i is a sample
//  from the linear model
//
//  Y_i = \sum_k \sum_j x_{ijk} \beta_{jk} + Z_{ik}^T U_k + e_i
//
//  with random variables
//
//  U_k     \sim Normal(0, \sigma_u^2 I)
//  e_i     \sim Normal(0, \sigma_e^2 I)
//
//  and fixed variables
//
//  x_{ijk}      is the expected count of founder haplotype k, for 
//               sample i, and genomic position j.
//  \beta_{jk}   effect size of founder haplotype k on genomic position j
//  Z_{ik}^T     1 x M matrix of expected founder haplotype counts of 
//               founder k, in sample i, at genomic positions 1, 2, ..., M.
//  
//
//  I do not know how to rigorously include the QTL effects, so I have
//  an ad hoc procedure.  The consequence is that the input heritability
//  should be a lower bound to the heritability of the simulation 
//  results as the results will have the contributions of both the
//  polygenic effects and fixed QTL effects.  
//
// @param[in] bconn open connection to bcf file of expected haplotype
//  counts.
// @param[in] qtl_n is the number of qtls to simulate.
// @param[in] qtl_effect_size the effect size of qtl's.  All have the
//  same size.
// @param[in] heritability is the narrow sense heritability on the
//  interval [0,1] 
// @return N length vector of phenotype values per sample
Rcpp::NumericVector pg_sim_qtl(const bconn_t bconn,
        const Rcpp::NumericMatrix& grmatrix,
        const float qtl_n,
        const float qtl_effect_mf,
        const float heritability,
        const char* id);


Rcpp::NumericVector pg_sim(const Rcpp::NumericMatrix& grmatrix,
        const float heritability);


// template <typename T>
// Rcppp::NumericVector pg_sim_qtl_(bcfio::Bcf* bid, 
//         const Rcpp::NumericMatrix& grmatrix,
//         const PgSimParams& params,
//         const char* id) {
//     // convert grmatrix to UnnormalizedGrm
// }

#endif
