#ifndef RPGSIM_H
#define RPGSIM_H


#include <Rcpp.h>
#include <Rbcfio.h>


namespace rpgsim {
enum struct Status : int {
    Success                     = 0,
    ErrNotSymmetricMatrix       = -2,
    ErrHeritabilityOutOfRange   = -3
};

const char* status_msg(Status status);
}


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
    double var_u;   // variance of random polygenic effects
    double var_e;   // environment variance
};


//  Consequently, to simulate phenotype data we need to provide expected
//  founder haplotype counts, the variance of random effect sizes 
//  (\sigma_u), the variance of random effect sizes(\sigma_e), and fixed
//  effect sizes.  These parameters alone do not provide any intution on
//  how the simulated data will look.  It would be better provide the
//  heritability, and solve for these parameters in terms of the input
//  heritability.  This can be done be computing the heritability under
//  the model.  By definition of narrow sense heritability and a polygenic
//  random effects model
//
//  h^2 = var(G) / (var(G) + var(E))
//
//  G_i = \sum_k Z_{ik}^T U_k is simply the polygenic prediction of 
//  phenotype from genotypes.  The variance of G under the model can be
//  computed as the expected value of the empirical variance over the
//  population of N samples
//
//  E(var(G)) = E( \frac{1}{N-1} \sum_{i} (G_i - \bar{G})^2 )
//
//  which after some arithmetic and plugging in model values and applying
//  definitions
//
//  E(var(G)) = \frac{\sigma_u^2}{N-1} [
//                  Trace(A) - \frac{1}{N} S
//              ]
//  
//  where A = \sum_k Z_k Z_k^T is the GRM and S = \sum_i\sum_j A_{ij} is 
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
rpgsim::Status calc_pgsim_pars (const double heritability, 
        const Rcpp::NumericMatrix& grm,
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
        const Rcpp::NumericMatrix& grm,
        const float qtl_n,
        const float qtl_effect_mf,
        const float heritability,
        const char* id);

#endif
