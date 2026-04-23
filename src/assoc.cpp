// Iterate
//


#include <Rcpp.h>
#include <string>
#include <io.h>
#include <bcfio.h>
#include <grm.h>
#include "fit1_pg.h"

#define HAPCONST "HD"

// [[Rcpp::export]]
Rcpp::List assoc(const std::string bcfname,
        const char* position_file,
        const Rcpp::NumericVector& phenotypes,
        const Rcpp::NumericMatrix& covariates, 
        const std::string grm_filename,
        const std::string pos_file_name) {
    // Load GRM and eigen decomposition
    
    grm::Grm grmatrix {};
    io.FileIO fio = 
    grm::read(

    // initialize the genetics data
    bcfio::ReadBcf bid = bcfio::open(bcfname, "rb");
    bid.set_samples(pos_file_name.c_str());

    bcfio::BcfFloatRecord rec {};

    Rcpp::NumericMatrix haplotypes = Rcpp::NumericMatrix::create();

    while (bid.next_record(&rec, HAPCONST)) {
        rec_to_rcpp(&rec, &haplotypes);

        rqtl::fit1_pg_addcovar(haplotypes, phenotypes,
                 )   
    }
}
