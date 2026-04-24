// Palmer Lab at UCSD 2026
//
#include <hwas_types.h>

// @return -1 on error and 0 or success
static int update(Rcpp::NumericMatrix& grmatrix, bcfio::BcfFloatRecord* rec) {

    // instantiate indexing variables used in for loops
    // use static to prevent construction and destruction of variables
    // between function calls
    uint64_t grow = 0;
    uint64_t gcol = 0;
    uint64_t k_hap = 0;
    std::optional<float> val_i = std::nullopt;
    std::optional<float> val_j = std::nullopt;
    float val = 0;
    uint64_t n_samples = rec->nrows();
    uint64_t k_haps = rec->ncols();

    // only iterate over upper triangle
    // remember that record data is an n_sample by k haplotype matrix 
    for (grow = 0; grow < n_samples; grow++) {
        for (gcol = grow; gcol < n_samples; gcol++) {

            val = 0;

            for (k_hap = 0; k_hap < k_haps; k_hap++) {

                if ((val_i = rec->get(grow, k_hap)) == std::nullopt) {
                    Rprintf("ERROR: grow, k_hap: %lu, %lu\n", grow, k_hap);
                    return -1;
                }

                if ((val_j = rec->get(gcol, k_hap)) == std::nullopt) {
                    Rprintf("ERROR: gcol, k_hap: %lu, %lu\n", gcol, k_hap);
                    return -1;
                }

                val += val_i.value() * val_j.value();
            } 

            grmatrix(grow, gcol) += val;
        }
    }

    return 0;
}


// [[Rcpp::export]]
Rcpp::RObject calc_grm(Rcpp::XPtr<bcfio::Bcf> bid, const char* id) {

    // instantiate matrices to hold calculations
    int32_t k { 0 };
    if ((k = k_fmt(bid, id)) < 0) {
        return R_NilValue;
    }
    uint64_t nsamps = num_samples(bid);

    bcfio::BcfFloatRecord rec {};
    Rcpp::NumericMatrix grmatrix(nsamps, nsamps);

    size_t rec_count = 0;
    while (bcfio::next_record(bid.get(), &rec, id) == 0) {

        if (update(grmatrix, &rec) != 0) {
            return R_NilValue;
        }

        if (++rec_count % 1000 == 0)
            Rprintf("Processed %zu records\n", rec_count);
    }


    // fill in lower diagonal elements
    for (uint64_t i = 0; i < nsamps; i++)
        for (uint64_t j = 0; j < i; j++)
            grmatrix(i, j) = grmatrix(j, i);
            
    Rprintf("Done\n")
    return grmatrix;
}
