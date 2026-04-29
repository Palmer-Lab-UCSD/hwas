
#include <hwas_types.h>

volatile std::sig_atomic_t signal_received = 0;

void signal_handler(int signal) {
    signal_received = 1;
}

// [[Rcpp::export]]
Rcpp::RObject calc_grm(Rcpp::XPtr<bcfio::Bcf> bid, const char* id) {

    std::signal(SIGINT, signal_handler);

    // instantiate matrices to hold calculations
    int32_t k { 0 };
    if ((k = k_fmt(bid, id)) < 0) {
        return R_NilValue;
    }
    uint64_t nsamps = num_samples(bid);

    grm::Grm grmat { nsamps };
    bcfio::BcfRecord<float> rec {};

    using clock = std::chrono::steady_clock;
    auto t_io = clock::duration::zero();
    auto t_compute = clock::duration::zero();

    auto start_interval = clock::now();
    auto end_interval= clock::now();
    int s = 0;

    size_t rec_count = 0;
    while (!signal_received) {
        start_interval = clock::now();
        s = bcfio::next_record(bid.get(), &rec, id);
        end_interval = clock::now();
        t_io += end_interval - start_interval;

        if (s != 0) break; 

        start_interval = clock::now();
        if(grm::hap_update_kernel(&grmat, &rec) != 0) return R_NilValue;
        end_interval = clock::now();
        t_compute += end_interval - start_interval;

        if (++rec_count % 1000 == 0)
            Rprintf("Processed %zu records\n", rec_count);
    }

    Rcpp::Rcout 
        << "I/O time: " 
        << std::chrono::duration_cast<std::chrono::seconds>(t_io).count()
        << std::endl;

    Rcpp::Rcout 
        << "Compute time: " 
        << std::chrono::duration_cast<std::chrono::seconds>(t_compute).count()
        << std::endl;

    Rcpp::NumericMatrix grmatrix(nsamps, nsamps);
    // fill in lower diagonal elements
    for (uint64_t i = 0; i < nsamps; i++)
        for (uint64_t j = 0; j < nsamps; j++)
            grmatrix(i, j) = grmat(j, i);
            
    Rprintf("Done\n");
    return grmatrix;
}
