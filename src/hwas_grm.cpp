
#include <hwas.h>


// For signal handling see
// https://en.cppreference.com/cpp/utility/program/signal
// I need to review this more closely.  The use of atomic
// indicates that this is using mutex free concurrency.
volatile std::sig_atomic_t signal_received = 0;

void signal_handler(int signal) {
    signal_received = 1;
}


// [[Rcpp::export]]
Rcpp::RObject calc_grm(bcf_conn_t bid, const char* id) {

    // reset signal for subsequent usage.
    signal_received = 0;

    std::signal(SIGINT, signal_handler);

    // instantiate matrices to hold calculations
    uint16_t k = 0;
    if (bcfio::k_fmt(bid, id, &k) != 0)
        return R_NilValue;

    uint32_t nsamps = 0;
    if (num_samples(bid, &nsamps) != 0)
        return R_NilValue;

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
        s = bcfio::next_record(bid.checked_get(), &rec, id);
        end_interval = clock::now();
        t_io += end_interval - start_interval;

        if (s != 0) break; 

        start_interval = clock::now();
        if (grm::hap_update_kernel(&grmat, &rec) != 0) 
            return R_NilValue;

        end_interval = clock::now();
        t_compute += end_interval - start_interval;

        if (++rec_count % 1000 == 0)
            Rprintf("Processed %zu records\n", rec_count);
    }

    if (signal_received)
        return R_NilValue;

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
