library(hwas)

bid <- hwas::bopen("../data/chr12.bcf", "rb")
success <- hwas::subset_samples(bid, "../data/sample_subset")
if (success != 0)
    stop("Subset vcf by sample failed.")

if (hwas::set_threads(bid, 8) != 0)
    stop("multithreading failure")

start_time <- Sys.time()
grm <- hwas::calc_grm(bid, "HD")
print(paste("Time to compute GRM elapsed:", Sys.time() - start_time))
