
library(hwas)

source(system.file("utils/unittest.R", package = "hwas"))


set_up <- function() {
    bcf_filename <- system.file("exdata", "geno_test_data.bcf", 
                                package = "hwas",
                                mustWork = TRUE)


    samples <- readLines(system.file("exdata", "samples",
                                     package = "hwas",
                                     mustWork = TRUE))

    loci <- read.table(system.file("exdata", "pos.tsv",
                                   package = "hwas",
                                   mustWork = TRUE),
                        row.names = NULL,
                        header = TRUE)

    tmp <- 42
    return(structure(tmp,
                     samples = samples,
                     chrom = loci$chrom,
                     pos = loci$pos,
                     bcf = bcf_filename,
                     class = "data"))
}


test_bopen <- function() {
    test_data <- set_up()
    bid <- hwas::bopen(attr(test_data, "bcf"), "r")

    TEST_NOTNULL(bid);
    
    status <- hwas::bclose(bid)
    
    TEST_EQ(0, status)
}


if (!interactive()) {
    test_bopen()
    
    ERROR_REPORT()
}

