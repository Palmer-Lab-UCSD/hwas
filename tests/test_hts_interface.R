
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


unittests$TEST("TestHts", "test_open", function() {
    test_data <- set_up()
    bid <- hwas::bopen(attr(test_data, "bcf"), "r")

    Expect$not_null(bid);
    
    status <- hwas::bclose(bid)
    cat(sprintf("\n\nBopen status: %d\n\n\n", status))
    
    Expect$eq(0, status)
})


unittests$TEST("TestHts",
               "test_is_bcf",
               function() {

    bcf_compressed <- system.file("exdata", "geno_test_data_compressed.bcf",
                                  package = "hwas",
                                  mustWork = TRUE)
    bcf_uncompressed <- system.file("exdata", "geno_test_data_uncompressed.bcf",
                                    package = "hwas",
                                    mustWork = TRUE)
    vcf_uncompressed <- system.file("exdata", "geno_test_data.vcf",
                                    package = "hwas",
                                    mustWork = TRUE)
    vcf_compressed <- system.file("exdata", "geno_test_data.vcf.gz",
                                  package = "hwas",
                                  mustWork = TRUE)
    

    Expect$true(hwas::is_bcf(bcf_compressed))
    Expect$true(hwas::is_bcf(bcf_uncompressed))
    Expect$false(hwas::is_bcf(vcf_uncompressed))
    Expect$false(hwas::is_bcf(vcf_compressed))
})


if (!interactive()) {
    print("Running as not interactive")
    unittests$main()
}

