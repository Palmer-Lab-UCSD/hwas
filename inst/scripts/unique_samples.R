
main <- function(geno_filename, exclude_filename) {
    gsamples <- readLines(geno_filename)
    exsamples <- readLines(exclude_filename)

    samples <- rep(NA, length(gsamples))

    j <- 1
    for (samp in gsamples) {
        if (samp %in% exsamples)
            next

        samples[j] <- samp
        j <- j + 1
    }

    writeLines(
}


if (!interactive()) {
    args <- commandArgs(trailingOnly=TRUE)

    if (length(args) != 2)
        stop("insufficient input")

    for (arg in args)
        if (!file.exists(arg))
            stop(paste("File,", arg, "Does not exist"))

    main(args[1], args[2])
}
