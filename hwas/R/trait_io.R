
BUFFER_SIZE <- 100
META_PREFIX <- "##"
HEADER_PREFIX <- "#"
CSV_DELIMITER <- ","

read_data <- function(filename, delimiter = CSV_DELIMITER) {

    if (!file.exists(filename))
        stop(sprintf("File %s not found", filename))


    fid <- file(filename, "r")
    meta <- vector(mode="character", length=BUFFER_SIZE)

    i <- 1
    tline <- readLines(fid, n = 1)
    while (startsWith(tline, META_PREFIX)) {
        meta[i] <- tline
        tline <- readLines(fid, n = 1)
        i <- i + 1
    }

    meta <- meta[1:i-1]

    if (startsWith(tline, HEADER_PREFIX))
        header <- strsplit(sub(HEADER_PREFIX, "", tline), delimiter)[[1]]

    data <- read.table(fid, header=FALSE, col.names=header, sep=delimiter)



    if (isOpen(fid))
        close(fid)

    return(list("meta" = meta, "data" = data))
}


write_data <- function(filename, data_set) {

    fid <- file(filename, 'w')
    cat(data_set$meta, file = fid, sep = '\n')

    cnames <- colnames(data_set$data)
    cnames[1] <- paste0(HEADER_PREFIX, cnames[1])
    colnames(data_set$data) <- cnames
    write.table(data_set$data,
                file = fid,
                row.names = FALSE,
                quote = FALSE,
                append = TRUE,
                sep = CSV_DELIMITER)

    if (isOpen(fid))
        close(fid)
}
