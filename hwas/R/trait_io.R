
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
