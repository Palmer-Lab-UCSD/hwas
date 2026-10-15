


# @brief Accumulate genetic similarity matrices
# @param ... gsim matrices
# @return A gsim matrix instance that is the sum of the input 
#   gsim matrices
accumulate_gsims <- function(...) {
    gsims <- list(...)
    if (length(gsims) == 0)
        stop("Need to specify similarity matrices to combine")

    format_key <- attr(gsims[[1]], "format_key")

    samp_names <- rownames(gsims[[1]])
    nsamps <- length(samp_names)
    
    agsim <- matrix(0, nsamps, nsamps)

    rownames(agsim) <- samp_names
    colnames(agsim) <- samp_names

    m_pos <- 0
    for (w in gsims) {
        if (attr(w, "class") != "gsim")
            stop("All inputs matrices must be of class gsim")

        if (any(samp_names != rownames(w)) ||
            any(samp_names != colnames(w))) {
            stop("Sample names of matrices do not align")
        }

        if (attr(w, "format_key") != format_key)
            stop("Similarity matrices are of different data types")

        agsim <- agsim + w
        m_pos <- m_pos + attr(w, "m_pos")
    }

    attr(agsim, "class") <- "gsim"
    attr(agsim, "m_pos") <- m_pos
    attr(agsim, "format_key") <- format_key

    return(agsim)
}
