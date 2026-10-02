


# @brief compute the loco GRM
# @param 
mk_loco_grm <- function(...) {
    gsims <- list(...)
    if (length(gsims) == 0)
        stop("Need to specify similarity matrices to combine")

    for (w in gsims) {
        if (attr(w, "class") != "gsim")
            stop("All inputs matrices must be of class gsim")


    }

    nsamps <- dim(gsims[[1]])[1]
    loco <- matrix(0, nsamps, nsamps)
    for (w in gsims) {
        loco <- loco + w
    }

    gsim_to_grm(loco)
    attr(loco, "class") <- "loco_grm"
    return(loco)
}
