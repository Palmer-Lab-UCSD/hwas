#
# Copyright (C) 2020 Karl Browman 
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# This code was developed following study of Nick Furlotte's pylmm code
# (https://github.com/nickFurlotte/pylmm).
#
# ###################################################################
# DISCLAIMER ON COPYRIGHT
# ###################################################################
#
# The above copyright notice was added by Robert Vogel to clearly
# distinguish the author's work, and license for which it is published, 
# from my own.  The source was copied from GitHub repository 
#
# https://github.com/rqtl/qtl2
#
# The DESCRIPTION file states that the author and current (April 2026)
# maintainer is Karl Broman, as such I have attributed the copyright
# to him despite not definitively seeing a statement of copyright.
# Moreover, I specified the year 2020 as, according to GitHub, the file
# was last modified 6 years ago.
#
#####################################################################
# END DISCLAIMER
#####################################################################

#' Batch columns by pattern of missing values
#'
#' Identify batches of columns of a matrix that have the same pattern
#' of missing values.
#'
#' @param mat A numeric matrix
#' @param max_batch Maximum batch size
#'
#' @return A list containing the batches, each with two components:
#' `cols` containing numeric indices of the columns in the
#' corresponding batch, and `omit` containing a vector of row indices
#' that have missing values in this batch.
#'
#' @export
#' @seealso [batch_vec()]
#' @examples
#' x <- rbind(c( 1,  2,  3, 13, 16),
#'            c( 4,  5,  6, 14, 17),
#'            c( 7, NA,  8, NA, 18),
#'            c(NA, NA, NA, NA, 19),
#'            c(10, 11, 12, 15, 20))
#' batch_cols(x)
batch_cols <-
    function(mat, max_batch=NULL)
{
    if(is.null(mat)) stop("mat is NULL")
    if(!is.matrix(mat)) mat <- as.matrix(mat)

    if(is.null(max_batch) || max_batch<=0)
        max_batch <- ncol(mat)

    mat <- !is.finite(mat)
    n <- nrow(mat)
    all_true <- rep(TRUE, n)
    result <- NULL

    # first pull out the columns that have 0 NAs
    n_na <- colSums(mat)
    no_na <- (n_na==0)
    if(any(no_na))
        result <- list(list(cols=which(no_na),
                                omit=numeric(0)))

    # now deal with the columns that have exactly 1 NA
    one_na <- (n_na==1)
    if(any(one_na)) {
        wh <- apply(mat[,one_na,drop=FALSE], 2, which)
        spl <- split(which(one_na), wh)
        part2 <- lapply(seq_along(spl), function(i)
                        list(cols=as.numeric(spl[[i]]),
                             omit=as.numeric(names(spl)[i])))

        if(is.null(result)) result <- part2
        else result <- c(result, part2)
    }

    # now the columns with >1 NA
    other_cols <- !(no_na | one_na)
    if(any(other_cols)) {
        other_cols <- seq_len(ncol(mat))[other_cols]

        # pattern of missing data (as character string with 0's and 1's
        pat <- apply(mat[,other_cols,drop=FALSE], 2, function(a) paste(which(a), collapse=":"))

        # the unique patterns
        u <- unique(pat)

        # the result
        part3 <- lapply(u, function(a)
                        list(cols=other_cols[pat==a],
                             omit=as.numeric(strsplit(a, ":")[[1]])))

    } else part3 <- NULL

    if(is.null(result)) result <- part3
    else result <- c(result, part3)

    # reduce to max_batch
    sizes <- vapply(result, function(a) length(a$cols), 1)
    big <- sizes > max_batch
    if(!any(big)) return(result)

    toreduce <- which(big)
    for(i in toreduce) {
        omit <- result[[i]]$omit
        col_list <- batch_vec(result[[i]]$cols, max_batch)
        this <- lapply(col_list, function(a) list(cols=a, omit=omit))
        result <- c(result, this)
    }

   result[-toreduce]
}
