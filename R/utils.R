# utilities from several qtl2 files
#
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

############################################################
# weight_util 
############################################################
# are the weights NULL or all 0's?
is_null_weights <-
    function(weights, tol=1e-12)
{
    if(is.null(weights) || max(abs(weights-1)) < tol)
        return(TRUE)

    FALSE
}


# check vector of weights & take square-roots
#
# - all positive?
#
# - if no missing values and all close to 1,
#   just use NULL rather than the weights

sqrt_weights <-
    function(weights, tol=1e-12)
{
    if(is.null(weights)) return(weights)

    if(any(weights <= 0))
        stop("weights must all be positive")

    if(all(!is.na(weights) & abs(weights - 1)<tol))
        return(NULL)

    weights <- stats::setNames( as.numeric(weights), names(weights) )

    return(sqrt(weights))
}


# multiply a vector by a set of weights
#' @importFrom stats setNames
weight_vector <-
    function(vec, weights, tol=1e-12)
{
    if(is_null_weights(weights, tol) || is.null(vec)) return(vec)

    # align and multiply
    id <- get_common_ids(setNames(names(vec), NULL), setNames(names(weights), NULL))
    vec[id] * weights[id]
}


# multiply a matrix by a set of weights
#' @importFrom stats setNames
weight_matrix <-
    function(mat, weights, tol=1e-12)
{
    if(is_null_weights(weights, tol) || is.null(mat)) return(mat)

    # force the input mat to be a matrix
    if(!is.matrix(mat)) mat <- as.matrix(mat)

    # align and multiply
    id <- get_common_ids(setNames(rownames(mat), NULL), setNames(names(weights), NULL))
    mat[id,,drop=FALSE] * weights[id]
}

# multiply an array by a set of weights
#' @importFrom stats setNames
weight_array <-
    function(arr, weights, tol=1e-12)
{
    if(is_null_weights(weights, tol) || is.null(arr)) return(arr)

    # force the input mat to be a matrix
    if(!is.array(arr) || length(dim(arr)) != 3)
        stop("arr should be a 3-dimensional array")

    # align and multiply
    id <- get_common_ids(setNames(rownames(arr), NULL), setNames(names(weights), NULL))
    arr[id,,,drop=FALSE] * weights[id]
}


############################################################
# arg_util
############################################################

#
# dotargs = list(...) from function call
# argname = character string of the argument to grab
# default = default value for argument
# values  = optional vector of character strings with possible values
grab_dots <-
    function(dotargs, argname, default, values=NULL)
{
    if(argname %in% names(dotargs)) {
        arg <- dotargs[[argname]]
        if(!is.null(values) && !(arg %in% values)) {
            warning(argname, ' "', arg, '" not valid; using "',
                    default, '".')
            arg <- default
        }
    }
    else arg <- default
    arg
}

# give warning if there are extra named arguments
#
# dotargs = list(...) from function call
# expected = vector of character strings of the names of possible args anticipated
check_extra_dots <-
    function(dotargs, expected)
{
    args <- names(dotargs)
    if(!all(args %in% expected)) {
        extra <- args[!(args %in% expected)]
        if(length(extra) == 1)
            warning("Extra argument ignored: ", extra)
        else
            warning("Extra arguments ignored: ",
                    paste(extra, collapse=" "))
    }
}

############################################################
# cluster_util
############################################################

# code related to clusters

# test if input is a prepared cluster (vs. just a number)
is_cluster <-
    function(cores)
{
    inherits(cores, "cluster") && inherits(cores, "SOCKcluster")
}

# number of cores being used
n_cores <-
    function(cores)
{
    if(is_cluster(cores)) return( length(cores) )
    cores
}

# set up a cluster
setup_cluster <-
    function(cores, quiet=TRUE)
{
    if(is_cluster(cores)) return(cores)

    if(is.null(cores) || is.na(cores)) cores <- 1
    if(cores==0) cores <- parallel::detectCores() # if 0, detect cores
    if(is.na(cores)) cores <- 1

    if(cores > 1 && Sys.info()[1] == "Windows") { # windows doesn't support mclapply
        cores <- parallel::makeCluster(cores)
        # the following calls on.exit() in the function that called this one
        # see http://stackoverflow.com/a/20998531
        do.call("on.exit",
                list(quote(parallel::stopCluster(cores))),
                envir=parent.frame())
    }
    cores
}

# run code by cluster (generalizes lapply, parLapply, and mclapply)
# (to deal with different methods on different architectures)
# if cores==1, just use lapply
cluster_lapply <-
    function(cores, ...)
{
    if(is_cluster(cores)) { # cluster object; use mclapply
        return( parallel::parLapply(cores, ...) )
    } else {
        if(cores==1) return( lapply(...) )
        return( parallel::mclapply(..., mc.cores=cores) )
    }
}

############################################################
# kinship_util
############################################################

# returns vector of IDs
check_kinship <-
    function(kinship, n_chr)
{
    if(is_kinship_decomposed(kinship)) {
        if(is_kinship_list(kinship)) {
            if(length(kinship) != n_chr)
                stop("length(kinship) != no. chromosomes (", n_chr, ")")
            return(rownames(kinship[[1]]$vectors))
        }
        else {
            return(rownames(kinship$vectors))
        }
    }

    if(!is_kinship_list(kinship)) { # one kinship matrix
        stopifnot(nrow(kinship) == ncol(kinship))
        stopifnot( all(rownames(kinship) == colnames(kinship)) )
        return(rownames(kinship))
    } else {
        if(length(kinship) != n_chr)
            stop("length(kinship) != no. chromosomes (", n_chr, ")")

        kinship_square <- vapply(kinship, function(mat) nrow(mat) == ncol(mat), TRUE)
        stopifnot( all(kinship_square) )

        kinship_sameIDs <- vapply(kinship, function(mat) (nrow(mat) == nrow(kinship[[1]])) &&
                                  all((rownames(mat) == rownames(kinship[[1]])) &
                                      (colnames(mat) == colnames(kinship[[1]])) &
                                      (rownames(mat) == colnames(kinship[[1]]))), TRUE)
        if(!all(kinship_sameIDs))
            stop("All kinship matrices should be the same size ",
                 "and have the same row and column names")

        return(rownames(kinship[[1]]))
    }
}

# multiply kinship matrix by 2
# see Almasy & Blangero (1998) https://doi.org/10.1086/301844
#
# This can also handle the case of "loco", and of having eigen decomposition pre-computed
double_kinship <-
    function(kinship)
{
    if(is.null(kinship)) return(NULL)

    if(is_kinship_decomposed(kinship)) { # already did decomposition
        if(is_kinship_list(kinship)) { # list of decomposed kinship matrices
            kinship <- lapply(kinship, function(a) { a$values <- 2*a$values; a })
        }
        else { # one decomposed kinship matrix
            kinship$values <- 2*kinship$values
        }
    }
    else {
        if(is.list(kinship))
            kinship <- lapply(kinship, function(a) a*2)
        else kinship <- 2*kinship
    }

    kinship
}


# check if alread decomposed
is_kinship_decomposed <-
    function(kinship)
{
    decomp <- attr(kinship, "eigen_decomp")

    (!is.null(decomp) && decomp) || # should have attribute
        (length(kinship)==2 && all(names(kinship) == c("values", "vectors"))) || # single-chr case missing attribute
        (is.list(kinship) && all(vapply(kinship, length, 1)==2) && all(vapply(kinship, function(a) all(names(a)==c("values", "vectors")), TRUE))) # multi-chr case
}

# is kinship a list with (potentially) multiple chromosomes
is_kinship_list <-
    function(kinship)
{
    if(is_kinship_decomposed(kinship)) {
        if(length(kinship)==2 && all(names(kinship) == c("values", "vectors"))) { # just one chromosome
            return(FALSE)
        }
        else return(TRUE)
    }
    else {
        return(is.list(kinship))
    }
}

# check that kinship concerns one chromosome
# and remove outer list if necessary
check_kinship_onechr <-
    function(kinship)
{
    if(is_kinship_list(kinship)) {
        if(length(kinship) > 1)
            stop("kinship should concern just one chromosome")
        decomp <- attr(kinship, "eigen_decomp") ## preserve attribute
        kinship <- kinship[[1]]
        attr(kinship, "eigen_decomp") <- decomp
    }

    kinship
}

# multiply kinship by weights (from and back)
# assuming weights are really square-root weights
#' @importFrom stats setNames
weight_kinship <-
    function(kinship, weights=NULL, tol=1e-8)
{
    # if null weights are all 1's, just skip it
    if(is.null(weights) || max(abs(weights-1)) < tol) return(kinship)

    if(is_kinship_list(kinship)) {
        for(i in seq_along(kinship)) {
            kinship[[i]] <- weight_kinship(kinship[[i]], weights)
        }
        return(kinship)
    }

    # if kinship was decomposed, expand it and then decompose it later
    do_decomp <- FALSE
    if(is_kinship_decomposed(kinship)) {
        do_decomp <- TRUE
        # expand out the decomposition
        kinship <- t(kinship$vectors) %*% diag(kinship$values) %*% kinship$vectors
    }

    # line them up
    ind2keep <- get_common_ids(setNames(rownames(kinship), NULL), setNames(names(weights), NULL))
    weights <- weights[ind2keep]
    kinship <- kinship[ind2keep, ind2keep, drop=FALSE]

    # multiply kinship matrix by weights
    kinship <- kinship * outer(weights, weights, "*")

    # decompose kinship again
    if(do_decomp) kinship <- decomp_kinship(kinship)

    kinship
}

############################################################
# get_common_ids.R
############################################################
#' Get common set of IDs from objects
#'
#' For a set objects with IDs as row names (or, for a vector, just
#' names), find the IDs that are present in all of the objects.
#'
#' @param ... A set of objects: vectors, lists, matrices, data frames,
#' and/or arrays. If one is a character vector with no names
#' attribute, it's taken to be a set of IDs, itself.
#' @param complete.cases If TRUE, look at matrices and non-character
#' vectors and keep only individuals with no missing values.
#'
#' @return A vector of character strings for the individuals that are
#' in common.
#'
#' @details This is used (mostly internally) to align phenotypes,
#' genotype probabilities, and covariates in preparation for a genome
#' scan. The `complete.cases` argument is used to omit
#' individuals with any missing covariate values.
#'
#' @examples
#' x <- matrix(0, nrow=10, ncol=5); rownames(x) <- LETTERS[1:10]
#' y <- matrix(0, nrow=5, ncol=5);  rownames(y) <- LETTERS[(1:5)+7]
#' z <- LETTERS[5:15]
#' get_common_ids(x, y, z)
#'
#' x[8,1] <- NA
#' get_common_ids(x, y, z)
#' get_common_ids(x, y, z, complete.cases=TRUE)
#'
#' @export
get_common_ids <-
    function(..., complete.cases=FALSE)
{
    args <- list(...)
    if(length(args)==0) {
        return(character(0))
    }

    # find the IDs in common across all
    id <- NULL
    for(i in seq_along(args)) {
        if(is.null(args[[i]])) next

        if(is.matrix(args[[i]]) || is.data.frame(args[[i]]) || is.array(args[[i]])) {
            if(length(dim(args[[i]])) > 3)
                stop("Can't handle arrays with >3 dimensions")
            these <- rownames(args[[i]])
            if(complete.cases && (is.matrix(args[[i]]) || is.data.frame(args[[i]])))
                these <- these[rowSums(!is.finite(args[[i]]))==0]
        }
        else if(is.list(args[[i]]) && !is.null(rownames(args[[i]][[1]]))) {
            these <- rownames(args[[i]][[1]])
        }
        else if(is.vector(args[[i]])) {
            if(is.character(args[[i]])) {
                if(is.null(names(args[[i]]))) {
                    these <- args[[i]]
                } else {
                    these <- names(args[[i]])
                    if(complete.cases) {
                        these <- these[!is.na(args[[i]])]
                    }
                }
            }
            else {
                these <- names(args[[i]])
                if(complete.cases) {
                    these <- these[is.finite(args[[i]])]
                }
            }
        }
        else if(is.character(args[[i]])) { # character but not vector
            if(is.null(names(args[[i]]))) {
                these <- args[[i]]
            } else {
                these <- names(args[[i]])
                if(complete.cases) {
                    these <- these[!is.na(args[[i]])]
                }
            }
        }
        else if(!is.null(names(args[[i]]))) { # not a vector but has names
            these <- names(args[[i]])
        }
        else {
            stop("Not sure what to do with object of class ", class(args[[i]]))
        }

        if(length(unique(these)) != length(these))
            stop("Duplicate names in argument ", i)

        if(is.null(id)) id <- these
        else id <- id[id %in% these]
    }

    names(id) <- NULL # just in case; should be bare vector of character strings
    id
}
