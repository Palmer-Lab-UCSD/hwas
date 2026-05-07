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

#' Calculate eigen decomposition of kinship matrix
#'
#' Calculate the eigen decomposition of a kinship matrix, or of a list of such matrices.
#'
#' @param kinship A square matrix, or a list of square matrices.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return The eigen values and the **transposed** eigen vectors,
#' as a list containing a vector `values` and a matrix
#' `vectors`.
#'
#' @details The result contains an attribute `"eigen_decomp"`.
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[1:30,18:19] # subset to 30 individuals and two chromosomes}
#' map <- insert_pseudomarkers(iron$gmap, step=1)
#' probs <- calc_genoprob(iron, map, error_prob=0.002)
#' K <- calc_kinship(probs)
#'
#' Ke <- decomp_kinship(K)
#'
#' @export
decomp_kinship <-
    function(kinship, cores=1)
{
    # already done?
    if(is_kinship_decomposed(kinship))
        return(kinship) # no need to do it again

    if(is.matrix(kinship)) {
        if(ncol(kinship) != nrow(kinship))
            stop("matrix must be square")
        if(ncol(kinship) == 0)
            stop("matrix has dimension (0,0)")
        return(Rcpp_eigen_decomp(kinship))
    }

    if(!is.list(kinship))
        stop("kinship should be either a square matrix or a list of square matrices")

    cores <- setup_cluster(cores)

    result <- cluster_lapply(cores, kinship, Rcpp_eigen_decomp)
    attr(result, "eigen_decomp") <- TRUE
    result
}
