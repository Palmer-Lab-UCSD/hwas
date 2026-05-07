# functions to deal with covariates
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

# force interactive covariates into the additive covariate matrix
#
# addcovar and intcovar are two matrices
# intcovar columns should all be within the addcovar columns
# tol is tolerance for determining matching columns
#' @importFrom stats complete.cases
force_intcovar <-
    function(addcovar=NULL, intcovar=NULL, tol=1e-12)
{
    if(is.null(intcovar)) # no intcovar, so return addcovar w/o change
        return(addcovar)

    if(is.null(addcovar)) # no addcovar, so just return intcovar
        return(intcovar)

    if(!is.matrix(addcovar)) addcovar <- as.matrix(addcovar)
    if(!is.matrix(intcovar)) intcovar <- as.matrix(intcovar)

    # IDs in both; omitting any individuals with missing values
    ids <- get_common_ids(addcovar[complete.cases(addcovar),,drop=FALSE],
                          intcovar[complete.cases(intcovar),,drop=FALSE])

    # look for matching columns, having reduced to common individuals
    full <- cbind(addcovar[ids,,drop=FALSE], intcovar[ids,,drop=FALSE])
    has_match <- find_matching_cols(full, tol)
    if(any(has_match > 0))
        full <- full[, has_match<0, drop=FALSE]

    if(ncol(full)==0) return(NULL)

    full
}

# drop linearly dependent columns
# if intercept=TRUE, add intercept before checking and then remove afterwards
#' @importFrom stats complete.cases
drop_depcols <-
    function(covar=NULL, add_intercept=FALSE, tol=1e-12)
{
    if(is.null(covar)) return(covar)

    if(!is.matrix(covar)) covar <- as.matrix(covar)
    if(add_intercept) covar <- cbind(rep(1, nrow(covar)), covar)

    if(ncol(covar) <= 1) return(covar)

    X <- covar[complete.cases(covar),,drop=FALSE]

    # deal with NAs by omitting those rows before
    indep_cols <- sort(find_lin_indep_cols(X, tol))

    if(add_intercept) {

        target_ncol <- length(indep_cols)

        n_it <- 0
        while(!(1 %in% indep_cols)) {
            # don't want to omit the first column (the intercept)
            # need to work harder...
            #  - drop one column at a time other the intercept
            #  - when you find a column that doesn't reduce the target number of columns, drop it
            #  - check again if the intercept is being retained; if not, repeat

            for(i in seq_len(ncol(X))[-1]) { # loop over all but the first column (the intercept)
                indep_cols <- find_lin_indep_cols(X[,-i,drop=FALSE], tol)
                if(length(indep_cols) == target_ncol) { # this one definitely dependent; omit it
                    X <- X[,-i,drop=FALSE]
                    break
                }
            }

            # new set of linearly independent columns
            indep_cols <- sort(find_lin_indep_cols(X, tol))

            n_it <- n_it+1 # count number of iterations and bail if it's large
            if(n_it > ncol(covar)+5) { # something is seriously messed up
                stop("Not above to find set of linearly independent columns")
            }
        }

        # now drop the intercept
        indep_cols <- indep_cols[-1]
    }
    if(length(indep_cols)==0) return(NULL)

    covar[, indep_cols, drop=FALSE]
}

# drop columns from X covariates that are already in addcovar
#' @importFrom stats complete.cases
drop_xcovar <-
    function(covar=NULL, Xcovar=NULL, tol=1e-12)
{
    if(is.null(Xcovar) || is.null(covar)) return(Xcovar)

    if(!is.matrix(covar)) covar <- as.matrix(covar)
    if(!is.matrix(Xcovar)) Xcovar <- as.matrix(Xcovar)

    # IDs in both; omitting any individuals with missing values
    ids <- get_common_ids(covar[complete.cases(covar),,drop=FALSE],
                          Xcovar[complete.cases(Xcovar),,drop=FALSE])

    # look for matching columns, having reduced to common individuals
    matches <- find_matching_cols(cbind(covar[ids,], Xcovar[ids,]), tol)[-seq_len(ncol(covar))]

    if(all(matches > 0)) return(NULL)

    # drop the columns with matches
    Xcovar[,matches<0,drop=FALSE]
}
