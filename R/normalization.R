
# rbiqt: rank based inverse normal quantile transformation
NORMS <- c("rbiqt", "zscore")
TYPES <- c("character", "float")


is_valid_type <- function(s) {
    return(s %in% TYPES)
}

is_norm <- function(s) {
    return(s %in% NORMS)
}

#' @title inverse normal quantile transformation
#'
#' @description Implements a deterministic rank based inverse normal
#'  quantile transformation accoring Equation 1 in ref [1],
#'
#'  y_i = \frac{r_i - c}{N - 2c + 1}
#'
#'  with r_i being the rank of sample i, and c a constant set by the 
#'  user.  Note that according to [1] values commonly take values:
#'  c = 1/3 Tukey 1962, c = 1/2 Bliss 1967, c = 0 van der Waerden 
#'  1952, and c = 3/8 Bloom 1958.  I have not vetted this, for example
#'  I was able to find reference to van der Waerden paper, but wasn't
#'  able to find the paper itself.  We set c = 0 to the default as there
#'  is both precedent in the literature and the labs methodology.
#'
#'  To compute the normalized data, we simply apply the inverse 
#'  distribution function of the standard normal,
#'
#'  z_i = \Phi^{-1} (y_i).
#'
#'  Important, when ties are present the resulting transformed data
#'  are not guaranteed to have equivalent quantiles as the standard
#'  normal.  Instead, I've empirically shown that the quantiles are
#'  correct on average.
#'
#'  The rank of ties are handled by the average method.  To illustrate
#'  this method consider the folllowing example. Suppose that data values
#'  x_{i-1} = x_i = x_{i+1} are all equal.  Let's assign a unique rank
#'  based upon the order of the data r_{i-1}, r_i = r_{i-1}+1, and 
#'  r_{i+1} = r_{i-1} + 2.  The final ranks of these samples are their
#'  average, that is: 
#'
#'  r_{i-1} = r_i = r_{i+1} = \bar{r}_i = 1/3(r_{i-1} + r_i + r_{i+1}).
#'
#' @references 
#'  [1] Beasley, Erickson, and Allison, "Rank-Based Inverse Normal
#'      Transformations are Increasingly Used, But are They Merited?" 
#'      Behav Genet. 2009.
#' 
#' @export
rbiqt <- function(x, cval = 0) {

    if(!is.null(dim(x)))
        return(err(numeric(0), "Requires vector input"))

    n <- length(x)
    if (n <= 0)
        return(err(numeric(0), "empty data"))

    if (2 * cval > n || cval < 0)
        return(err(numeric(0), "cval must be on interval [0, n/2]"))

    if (!is.numeric(x))
        return(err(numeric(0), "requires a numeric vector input"))

    if (any(is.na(x)) || any(is.nan(x)))
        return(err(numeric(0), 
                   "NA or nan value(s) detected and unsupported."))

    r <- rank(x, ties.method = "average")

    return(qnorm((r - cval) / (n - 2*cval + 1)))
}


#'
#' @export
zscore <- function(x) {
    if(!is.null(dim(x)))
        return(err(numeric(0), "Requires vector input"))

    if (length(x) <= 0)
        return(err(numeric(0), "empty data"))

    if (!is.numeric(x))
        return(err(numeric(0), "requires a numeric vector input"))

    if (any(is.na(x)) || any(is.nan(x)))
        return(err(numeric(0), 
                   "NA or nan value(s) detected and unsupported."))

    return((x - mean(x)) / sd(x))
}
