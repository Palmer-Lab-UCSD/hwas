\name{pg_sim}
\alias{pg_sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Simulate phenotype data under the haplotype polygenic model.
}
\description{
    The polygenic model is a quantitative description of phenotype
    as the summation of independent random genetic and environment 
    effects.
}
\usage{
pg_sim(bcf_filename, heritability)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bcf_filename}{Name 
%%     ~~Describe \code{object} here~~
}
  \item{heritability}{
      Heritability of the trait under the haplotype random effects
      model.
}
  \item{seed}{Random number generator seed (default NULL)}
}
\details{
    distinct from modeling a quantitative trait locus, where a 
    the QTL effect is fixed.  Mathematically this distinction
    between modeling genetic effects on phenotypes as random or
    fixed effects 
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (object, nsim = 1, seed = NULL, ...) 
UseMethod("simulate")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
