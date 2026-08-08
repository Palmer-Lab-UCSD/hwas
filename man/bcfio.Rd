% File src/library/base/man/load.Rd
\name{bcfio}
\alias{bcfio}
\title{Parse vcf, vcf.gz, and bcf files}
\description{
    The \code{vcf}, \code{vcf.gz}, and \code{bcf} file formats,
    collectively referred here has \code{hts} formats, are standard
    in bioinformatics for storing genomic data.  Programmatic 
    access to data in these formats use the C library \code{htslib}.
    These functions create interface between R and \code{htslib} to
    give high level access in R to genomic data.
}
\usage{
hwas::bopen(filename, mode)
hwas::bclose(bid)
hwas::is_open(bid)
hwas::is_bcf(filename)
hwas::num_positions(bid)
hwas::num_samples(bid)
hwas::sample_names(bid)
hwas::k_fmt(bid, id)
hwas::subset_samples(bid, samplesFile)
hwas::set_threads(bid, nthreads)
}
\arguments{
  \item{filename}{a string containing the path and name of \code{hts}
      file to establish a connection.}
  \item{mode}{the mode in which the file is open, currently read only
      (\code{r}) is supported.}
  \item{bid}{a file connection object.  Specifically, a \code{Rcpp}
    nullable pointer to the C++ class instance managing the file
    connection.  \code{bid} is \code{NULL} if establishing a file
    connection fails.}
  \item{id}{the data \code{hts} ``FORMAT'' ``id'' to be retrieved.  
    The ``FORMAT'' is defined in the \code{VCF} file specification.}
  \item{samplesFile}{file containing newline delimited names of samples 
    that should be included when querying data.}
  \item{nthreads}{the number of threads to divide data processing over.}
}
\value{
  A character vector of the names of objects created, invisibly.
}
%\seealso{
%  \code{\link{save}}.
%}
%\examples{
%hwas::init(geno_path, )
%}
\keyword{file}
