#/bin/bash
#
#

R -e "library(Rcpp); Rcpp::compileAttributes()"

R CMD build .
R CMD check "hwas_0.0.1.tar.gz"
