
#include<Rcpp.h>

// [[Rcpp::export]]
Rcpp::List hello_world() {
    Rcpp::List output = Rcpp::List::create("hello world", 42);
    return output;
}

