#ifndef __GDFMMEXPORT_HPP__
#define __GDFMMEXPORT_HPP__

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>

//' Title Rcpp function
//'
//' @export
// [[Rcpp::export]]
int try_rcpp(int x){
	Rcpp::Rcout<<"Inside test function"<<std::endl;
	return x+10;
}

#endif
