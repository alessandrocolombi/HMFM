#ifndef __GDFMMEXPORT_HPP__
#define __GDFMMEXPORT_HPP__

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>

using MatRow        = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatCol        = Eigen::MatrixXd;
using MatUnsRow     = Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatUnsCol     = Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using VecUnsCol     = Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>;
using VecCol        = Eigen::VectorXd;
using VecRow        = Eigen::RowVectorXd;

//' Title Rcpp test function
//'
//' @export
// [[Rcpp::export]]
int try_Rcpp(int x){
  Rcpp::Rcout<<"Inside first test function. x = "<<x<<std::endl;
  return (x+10);
}

#endif
