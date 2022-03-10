#ifndef __SPECIES_HPP__
#define __SPECIES_HPP__

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>
#include "include_headers.h"
#include "recurrent_traits.h"
#include "GSL_wrappers.h"
#include "utils.h"
#include <gsl/gsl_sf.h>
#include "ComponentPrior_factory.h"

//' Title Rcpp function
//'
//' @export
// [[Rcpp::export]]
int try_rcpp(int x);

//' Raising Factorial
//'
//' \loadmathjax This function computes the rising factorial \mjseqn{(a)^{n}} using the gsl code for the Pochhammer symbol, i.e
//' \mjsdeqn{(a)^{n} = \frac{\Gamma(a+n)}{\Gamma(a)}}.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
//' @export
// [[Rcpp::export]]
double raising_factorial(const unsigned int& n, const double& a);

//' log Raising Factorial
//'
//' \loadmathjax This function computes the logarithm of the rising factorial \mjseqn{(a)^{n}} using the gsl code for the log of Pochhammer symbol.
//' See \code{\link{raising_factorial}} and \code{\link{compute_Pochhammer}} for details.
//' @export
// [[Rcpp::export]]
double log_raising_factorial(const unsigned int& n, const double& a);

//' Falling Factorial
//'
//' \loadmathjax This function computes the falling factorial \mjseqn{ a_{n} }. See \code{\link{raising_factorial}} for details.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
//' @export
// [[Rcpp::export]]
double my_falling_factorial(const unsigned int& n, const double& a);

//' log Falling Factorial
//'
//' \loadmathjax This function computes the logarithm of the falling factorial \mjseqn{ a_{n} } using the gsl code for the log of Pochhammer symbol.
//' See \code{\link{my_falling_factorial}} and \code{\link{compute_Pochhammer}} for details.
//' @export
// [[Rcpp::export]]
double my_log_falling_factorial(const unsigned int& n, const double& a);

//' Pochhammer Symbol
//'
//' \loadmathjax This function computes the Pochhammer symbol,
//' \mjsdeqn{(a)^{x} = \frac{\Gamma(a+x)}{\Gamma(a)}}.
//' Where \code{x} is a real number. When x is an integer, such a function coincides with the rising factorial defined in \code{\link{raising_factorial}}.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
//' @export
// [[Rcpp::export]]
double compute_Pochhammer(const unsigned int& x, const double& a);

//' Pochhammer log Symbol
//'
//' \loadmathjax This function computes the Pochhammer symbol in log form. See \code{\link{compute_Pochhammer}} for details.
//' @export
// [[Rcpp::export]]
double compute_log_Pochhammer(const unsigned int& x, const double& a);

//' Build matrix of logC numbers
//'
//' This is the recursive function called by \code{\link{my_logC}} that builds the matrix containing all the log(|C(n,k)|) numbers.
//' It gets as input the element (n,k) to build, the scale s and location r (here defined as positive and non-negative numbers) and the
//' matrix res to be build. The matrix is constructed diagonal by diagonal, starting from the bottom.
//' Important remark, note that log(|C(n,0)|) = log(raising factorial (n,r)).
//'
// [[Rcpp::export]]
void build_logC_matrix(const unsigned int& n, const unsigned int& k, const double& s, const double& r, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& res);


//' My logC
//'
//' This function is the recursive formula 2.69 in "Combinatorial methods in discrete distributions" book by Charalambides.
//' It returns an (n+1 x n+1) matrix containing all the log(|C(nn,k)|) numbers, for nn = 0,...,n+1 and k = 0,...,nn.
//' scale and location must be negative and non-positive, respectively.
//' As a consequence, it is memory expensive.
//' @export
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
my_logC(const unsigned int& n, const double& scale, const double& location);


//' Compute log of absolute values of non Central C number
//'
//' \loadmathjax This is the main function in the computation of C numbers. It uses the (2.69) formula in the "Combinatorial methods in discrete distributions" book.
//' It computes \mjseqn{log(|C(n,k; scale, location)|)} for each k=0,...,n.
//' scale and location must be negative and non-positive, respectively.
//' It uses eigen objects, apparetly it is slower than using Rcpp vectors.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd my_logC2(const unsigned int& n, const double& scale, const double& location);

//' My logC2 - central
//'
//' This function is the recursive formula 2.69 in "Combinatorial methods in discrete distributions" book by Charalambides for central numbers.
//' This is the specialization of \code{\link{my_logC2}} for central C numbers.
//' scale must be negative.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd my_logC2_central(const unsigned int& n, const double& scale);



//' compute_logC - Compute log of absolute values of non Central C number
//'
//' \loadmathjax This is the main function in the computation of C numbers. It uses the (2.69) formula in the "Combinatorial methods in discrete distributions" book.
//' It computes \mjseqn{log(|C(n,k; scale, location)|)} for each k=0,...,n.
//' This implementation uses Rcpp vectors.
//' @param scale must be strictly negative.
//' @param locatio must be non positive. Set to 0 for central C numbers.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector compute_logC(const unsigned int& n, const double& scale, const double& location);


// Da Testare
double compute_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max = 100 );



// Tutta da testare
double compute_log_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max = 100 );


//questa è sola per 1 o 2 gruppi
double compute_Kprior_unnormalized(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma);


//' 
// [[Rcpp::export]] 
double p_distinct_prior_c(const unsigned int& k, const Rcpp::NumericVector& n_groups, const Rcpp::NumericVector& gamma_groups, const Rcpp::String& prior, 
						  const Rcpp::List& prior_param, unsigned int M_max  );

//' Test ComponentPrior
//' @export
// [[Rcpp::export]]
void Test_Prior();


#endif
