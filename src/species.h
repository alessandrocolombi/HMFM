#ifndef __SPECIES_HPP__
#define __SPECIES_HPP__

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppParallel)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppGSL.h>
#include "include_headers.h"
#include "recurrent_traits.h"
#include "GSL_wrappers.h"
//#include "utils.h"
#include <gsl/gsl_sf.h>
#include "ComponentPrior_factory.h"

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	log_stable_sum
//------------------------------------------------------------------------------------------------------------------------------------------------------

// These functions computes log(sum_i(a_i)) using a stable formula for log values. Let us denote a* to the the maximum value of vector a which is attained when i = i*.
// ---> log(sum_i(a_i)) = log(a*) + log[ 1 + sum_{i not i*}(exp{log(a_i) - log(a*)}) ]
// See that only the logarithm of the elements of a are needed. Hence, it is likely that one has already computed them in log scale. If so, set is_log = T
// In this version of the function, the maximum and its position are passed as parameters. No check is performed to be sure that such information were true or not.
// It is assumed that the value in val_max is on the same scale of the values of a, i.e it is in log scale if is_log is set to true.
double log_stable_sum(const std::vector<double>& a, const bool is_log, const double& val_max, const unsigned int& idx_max);

// As before but gets an iterator poining to the maximum value
double log_stable_sum(const std::vector<double>& a, const bool is_log, const GDFMM_Traits::vector_d_citerator& it_max);

// In this specialized version of the function, the position of the max value is not provided. Hence, one additional operation is done. 
// Since exp( log(a*) - log(a*)) = 1, the previous formula becomes 
// ---> log(sum_i(a_i)) = log(a*) + log[ sum_{i}(exp{log(a_i) - log(a*)}) ]
double log_stable_sum(const std::vector<double>& a, const bool is_log, const double& val_max);

// In this version of the formula, the maximum value is computed
double log_stable_sum(const std::vector<double>& a, const bool is_log);



//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Factorials and Pochammer
//------------------------------------------------------------------------------------------------------------------------------------------------------


//' Raising Factorial
//'
//' \loadmathjax This function computes the rising factorial \mjseqn{(a)^{n}} using the gsl code for the Pochhammer symbol, i.e
//' \mjsdeqn{(a)^{n} = \frac{\Gamma(a+n)}{\Gamma(a)}}.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
double raising_factorial(const unsigned int& n, const double& a);

//' log Raising Factorial
//'
//' \loadmathjax This function computes the logarithm of the rising factorial \mjseqn{(a)^{n}} using the gsl code for the log of Pochhammer symbol.
//' See \code{\link{raising_factorial}} and \code{\link{compute_Pochhammer}} for details.
double log_raising_factorial(const unsigned int& n, const double& a);


//' Raising Factorial
//'
//' \loadmathjax This function computes the rising factorial \mjseqn{(a)^{n}} using the gsl code for the Pochhammer symbol, i.e
//' \mjsdeqn{(a)^{n} = \frac{\Gamma(a+n)}{\Gamma(a)}}.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
double raising_factorial_poch(const unsigned int& n, const double& a);

//' log Raising Factorial
//'
//' \loadmathjax This function computes the logarithm of the rising factorial \mjseqn{(a)^{n}} using the gsl code for the log of Pochhammer symbol.
//' See \code{\link{raising_factorial}} and \code{\link{compute_Pochhammer}} for details.
double log_raising_factorial_poch(const unsigned int& n, const double& a);


//' Falling Factorial
//'
//' \loadmathjax This function computes the falling factorial \mjseqn{ a_{n} }. See \code{\link{raising_factorial}} for details.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
double my_falling_factorial(const unsigned int& n, const double& a);

//' log Falling Factorial
//'
//' \loadmathjax This function computes the logarithm of the falling factorial \mjseqn{ a_{n} } using the gsl code for the log of Pochhammer symbol.
//' See \code{\link{my_falling_factorial}} and \code{\link{compute_Pochhammer}} for details.
double my_log_falling_factorial(const unsigned int& n, const double& a);

//' Falling Factorial
//'
//' \loadmathjax This function computes the falling factorial \mjseqn{ a_{n} }. See \code{\link{raising_factorial}} for details.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
double my_falling_factorial_old(const unsigned int& n, const double& a);

//' log Falling Factorial
//'
//' \loadmathjax This function computes the logarithm of the falling factorial \mjseqn{ a_{n} } using the gsl code for the log of Pochhammer symbol.
//' See \code{\link{my_falling_factorial}} and \code{\link{compute_Pochhammer}} for details.
double my_log_falling_factorial_old(const unsigned int& n, const double& a);

//' Pochhammer Symbol
//'
//' \loadmathjax This function computes the Pochhammer symbol,
//' \mjsdeqn{(a)^{x} = \frac{\Gamma(a+x)}{\Gamma(a)}}.
//' Where \code{x} is a real number. When x is an integer, such a function coincides with the rising factorial defined in \code{\link{raising_factorial}}.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
double compute_Pochhammer(const unsigned int& x, const double& a);

//' Pochhammer log Symbol
//'
//' \loadmathjax This function computes the Pochhammer symbol in log form. See \code{\link{compute_Pochhammer}} for details.
double compute_log_Pochhammer(const unsigned int& x, const double& a);


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	C numbers
//------------------------------------------------------------------------------------------------------------------------------------------------------



//' Build matrix of logC numbers
//'
//' This is the recursive function called by \code{\link{my_logC}} that builds the matrix containing all the log(|C(n,k)|) numbers.
//' It gets as input the element (n,k) to build, the scale s and location r (here defined as positive and non-negative numbers) and the
//' matrix res to be build. The matrix is constructed diagonal by diagonal, starting from the bottom.
//' Important remark, note that log(|C(n,0)|) = log(raising factorial (n,r)).
void build_logC_matrix(const unsigned int& n, const unsigned int& k, const double& s, const double& r, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& res);


//' My logC
//'
//' This function is the recursive formula 2.69 in "Combinatorial methods in discrete distributions" book by Charalambides.
//' It returns an (n+1 x n+1) matrix containing all the log(|C(nn,k)|) numbers, for nn = 0,...,n+1 and k = 0,...,nn.
//' scale and location must be negative and non-positive, respectively.
//' As a consequence, it is memory expensive.
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>
my_logC(const unsigned int& n, const double& scale, const double& location);


//' Compute log of absolute values of non Central C number
//'
//' \loadmathjax This is the main function in the computation of C numbers. It uses the (2.69) formula in the "Combinatorial methods in discrete distributions" book.
//' It computes \mjseqn{log(|C(n,k; scale, location)|)} for each k=0,...,n.
//' scale and location must be negative and non-positive, respectively.
//' It uses eigen objects, apparetly it is slower than using Rcpp vectors.
Eigen::VectorXd my_logC2(const unsigned int& n, const double& scale, const double& location);

//' My logC2 - central
//'
//' This function is the recursive formula 2.69 in "Combinatorial methods in discrete distributions" book by Charalambides for central numbers.
//' This is the specialization of \code{\link{my_logC2}} for central C numbers.
//' scale must be negative.
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


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	A priori functions
//------------------------------------------------------------------------------------------------------------------------------------------------------

// Da Testare
double compute_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max = 100 );

// Tutta da testare
double compute_log_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max = 100 );

//Direct formula per d=1 or d=2
double compute_Kprior_unnormalized(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma);

//Recursive formula for d>2
double compute_Kprior_unnormalized_recursive(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma);

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Rcpp call functions
//------------------------------------------------------------------------------------------------------------------------------------------------------

// The same code was repeated multiple times in different functions. Hence it has been collected in a single function.
// Wrap the call for the ComponentPrior from Rcpp objects to c++ object
std::unique_ptr< ComponentPrior > Wrapper_ComponentPrior(const Rcpp::String& prior, const Rcpp::List& prior_param);

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Rcpp call functions for computing probabilities
//------------------------------------------------------------------------------------------------------------------------------------------------------

//' 
// [[Rcpp::export]] 
double p_distinct_prior_c(const unsigned int& k, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, 
							  const Rcpp::List& prior_param, unsigned int M_max  );

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Rcpp call functions for computing expected values
//------------------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Tests
//------------------------------------------------------------------------------------------------------------------------------------------------------


// This function computes prod_{i=1}^{n}( f(a_i*b_i) )
// Questo è solo un test, sarebbe carinio farla che prende in input anche la funzione f() generica da applicare
double combined_product(const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const unsigned int& Mstar, const unsigned int& k){

	return (
				std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 1.0, std::multiplies<>(),
	    					   		[&Mstar, &k](const double& nj, const double& gamma_j){return (  (double)nj + gamma_j*(double)(Mstar + k)  );} 
	    					   	  )
			);
}

// This function computes sum_{i=1}^{n}( f(ni_i*gamma_i) )
// Questo è solo un test, sarebbe carinio farla che prende in input anche la funzione f() generica da applicare
double combined_sum(const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const unsigned int& Mstar, const unsigned int& k){

	return (
				std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
	    					   		[&Mstar, &k](const double& nj, const double& gamma_j){return (  (double)nj + gamma_j*(double)(Mstar + k)  );} 
	    					   	  )
			);
}



// Test ComponentPrior
void Test_Prior();


void Test_prod_sum();


#endif
