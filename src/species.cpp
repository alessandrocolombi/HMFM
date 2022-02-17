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

//' Title Rcpp function
//'
//' @export
// [[Rcpp::export]]
int try_rcpp(int x){
	Rcpp::Rcout<<"Inside test function"<<std::endl;
	return x+10;
}

//' Raising Factorial
//'
//' \loadmathjax This function computes the rising factorial \mjseqn{(a)^{n}} using the gsl code for the Pochhammer symbol, i.e
//' \mjsdeqn{(a)^{n} = \frac{\Gamma(a+n)}{\Gamma(a)}}.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
//' @export
// [[Rcpp::export]]
double raising_factorial(const unsigned int& n, const double& a)
{
	return( gsl_sf_poch(a, (double)n ) );
}

//' Falling Factorial
//'
//' \loadmathjax This function computes the falling factorial \mjseqn{ a_{n} }. See \code{\link{raising_factorial}} for details.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
//' @export
// [[Rcpp::export]]
double my_falling_factorial(const unsigned int& n, const double& a)
{
	if(n%2 == 0) //n is even
		return( gsl_sf_poch(-a, (double)n ) );
	else //n is odd, change sign
		return( -1*gsl_sf_poch(-a, (double)n ) );
}

//' Pochhammer Symbol
//'
//' \loadmathjax This function computes the Pochhammer symbol, 
//' \mjsdeqn{(a)^{x} = \frac{\Gamma(a+x)}{\Gamma(a)}}.
//' Where \code{x} is a real number. When x is an integer, such a function coincides with the rising factorial defined in \code{\link{raising_factorial}}.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
//' @export
// [[Rcpp::export]]
double compute_Pochhammer(const unsigned int& x, const double& a)
{
	return( gsl_sf_poch(a,x ) );
}


//' Build matrix of logC numbers
//'
//' This is the recursive function called by \code{\link{my_logC}} that builds the matrix containing all the log(|C(n,k)|) numbers.
//' It gets as input the element (n,k) to build, the scale s and location r (here defined as positive and non-negative numbers) and the
//' matrix res to be build. The matrix is constructed diagonal by diagonal, starting from the bottom.
//' Important remark, note that log(|C(n,0)|) = log(raising factorial (n,r)). 
//'
// [[Rcpp::export]]
void build_logC_matrix(const unsigned int& n, const unsigned int& k, const double& s, const double& r, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& res){

	//Rcpp::Rcout<<"Costruisco n = "<<n<<", k = "<<k<<std::endl;
	if(k>n)
		throw std::runtime_error("Can not call build_logC_matrix if k > n");

	// Boundary conditions
	if(n == 0 && k == 0)
		res(n,k) = 0;
	else if(k == 0) 
		res(n,k) = std::log(raising_factorial(n,r));
	else if(n == k){
		build_logC_matrix(n-1, k-1, s, r, res);
		res(n,k) = std::log(s) + res(n-1,k-1);
	}
	else{
		double coef(s*k + r + n - 1);
		build_logC_matrix(n-1,k-1,s,r,res); // the matrix is always constructed diagonal by diagonal
		res(n,k) = std::log(coef) + res(n-1,k) + std::log( 1 + s/coef * std::exp( res(n-1,k-1) - res(n-1,k) ) );
	}
}

//' My logC
//'
//' This function is the recursive formula 2.69 in "Combinatorial methods in discrete distributions" book by Charalambides. 
//' It returns an (n+1 x n+1) matrix containing all the log(|C(nn,k)|) numbers, for nn = 0,...,n+1 and k = 0,...,nn. 
//' scale and location must be negative and non-positive, respectively.
//' As a consequence, it is memory expensive.
//' @export
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>	
my_logC(const unsigned int& n, const double& scale, const double& location){

	if(!( (scale<0) & (location<=0) ) )
		throw std::runtime_error("Error in my_logC. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative and location in non positive");
	
	const double& s = -scale; //s is strictly positive
	const double& r = -location; //r is non-negative

	MatCol res(MatCol::Constant(n+1,n+1,0.0)); //nxn matrix, n is in standard notation, i.e it starts from 1
	for(int k = n; k>=0; k--){
		//Rcpp::Rcout<<"Nel ciclo for, n = "<<n<<" , k = "<<k<<std::endl;
		build_logC_matrix(n,k,s,r,res);
		//Rcpp::Rcout<<"res finito k = "<<std::endl<<res<<std::endl;
	}
	return(res);
}


//' Compute log of absolute values of non Central C number 
//'
//' \loadmathjax This is the main function in the computation of C numbers. It uses the (2.69) formula in the "Combinatorial methods in discrete distributions" book.
//' It computes \mjseqn{log(|C(n,k; scale, location)|)} for each k=0,...,n. 
//' scale and location must be negative and non-positive, respectively.
//' It uses eigen objects, apparetly it is slower than using Rcpp vectors.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd my_logC2(const unsigned int& n, const double& scale, const double& location){

	if(!( (scale<0) & (location<=0) ) )
		throw std::runtime_error("Error in my_logC. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative and location in non positive");
	
	const double& s = -scale; //s is strictly positive
	const double& r = -location; //r is non-negative

	
	VecCol LogC_old(VecCol::Constant(n+1,0.0));

	if(n == 0)
		return LogC_old; // nothing to do in this case

	// Compute the first row
	LogC_old(0) = std::log(raising_factorial(1,r));
	LogC_old(1) = std::log(s);

	VecCol LogC_update(VecCol::Constant(n+1,0.0));
	double coef(0.0);
	for(std::size_t nn = 2; nn <= n; ++nn ){ //for each row
		LogC_update(0) = std::log(raising_factorial(nn,r));
		for(std::size_t k = 1; k < nn; ++k){ //for each column but the first and the last one
			coef = s*k + r + nn - 1;
			LogC_update(k) = std::log(coef) + LogC_old(k) + std::log( 1 + s/coef*std::exp( LogC_old(k-1) - LogC_old(k) ) );
		}
		LogC_update(nn) = nn*std::log(s); //update last element
		LogC_old.swap(LogC_update); //avoid copy, LogC_update has to be overwritten but in this way no copy is performed to update LogC_old.
	} 
	return (LogC_old);
}


//' My logC2 - central
//'
//' This function is the recursive formula 2.69 in "Combinatorial methods in discrete distributions" book by Charalambides for central numbers.
//' This is the specialization of \code{\link{my_logC2}} for central C numbers.
//' scale must be negative.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd my_logC2_central(const unsigned int& n, const double& scale){

	if(!( scale<0 ) )
		throw std::runtime_error("Error in my_logC2. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative.");
	
	const double& s = -scale; //s is strictly positive
	const double inf = std::numeric_limits<double>::infinity();

	VecCol LogC_old(VecCol::Constant(n+1,0.0));

	if(n == 0)
		return LogC_old; // nothing to do in this case

	// Compute the first row
	LogC_old(0) = -inf;
	LogC_old(1) = std::log(s);

	VecCol LogC_update(VecCol::Constant(n+1,0.0));
	LogC_update(0) = -inf;
	double coef(0.0);
	for(std::size_t nn = 2; nn <= n; ++nn ){ //for each row
		for(std::size_t k = 1; k < nn; ++k){ //for each column but the first and the last one
			coef = s*k + nn - 1;
			LogC_update(k) = std::log(coef) + LogC_old(k) + std::log( 1 + s/coef*std::exp( LogC_old(k-1) - LogC_old(k) ) );
		}
		LogC_update(nn) = nn*std::log(s); //update last element
		LogC_old.swap(LogC_update); //avoid copy, LogC_update has to be overwritten but in this way no copy is performed to update LogC_old.
	} 
	return (LogC_old);
}



//' compute_logC - Compute log of absolute values of non Central C number 
//'
//' \loadmathjax This is the main function in the computation of C numbers. It uses the (2.69) formula in the "Combinatorial methods in discrete distributions" book.
//' It computes \mjseqn{log(|C(n,k; scale, location)|)} for each k=0,...,n. 
//' This implementation uses Rcpp vectors.
//' @param scale must be strictly negative.
//' @param locatio must be non positive. Set to 0 for central C numbers.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector compute_logC(const unsigned int& n, const double& scale, const double& location){

	if(!( (scale<0) & (location<=0) ) )
		throw std::runtime_error("Error in my_logC. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative and location in non positive");
	
	const double& s = -scale; //s is strictly positive
	const double& r = -location; //r is non-negative

	
	Rcpp::NumericVector LogC_old(n+1, 0.0);

	if(n == 0)
		return LogC_old; // nothing to do in this case

	// Compute the first row
	LogC_old[0] = std::log(raising_factorial(1,r));
	LogC_old[1] = std::log(s);

	//Rcpp::NumericVector LogC_update(LogC_old);
	Rcpp::NumericVector LogC_update(n+1, 0.0);
	double coef(0.0);
	for(std::size_t nn = 2; nn <= n; ++nn ){ //for each row
		LogC_update[0] = std::log(raising_factorial(nn,r));
		for(std::size_t k = 1; k < nn; ++k){ //for each column but the first and the last one
			coef = s*k + r + nn - 1;
			LogC_update[k] = std::log(coef) + LogC_old[k] + std::log( 1 + s/coef*std::exp( LogC_old[k-1] - LogC_old[k] ) );
		}
		LogC_update[nn] = nn*std::log(s); //update last element
		
		//std::copy(LogC_update.begin(),LogC_update.end(),LogC_old.begin());
		std::swap(LogC_old, LogC_update); //avoid copy, LogC_update has to be overwritten but in this way no copy is performed to update LogC_old.
	} 
	return (LogC_old);
}

#endif
