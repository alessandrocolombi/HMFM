#include "ComponentPrior.h"

std::string ComponentPrior::showMe() const{
	return name;
}

double Poisson1::eval_prob(unsigned int const &  k) const{
	if(k == 0)
		return 0.0;
	else
		return gsl_ran_poisson_pdf(k-1, Lambda);
	//double res(std::exp(-Lambda)/gsl_sf_fact(k-1) );
	//for(std::size_t i = 1; i <= k-1; ++i)
		//res *= Lambda; 
	//return(res); 
	
} 

double Poisson1::log_eval_prob(unsigned int const & k) const{
	if(k == 0)
		return -std::numeric_limits<double>::infinity();
	double res(-Lambda + (k-1)*std::log(Lambda) );
	for(std::size_t i = 2; i <= k-1; ++i)
		res += std::log(i); 
	return(res); 
	
} 

double NegativeBinomial1::eval_prob(unsigned int const &  k) const{
	if(k == 0)
		return 0.0;
	else
		return gsl_ran_negative_binomial_pdf(k-1, p, n);
	
} 

double NegativeBinomial1::log_eval_prob(unsigned int const & k) const{
	if(k == 0)
		return -std::numeric_limits<double>::infinity();
	else{
		if( p==0 || p == 1)
			throw std::runtime_error("Error in NegativeBinomial1::log_eval_prob. It is not possible to compute the log probability if p is equal to 1 or 0");
		return ( std::lgamma(n+(double)k-1.0) - std::lgamma(k+1) - std::lgamma(n) + n*std::log(p) + (double)(k-1)*std::log(1-p) ); 
	}
	
} 