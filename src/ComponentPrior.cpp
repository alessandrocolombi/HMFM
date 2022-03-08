#include "ComponentPrior.h"

std::string ComponentPrior::showMe() const{
	return name;
}

unsigned int Poisson1::get_mode() const{
	return (static_cast<unsigned int>( std::floor(Lambda) ) + 1);
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
		return -std::numeric_limits<double>::infinity(); //sbagliata per k>1!!
	else
		return (-Lambda + (k-1)*std::log(Lambda) - gsl_sf_lnfact(k-1) );
} 

unsigned int NegativeBinomial1::get_mode() const{
	if(n <= 1)
		return 1;
	else
		return ( static_cast<unsigned int>( std::floor(p*(n-1)/(1-p)) ) + 1 );
}

double NegativeBinomial1::eval_prob(unsigned int const &  k) const{
	if(k == 0)
		return 0.0;
	else
		return gsl_ran_negative_binomial_pdf(k-1, p, (double)n);
	
} 

double NegativeBinomial1::log_eval_prob(unsigned int const & k) const{
	if(k == 0)
		return -std::numeric_limits<double>::infinity();
	else{
		
		if( p==0 || p == 1)
			throw std::runtime_error("Error in NegativeBinomial1::log_eval_prob. It is not possible to compute the log probability if p is equal to 1 or 0");

		return ( gsl_sf_lnchoose(k+n-2, k-1) + (double)n*std::log(p) + (double)(k-1)*std::log(1-p) );
		
	}
	
} 