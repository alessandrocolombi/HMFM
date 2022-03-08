#ifndef COMPONENTPRIOR_FACTORY_H
#define COMPONENTPRIOR_FACTORY_H


#include "include_headers.h"
#include "recurrent_traits.h"
#include "ComponentPrior.h"

// struct to collect all possibile parametes of the different types of priors. It is needed to handle distributions that have a different number of parameters
// that may have different types. 
struct ComponentPrior_Parameters{
	// For Poisson distribution
	double Lambda{1}; //Poisson parameter, mean and variance are equal to Lambda
	// For Negative Binomial Distribution
	double p{0.5}; // probability of successes of each experiment
	unsigned int n_succ{1}; // number of successes after which the experiments stops. It is not required to be integer
};


// Factory to create the polymorphic object
enum class Mpriors{
    Poi1,
    NegBin1
};

template< Mpriors qM = Mpriors::Poi1, typename... Args  >
std::unique_ptr< ComponentPrior > Create_qM(Args&&... args){

    static_assert(qM == Mpriors::Poi1 ||  qM == Mpriors::NegBin1,
                  "Error, only possible priors for the number of components are Poi1 and NegBin1");
    if constexpr(qM == Mpriors::Poi1)
        return std::make_unique< Poisson1 >(std::forward<Args>(args)...); 
    else if constexpr(qM == Mpriors::NegBin1)
        return std::make_unique< NegativeBinomial1 >(std::forward<Args>(args)...); 
    //else if constexpr(algo == GGMAlgorithm::DRJ)
        //return std::make_unique< DoubleReversibleJumpsMH<GraphStructure, T> > (std::forward<Args>(args)...);

}


// Creates the actual object using the string passed in the main
std::unique_ptr< ComponentPrior > 
Select_ComponentPrior(std::string const & namePr, ComponentPrior_Parameters const & qM_params){
	// Select Prior
	if(namePr == "Poisson")
		return Create_qM<Mpriors::Poi1>(namePr, qM_params.Lambda);
	else if(namePr == "NegativeBinomial")
		return Create_qM<Mpriors::NegBin1>(namePr, qM_params.p, qM_params.n_succ);
	else
		throw std::runtime_error("Error, only possible priors for the number of components are Poi1 and NegBin1");
}

#endif