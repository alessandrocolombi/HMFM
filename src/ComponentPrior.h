#ifndef GDFMM_COMPONENTPRIOR_H
#define GDFMM_COMPONENTPRIOR_H


#include "include_headers.h"
#include "recurrent_traits.h"
#include "GSL_wrappers.h"

class ComponentPrior{ 
public :
    ComponentPrior(const std::string& _name) : name(_name){};
    ComponentPrior() = default;
    std::string showMe() const;
    virtual std::unique_ptr<ComponentPrior> clone() const = 0;
    virtual unsigned int get_mode() const = 0;
    virtual double eval_prob(unsigned int const &  k) const = 0;
    virtual double log_eval_prob(unsigned int const & k) const = 0;
    virtual ~ComponentPrior() = default;
private :    
    std::string name;
};

class Poisson1 : public ComponentPrior{ 
public:
    Poisson1(const std::string& _name, const double& _Lambda) : ComponentPrior(_name), Lambda(_Lambda){
        if(Lambda <= 0)
            throw std::runtime_error("Error in Poisson1 constructor, Lambda parameter in Poisson distribution must be strictly positive");
    };
    std::unique_ptr<ComponentPrior> clone() const override{
        return std::make_unique<Poisson1>(*this);
    };
    unsigned int get_mode() const override;
    double eval_prob(unsigned int const & k) const override;
    double log_eval_prob(unsigned int const & k) const override;
private:    
    double Lambda;
};


// There is not an unique definition of the negative binomial distribution. 
// Here, we follow the one of the gls library,
// "This function returns a random integer from the negative binomial distribution, the number of failures occurring before n successes in independent trials with probability p of success."
// and it is shifted by 1, the probability of k=0 is 0 and the probability of k is equal to the one computed using the gsl function with k = k-1. 
class NegativeBinomial1 : public ComponentPrior{ 
public:
    NegativeBinomial1(const std::string& _name, const double& _p, const unsigned int& _n) : ComponentPrior(_name), p(_p), n(_n){
        if(n <= 0 || p<0 || p>1)
            throw std::runtime_error("Error in NegativeBinomial1 constructor, n parameter has to be strictly positive and integer and p has to be in [0,1]");
    };
    std::unique_ptr<ComponentPrior> clone() const override{
        return std::make_unique<NegativeBinomial1>(*this);
    };
    unsigned int get_mode() const override;
    double eval_prob(unsigned int const & k) const override;
    double log_eval_prob(unsigned int const & k) const override;
private:    
    double p;
    unsigned int n;
};



#endif
