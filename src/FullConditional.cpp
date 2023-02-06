#include "FullConditional.h"
bool FullConditional::binary_decision(double p1, const sample::GSL_RNG& engine) const{
    // Sample u from a Uniform(0,1) and verify if is less than p
    sample::runif unif;
    double u = unif(engine);
    return u<p1;
}

void FullConditional::print() const{
  std::cout<<name<<"\n";
}

double FullConditional::log_raising_factorial(const unsigned int& n, const double& a)const
{

  if(n==0)
    return 0.0;
  if(a<0)
    throw std::runtime_error("Error in my_log_raising_factorial, can not compute the raising factorial of a negative number in log scale"); 
  else if(a==0.0){
    return -std::numeric_limits<double>::infinity();
  }
  else{
    
    double val_max{std::log(a+n-1)};
    double res{1.0};
    if (n==1)
      return val_max;
    for(std::size_t i = 0; i <= n-2; ++i){
      res += std::log(a + (double)i) / val_max;
    }
    return val_max*res;

  }
}