#ifndef __GDFMMEXPORT_HPP__
#define __GDFMMEXPORT_HPP__

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>
#include "include_headers.h"
#include "recurrent_traits.h"
#include "GSL_wrappers.h"
#include "GibbsSampler.h"



//' Title Rcpp function
//'
//' @export
// [[Rcpp::export]]
int try_rcpp(int x){
	Rcpp::Rcout<<"Inside test function"<<std::endl;
	return x+10;
}


// [[Rcpp::export]]

Rcpp::List example_GDFMM_sampler_c( Eigen::MatrixXd const & dat, unsigned int n_iter, unsigned int burn_in, unsigned int thin , Rcpp::String P0_prior_name){

	// Note that there is not the //' @export command. The user can not call this function.
	// I am afraid that Rcpp can take only Column major matrices. (not sure)
	// Do not use deafult values here
	Rcpp::Rcout<<"This is the Rcpp function"<<std::endl;
	Rcpp::Rcout<<"In c++ environment you can create custom c++ classes"<<std::endl;
  GibbsSampler Gibbs;
  Gibbs.n_iter=n_iter;
  Gibbs.burn_in=burn_in;
  Gibbs.thin=thin;
  GS_data gsData;

  for (unsigned int j = 0; j < dat.rows(); ++j) {
        std::vector<double> v;
        for (unsigned int i = 0; i <dat.cols() ; ++i) {
            if(!isnan(dat(j,i))){
             v.push_back(dat(j,i));
            }
        }
        gsData.data.push_back(v);
  }

  std::map<string, std::vector<double>> out=Gibbs.sample();
  std::vector<double> M=out["M"];
  std::vector<double> Mstar=out["M*"];
  std::vector<double> K=out["K"];


    // Parameters param(niter, burnin, thin);   // example of another class that stores useful options
	if (P0_prior_name == "Normal-InvGamma")
	{
		// Hyperparameters hy(param1,param2,param3); 			//example of the creation of a custom c++ class
		//P0Prior P0(hy.nu, ...);
			    //--> example of the creation of a custom c++ class that defines the P0prior to be used in this case

		//SamplingStrategy sampler_obj(data, param, hy, P0);   //example of sampler object creation
		//sampler_obj(...);
				//--> run the sampler. prefer to use the call operator instead of a method called run().


		//Post-processing and return

		//you can mix types in Rcpp lists
		return Rcpp::List::create ( Rcpp::Named("M")= M,
                                  	Rcpp::Named("Mstar")= Mstar,
                                  	Rcpp::Named("K")= K,
                                  	Rcpp::Named("Data")=dat);

	}
	else if (P0_prior_name == "Normal-Gamma")
	{
		//Similar as before
		return Rcpp::List::create ( Rcpp::Named("return_1")= dat,
		                            Rcpp::Named("return_2")= 10,
		                            Rcpp::Named("return_3")= "string" );
	}
	else
		throw std::runtime_error("Runtime error. This avoid R session to abort in favour of a meaningful error.");
}



#endif
