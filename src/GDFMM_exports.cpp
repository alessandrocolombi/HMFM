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

Rcpp::List GDFMM_sampler_c( Eigen::MatrixXd const & dat, unsigned int n_iter, unsigned int burn_in,
			 unsigned int thin , unsigned int seed, Rcpp::String P0_prior_name, Rcpp::List option){

	// Note that there is not the //' @export command. The user can not call this function.
	// I am afraid that Rcpp can take only Column major matrices. (not sure)
	// Do not use deafult values here
	Rcpp::Rcout<<"This is the Rcpp function"<<std::endl;
	Rcpp::Rcout<<"In c++ environment you can create custom c++ classes"<<std::endl;
	  //GS_data g(dat, n_iter,burn_in,thin); non posso più inizializzare gsdata fuori
    GibbsSampler Gibbs(dat, n_iter, burn_in, thin, seed, option);
    //out_data out=Gibbs.sample();
    unsigned int N=dat.rows()*dat.cols();
    Rcpp::Rcout<<"N:"<<N<<std::endl;
    std::vector<int> Mstar=Gibbs.out.Mstar;
    std::vector<int> K=Gibbs.out.K;
    std::vector<double> lambda=Gibbs.out.lambda;
    std::vector<std::vector<double>> U=Gibbs.out.U;
    std::vector<std::vector<double>> gamma=Gibbs.out.gamma;
    std::vector<std::vector<std::vector<unsigned int>>> C=Gibbs.out.Ctilde;
    //we need a better structure for C
    std::vector<unsigned int> fprowvec;//final partition rowvec
    std::vector<std::vector<unsigned int>> fpmatr;  //final partition matrix

    for (unsigned k=0;k<n_iter-1;k++){
      fprowvec.clear();
      for(unsigned i=0; i<dat.rows(); i++){
        fprowvec.insert(fprowvec.end(), C[k][i].begin(), C[k][i].end());
      }
        fpmatr.push_back(fprowvec);

    }




    std::vector<std::vector<double>> mu=Gibbs.out.mu;
    std::vector<std::vector<double>> sigma=Gibbs.out.sigma;
    std::vector<GDFMM_Traits::MatRow> S=Gibbs.out.S;


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
   // Rcpp::Named("C")=C;
		//you can mix types in Rcpp lists
		return Rcpp::List::create ( Rcpp::Named("Mstar")= Mstar,
                                  	Rcpp::Named("K")= K,
                                  	Rcpp::Named("lambda")=lambda,
                                    Rcpp::Named("C")=C,
                                    Rcpp::Named("mu")=mu,
                                    Rcpp::Named("sigma")=sigma,
                                    Rcpp::Named("S")=S,
                                    Rcpp::Named("cout")=fpmatr
                                  	//Rcpp::Named("U")= U
                                  	//Rcpp::Named("gamma")=gamma
                                  	);

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
