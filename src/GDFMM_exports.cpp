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
#include "ComponentPrior_factory.h"

//' GDFMM sampler
// [[Rcpp::export]]
Rcpp::List GDFMM_sampler_c( Eigen::MatrixXd const & dat, unsigned int n_iter, unsigned int burn_in,
			 					unsigned int thin , unsigned int seed, Rcpp::String P0_prior_name, 
								bool FixPart, Rcpp::List option){

	// Note that there is not the //' @export command. The user can not call this function.
	// I am afraid that Rcpp can take only Column major matrices. (not sure)
	// Do not use deafult values here
	// Rcpp::Rcout<<"This is the Rcpp function"<<std::endl;
	// Rcpp::Rcout<<"In c++ environment you can create custom c++ classes"<<std::endl;
    // Create object GibbsSampler and sample
	GibbsSampler Gibbs(dat, n_iter, burn_in, thin, seed, P0_prior_name, FixPart, option);
    Gibbs.sample();
	// Take output data from the sample
	out_data out = Gibbs.out;
	auto n_j = Gibbs.get_nj();

	std::vector<std::vector<double>> mu = out.mu;
	std::vector<std::vector<double>> sigma = out.sigma;
    std::vector<GDFMM_Traits::MatRow> w_jk = out.w_jk;
	std::vector<double> lambda = out.lambda;
	Rcpp::NumericMatrix gamma(dat.rows(), n_iter, out.gamma.begin());
	Rcpp::NumericMatrix U(dat.rows(), n_iter, out.U.begin());

	if(FixPart){
		
		return Rcpp::List::create( Rcpp::Named("mu") = mu,
                                  	Rcpp::Named("sigma") = sigma,
									Rcpp::Named("w_jk") =  w_jk,
									Rcpp::Named("gamma") = gamma,
									Rcpp::Named("lambda") = lambda,
									Rcpp::Named("U") = U
									);
	}
	else{
		std::vector<unsigned int> K = out.K;
		std::vector<unsigned int> Mstar = out.Mstar;
		std::vector<std::vector< std::vector<unsigned int>>> C = out.Ctilde;

			//we need a better structure for C--> we save it as a vector instead of a matrix
		std::vector<unsigned int> fprowvec; //final partition rowvec
		
		// store total number of data (i.e. cardinality{y_ji})
		unsigned int n_data = std::accumulate(n_j.begin(), n_j.end(), 0);
		Rcpp::NumericMatrix fpmatr(n_iter, n_data );  //final partition matrix

		for (unsigned it = 0; it < n_iter; it++){
			fprowvec.clear();
			
			for(unsigned j=0; j<dat.rows(); j++){
				fprowvec.insert(fprowvec.end(), C[it][j].begin(), C[it][j].end());
			}
			
			fpmatr(it, Rcpp::_) = Rcpp::NumericMatrix( 1, n_data, fprowvec.begin() );
		}

		return Rcpp::List::create( Rcpp::Named("K") = K,
									Rcpp::Named("Mstar") = Mstar,
									Rcpp::Named("Partition") = fpmatr,
									Rcpp::Named("mu") = mu,
                                  	Rcpp::Named("sigma") = sigma,
									Rcpp::Named("gamma") = gamma,
									Rcpp::Named("lambda") = lambda,
									Rcpp::Named("U") = U
									);
	}
    
}


//' Test ComponentPrior
//' @export
// [[Rcpp::export]]
void Test_Prior(){
	Rcpp::String form = "Poisson";
	double lambda(5.0);
	ComponentPrior_Parameters param;
	param.Lambda = lambda;
	auto qM_ptr = Select_ComponentPrior(form, param);
	ComponentPrior& qM(*qM_ptr);
	std::string i_am = qM_ptr->showMe();
	std::string i_am2 = qM.showMe();
	Rcpp::Rcout<<"I am: "<<i_am<<std::endl;
	Rcpp::Rcout<<"I am - II - : "<<i_am2<<std::endl;
	form = "NegativeBinomial";
	param.p = 0.45;
	param.n_succ = 2.0;
	auto qM_ptr2 = Select_ComponentPrior(form, param);
	i_am = qM_ptr2->showMe();
	Rcpp::Rcout<<"I am: "<<i_am<<std::endl;
	// Devo testare se le densità sono giuste!
}

#endif
