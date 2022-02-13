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
    std::vector<GDFMM_Traits::MatRow> w_ji = out.w_ji;
	std::vector<double> lambda = out.lambda;
	Rcpp::NumericMatrix gamma(dat.rows(), n_iter, out.gamma.begin());
	Rcpp::NumericMatrix U(dat.rows(), n_iter, out.U.begin());

	if(FixPart){
		
		return Rcpp::List::create( Rcpp::Named("mu") = mu,
                                  	Rcpp::Named("sigma") = sigma,
									Rcpp::Named("w_ji") =  w_ji,
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


#endif
