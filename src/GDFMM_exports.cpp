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
#include "GibbsSamplerMarginal.h"

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
	out_data out = Gibbs.out; //copia inutile
	auto n_j = Gibbs.get_nj();
	// A cosa servono queste cose? Non sono copie inutili?
	//std::vector<std::vector<double>> mu = out.mu; // old version - secondo me non serve nemmeno
	//std::vector<std::vector<double>> sigma = out.sigma; // old version - secondo me non serve nemmeno
	std::vector< GDFMM_Traits::MatRow > S = out.S; //copia inutile
    //std::vector<GDFMM_Traits::MatRow> w_jk = out.w_jk; //POSSO TOGLIERE COMPLETAMENTE IL CALCOLO DEI w_jk??
	std::vector<double> lambda = out.lambda; //copia inutile
	Rcpp::NumericMatrix gamma(dat.rows(), n_iter, out.gamma.begin());
	Rcpp::NumericMatrix U(dat.rows(), n_iter, out.U.begin());
	if(FixPart){
		
		return Rcpp::List::create( Rcpp::Named("K") = out.K,
									Rcpp::Named("Mstar") = out.Mstar,
									Rcpp::Named("mu") = out.mu,
                                  	Rcpp::Named("sigma") = out.sigma,
									Rcpp::Named("S") =  out.S,  //Rcpp::Named("w_jk") =  w_jk, //POSSO TOGLIERE COMPLETAMENTE IL CALCOLO DEI w_jk??
									Rcpp::Named("gamma") = gamma,
									Rcpp::Named("lambda") = lambda,
									Rcpp::Named("U") = U,
									Rcpp::Named("log_sum") = out.log_prod_psiU //il parametro della Poisson di Mstar è lambda*exp(-log_sum)
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
			
			fpmatr(it, Rcpp::_) = Rcpp::NumericMatrix( 1, n_data, fprowvec.begin() ); //fpmatr[it,:] in R notation
		}

		return Rcpp::List::create( Rcpp::Named("K") = K,
									Rcpp::Named("Mstar") = Mstar,
									Rcpp::Named("Partition") = fpmatr,
									Rcpp::Named("mu") = out.mu,
                                  	Rcpp::Named("sigma") = out.sigma,
									Rcpp::Named("gamma") = gamma,
									Rcpp::Named("lambda") = lambda,
									Rcpp::Named("U") = U,
									Rcpp::Named("S") =  S, //aggiunto
									Rcpp::Named("log_sum") = out.log_prod_psiU //il parametro della Poisson di Mstar è lambda*exp(-log_sum)
									);
	}
    
}

//' GDFMM sampler
// [[Rcpp::export]]
Rcpp::List GDFMM_marginal_sampler_c( Eigen::MatrixXd const & dat, unsigned int n_iter, unsigned int burn_in,
			 						 unsigned int thin , unsigned int seed, Rcpp::String P0_prior_name, 
									 bool FixPart, Rcpp::List option)
{


    // Create object GibbsSamplerMarginal and sample
	GibbsSamplerMarginal Gibbs(dat, n_iter, burn_in, thin, seed, P0_prior_name, FixPart, option);
    Gibbs.sample();

	// Take output data from the sample
	return Rcpp::List::create( 	Rcpp::Named("K") = Gibbs.out.K,
								Rcpp::Named("mu") = Gibbs.out.mu,
	                            Rcpp::Named("sigma") = Gibbs.out.sigma,
								Rcpp::Named("gamma") = Gibbs.out.gamma,
								Rcpp::Named("lambda") = Gibbs.out.lambda,
								Rcpp::Named("U") = Gibbs.out.U,
								Rcpp::Named("Partition") = Gibbs.out.Partition
							);
    
}


//' Test
//' @export
// [[Rcpp::export]]
void Test_Rcpp(){

	std::vector<double> vec_c(5);
	std::iota(vec_c.begin(), vec_c.end(), 1.0);
	Rcpp::Rcout<<"Stampo vec_c: ";		
	for(auto __v : vec_c)
		Rcpp::Rcout<<__v<<", ";
	Rcpp::Rcout<<std::endl;

	Rcpp::NumericVector vec_rcpp1(vec_c.begin(),vec_c.end());
	Rcpp::Rcout<<"Stampo vec_rcpp1: ";		
	for(auto __v : vec_rcpp1)
		Rcpp::Rcout<<__v<<", ";
	Rcpp::Rcout<<std::endl;

	Rcpp::Rcout<<"Questo non va (uso std::copy)"<<std::endl;
	Rcpp::NumericVector vec_rcpp2;
	std::copy(vec_c.begin(), vec_c.end(), vec_rcpp2.begin()); 
	Rcpp::Rcout<<"Stampo vec_rcpp2: ";		
	for(auto __v : vec_rcpp2)
		Rcpp::Rcout<<__v<<", ";
	Rcpp::Rcout<<std::endl;

	Rcpp::NumericVector vec_rcpp3;
	Rcpp::NumericVector vec_rcpp_temp(vec_c.begin(),vec_c.end());
	Rcpp::Rcout<<"Stampo vec_rcpp_temp: ";		
	for(auto __v : vec_rcpp_temp)
		Rcpp::Rcout<<__v<<", ";
	Rcpp::Rcout<<std::endl;
	std::swap(vec_rcpp3, vec_rcpp_temp);
	//vec_rcpp3.reserve(5);
	Rcpp::Rcout<<"Stampo vec_rcpp3: ";		
	for(auto __v : vec_rcpp3)
		Rcpp::Rcout<<__v<<", ";
	Rcpp::Rcout<<std::endl;

	Rcpp::Rcout<<"Stampo vec_rcpp_temp: ";		
	for(auto __v : vec_rcpp_temp)
		Rcpp::Rcout<<__v<<", ";
	Rcpp::Rcout<<std::endl;
	
	
} 

#endif
