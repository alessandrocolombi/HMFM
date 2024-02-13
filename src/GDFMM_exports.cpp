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
#include "ConditionalSampler.h"
//' @importFrom RcppParallel RcppParallelLibs

# include "Individual.h"

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
								Rcpp::Named("Partition") = Gibbs.out.Partition,
								Rcpp::Named("log_q") = Gibbs.out.log_q
							);

}


//' Test new data
void Test_data(const Rcpp::List& data_list ){

	unsigned int n = Rcpp::as<unsigned int>(data_list["n"]);
	unsigned int d = Rcpp::as<unsigned int>(data_list["d"]);
	std::vector<unsigned int> n_j{ Rcpp::as<std::vector<unsigned int>>(data_list["n_j"]) };
	std::vector<std::string> ID_i{ Rcpp::as<std::vector<std::string>>(data_list["ID_i"]) };
	std::vector<unsigned int> s_i{ Rcpp::as<std::vector<unsigned int>>(data_list["s_i"]) };
	Rcpp::IntegerMatrix N_ji    = Rcpp::as<Rcpp::IntegerMatrix>(data_list["N_ji"]);
	Rcpp::NumericMatrix mean_ji = Rcpp::as<Rcpp::NumericMatrix>(data_list["mean_ji"]);
	Rcpp::NumericMatrix var_ji  = Rcpp::as<Rcpp::NumericMatrix>(data_list["var_ji"]);
	Rcpp::List obs = Rcpp::as<Rcpp::List>(data_list["observations"]);


	/*
	Rcpp::Rcout<<"n = "<<n<<std::endl;
	Rcpp::Rcout<<"d = "<<d<<std::endl;
	Rcpp::Rcout<<"Stampo n_j: ";
	for(auto __v : n_j)
		Rcpp::Rcout<<__v<<", ";
	Rcpp::Rcout<<std::endl;
	Rcpp::Rcout<<"Stampo ID_i: ";
	for(auto __v : ID_i)
		Rcpp::Rcout<<__v<<", ";
	Rcpp::Rcout<<std::endl;
	Rcpp::Rcout<<"Stampo s_i: ";
	for(auto __v : s_i)
		Rcpp::Rcout<<__v<<", ";
	Rcpp::Rcout<<std::endl;

	Rcpp::Rcout<<"N_ji is a "<< N_ji.nrow() <<" x "<< N_ji.ncol()<<" matrix "<<std::endl;
	Rcpp::Rcout<<"mean_ji is a "<< mean_ji.nrow() <<" x "<< mean_ji.ncol()<<" matrix "<<std::endl;
	Rcpp::Rcout<<"var_ji is a "<< var_ji.nrow() <<" x "<< var_ji.ncol()<<" matrix "<<std::endl;


	Rcpp::Rcout<<"Stampo N_ji: "<<std::endl;
	for (int i = 0; i < N_ji.nrow(); i++) {
  		for (int j = 0; j < N_ji.ncol(); j++) {
    		Rcpp::Rcout << N_ji(i, j) << " ";
  		}
  		Rcpp::Rcout << std::endl;
	}
	Rcpp::Rcout<<"Stampo mean_ji: "<<std::endl;
	for (int i = 0; i < mean_ji.nrow(); ++i) {
  		for (int j = 0; j < mean_ji.ncol(); ++j) {
    		Rcpp::Rcout << mean_ji(i, j) << " ";
  		}
  		Rcpp::Rcout << std::endl;
	}
	Rcpp::Rcout<<"Stampo var_ji: "<<std::endl;
	for (int i = 0; i < var_ji.nrow(); ++i) {
  		for (int j = 0; j < var_ji.ncol(); ++j) {
    		Rcpp::Rcout << var_ji(i, j) << " ";
  		}
  		Rcpp::Rcout << std::endl;
	}
	*/

	std::vector<std::vector<Individual>> data;
	data.resize(d);
	for(size_t j = 0; j < d; j++)
		data[j].reserve(n_j[j]);


	//Individual Adamo( ID_i[0], (unsigned int)N_ji(0,0), mean_ji(0,0), var_ji(0,0) );
	//Rcpp::Rcout<<"Adamo name = "<<Adamo.ID<<std::endl;
	//Rcpp::Rcout<<"Adamo n_ji = "<<Adamo.n_ji<<std::endl;
	//Rcpp::Rcout<<"Adamo mean_ji = "<<Adamo.mean_ji<<std::endl;
	//Rcpp::Rcout<<"Adamo var_ji = "<<Adamo.var_ji<<std::endl;

	Rcpp::Rcout<<"data.size() = "<<data.size()<<std::endl;
	for(size_t j = 0; j < d; j++){
		Rcpp::List obs_j = obs[j];
		for(size_t i = 0; i < n; i++){
			//Rcpp::Rcout<<j<<" , "<<i<<std::endl;
			if(N_ji(j,i) > 0){

				std::vector<double> obs_ji = Rcpp::as<std::vector<double>>(obs_j[i]);
				Rcpp::Rcout<<"Stampo obs_"<<j<<i<<std::endl;

				Rcpp::Rcout<<"Stampo obs_ji: ";
				for(auto __v : obs_ji)
				    Rcpp::Rcout<<__v<<", ";
				Rcpp::Rcout<<std::endl;


				Individual data_ji( ID_i[i], (unsigned int)N_ji(j,i), mean_ji(j,i), var_ji(j,i), obs_ji );
				data[j].push_back(data_ji);
			}
		}
		Rcpp::Rcout<<"data["<<j<<"].size() = "<<data[j].size()<<std::endl;
	}
	//Rcpp::Rcout<<"data[3][0] name = "<<data[3][0].ID<<std::endl;
	//Rcpp::Rcout<<"data[3][0] n_ji = "<<data[3][0].n_ji<<std::endl;
	//Rcpp::Rcout<<"data[3][0] mean_ji = "<<data[3][0].mean_ji<<std::endl;
	//Rcpp::Rcout<<"data[3][0] var_ji = "<<data[3][0].var_ji<<std::endl;
	//Rcpp::Rcout<<"--------------------------"<<std::endl;
	//Rcpp::Rcout<<"data[3][1] name = "<<data[3][1].ID<<std::endl;
	//Rcpp::Rcout<<"data[3][1] n_ji = "<<data[3][1].n_ji<<std::endl;
	//Rcpp::Rcout<<"data[3][1] mean_ji = "<<data[3][1].mean_ji<<std::endl;
	//Rcpp::Rcout<<"data[3][1] var_ji = "<<data[3][1].var_ji<<std::endl;
}

//' MCMC_conditional_c
// [[Rcpp::export]]
Rcpp::List MCMC_conditional_c( const Rcpp::List& data_list,
							   unsigned int n_iter, unsigned int burn_in, unsigned int thin ,
							   unsigned int seed, Rcpp::String P0_prior_name,
							   bool FixPart, Rcpp::String algorithm, Rcpp::List option){
	//Rcpp::Rcout<<"Dentro conditional sampler"<<std::endl;
    // Create object ConditionalSampler and sample
	ConditionalSampler Gibbs(data_list, n_iter, burn_in, thin, seed, P0_prior_name, FixPart, algorithm, option);
    Gibbs.sample();

	// Take output data from the sample
	Rcpp::NumericMatrix gamma( Gibbs.d, n_iter, Gibbs.out.gamma.begin() );
	Rcpp::NumericMatrix U( Gibbs.d, n_iter, Gibbs.out.U.begin() );

	if(FixPart){

		return Rcpp::List::create( Rcpp::Named("K") = Gibbs.out.K,
									Rcpp::Named("Mstar") = Gibbs.out.Mstar,
									Rcpp::Named("mu") = Gibbs.out.mu,
                                  	Rcpp::Named("sigma") = Gibbs.out.sigma,
									Rcpp::Named("S") =  Gibbs.out.S,
									Rcpp::Named("beta") =  Gibbs.out.beta,
									Rcpp::Named("gamma") = gamma,
									Rcpp::Named("lambda") = Gibbs.out.lambda,
									Rcpp::Named("U") = U,
									Rcpp::Named("log_sum") = Gibbs.out.log_prod_psiU //il parametro della Poisson di Mstar è lambda*exp(-log_sum)
									);
	}
	else{

		const std::vector<std::vector< std::vector<unsigned int>>>& C = Gibbs.out.Ctilde;

			//we need a better structure for C--> we save it as a vector instead of a matrix
		std::vector<unsigned int> fprowvec; //final partition rowvec

		// store total number of data (i.e. cardinality{y_ji})
		unsigned int n_data = std::accumulate(Gibbs.n_j.begin(), Gibbs.n_j.end(), 0);
		Rcpp::NumericMatrix fpmatr(n_iter, n_data );  //final partition matrix

		for (unsigned it = 0; it < n_iter; it++){
			fprowvec.clear();

			for(unsigned j=0; j<Gibbs.d; j++){
				fprowvec.insert(fprowvec.end(), C[it][j].begin(), C[it][j].end());
			}

			fpmatr(it, Rcpp::_) = Rcpp::NumericMatrix( 1, n_data, fprowvec.begin() ); //fpmatr[it,:] in R notation
		}

		return Rcpp::List::create( Rcpp::Named("K") = Gibbs.out.K,
									Rcpp::Named("Mstar") = Gibbs.out.Mstar,
									Rcpp::Named("Partition") = fpmatr,
									Rcpp::Named("mu") = Gibbs.out.mu,
                                  	Rcpp::Named("sigma") = Gibbs.out.sigma,
									Rcpp::Named("gamma") = gamma,
									Rcpp::Named("lambda") = Gibbs.out.lambda,
									Rcpp::Named("U") = U,
									Rcpp::Named("S") =  Gibbs.out.S,
									Rcpp::Named("beta") =  Gibbs.out.beta,
									Rcpp::Named("log_sum") = Gibbs.out.log_prod_psiU //il parametro della Poisson di Mstar è lambda*exp(-log_sum)
									);
	}

}



//' Test
void Test_Rcpp(const Rcpp::NumericMatrix& X){

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

	Rcpp::Rcout<<"test matrice"<<std::endl;
    Eigen::Map<Eigen::MatrixXd> X_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X));
    Rcpp::Rcout<<"X_eig:"<<std::endl<<X_eig<<std::endl;

}

#endif
