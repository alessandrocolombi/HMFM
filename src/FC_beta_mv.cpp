#include "FC_beta_mv.h"


void FC_beta_mv::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){

	//Retrive all data needed from gs_data
	const unsigned int& d = gs_data.d; // number of groups
	const unsigned int& r = gs_data.r; //number of covariates
	const std::vector<unsigned int>& n_j = gs_data.n_j; // number of observations per group
	const std::vector<double>& mu = gs_data.mu; // vector of the means for each component
	const std::vector<double>& sigma = gs_data.sigma; // vector of variances for each component
	const std::vector< std::vector<unsigned int>>& Ctilde = gs_data.Ctilde; // matrix of partition
	std::vector<std::vector<Individual>>& mv_data = gs_data.mv_data; //matrix of data we don't copy it since data can be big but we use a pointer
	GDFMM_Traits::MatRow& beta = gs_data.beta; // dxr matrix of regression coefficients

	sample::rmvnorm rmv; //Covariance parametrization

	for(size_t j = 0; j < d; j++){ // update regression coeff for j-th level

		// invSigma1 must contain a rxr matrix sum_i{ 1/sigma^2_ji * X_ji * X_ji^T }
		GDFMM_Traits::MatCol invSigma1( GDFMM_Traits::MatCol::Zero(r,r) ); 

		// beta1 must contain a r lenght vector sum_i{ 1/sigma^2_ji *  X_ji * (Y_ji - mu_ji 1_Ni)  }
		GDFMM_Traits::VecCol beta1( GDFMM_Traits::VecCol::Zero(r) );

		for(size_t i = 0; i < n_j[j]; i++){
			unsigned int C_ji = Ctilde[j][i]; //C_ji is the component mixture that defines mean-variance for obseration ji
			Individual& data_ji = mv_data[j][i]; // shortcut, just for notation

			if( data_ji.n_ji > 0 ){
				invSigma1 += 1.0/sigma[C_ji] * data_ji.X_jiX_jiT; //update invSigma1

				Eigen::Map<GDFMM_Traits::VecCol> eigen_data( &(data_ji.obs_ji[0]), data_ji.n_ji ); //cast observation into eigen form
				GDFMM_Traits::VecCol cl_means( GDFMM_Traits::VecCol::Constant(data_ji.n_ji, mu[C_ji]) ); // define a vector where each element is equal to mu_ji
				GDFMM_Traits::VecCol temp( eigen_data - cl_means ); // compute the difference
				beta1 += 1.0/sigma[C_ji] * data_ji.X_ji * temp; // update beta1
			}
			else
				throw std::runtime_error("Error in FC_beta_mv. è possibile che arrivi qui? Ho mv_data[j][i].n_ji == 0 ");
		}

		GDFMM_Traits::MatCol I( GDFMM_Traits::MatCol::Identity(r,r) );
		GDFMM_Traits::MatCol prec_post = invSigma0 + invSigma1 ; // compute the posterior precision matrix
		GDFMM_Traits::MatCol cov_post = prec_post.llt().solve(I) ; // invert the posterior precision matrix to obtain the posterior covariance matrix
		GDFMM_Traits::VecCol post_mean = cov_post * (invSigma0beta0 + beta1); // compute the posterior mean
		
		beta.row(j) = rmv(gs_engine, post_mean, cov_post); // sample and update beta_j

	}

	update_data_summary_statistics(gs_data);

}

void FC_beta_mv::update_data_summary_statistics(GS_data& gs_data){
	//Retrive all data needed from gs_data
	const unsigned int& d = gs_data.d; // number of groups
	const unsigned int& r = gs_data.r; //number of covariates
	const std::vector<unsigned int>& n_j = gs_data.n_j; // number of observations per group
	std::vector<std::vector<Individual>>& mv_data = gs_data.mv_data; //matrix of data we don't copy it since data can be big but we use a pointer
	GDFMM_Traits::MatRow& beta = gs_data.beta; // dxr matrix of regression coefficients

	for(size_t j = 0; j < d; j++){ 
		for(size_t i = 0; i < n_j[j]; i++){
			Individual& data_ji = mv_data[j][i]; // shortcut, just for notation
			if( data_ji.n_ji > 0 ){
				GDFMM_Traits::VecCol beta_s = beta.row(j);
				double beta_s_zji = data_ji.z_ji.dot( beta_s );
				data_ji.Ybar_star_ji = data_ji.mean_ji - 1.0/(double)data_ji.n_ji * beta_s_zji;
				if( data_ji.n_ji > 1 ){
					Eigen::Map<GDFMM_Traits::VecCol> eigen_data( &(data_ji.obs_ji[0]), data_ji.n_ji ); //cast observation into eigen form
					GDFMM_Traits::VecCol ones( GDFMM_Traits::VecCol::Constant(1.0,data_ji.n_ji) );
					GDFMM_Traits::VecCol first_vect( eigen_data - data_ji.mean_ji*ones );
					GDFMM_Traits::VecCol second_vect( data_ji.X_ji.transpose()*beta_s - 1.0/(double)data_ji.n_ji*beta_s_zji*ones );
					double double_product = first_vect.dot( second_vect ); 
					data_ji.Vstar_ji = 	data_ji.var_ji + 
										1.0/(double)(data_ji.n_ji - 1)*( beta_s.transpose()*data_ji.X_jiX_jiT*beta_s -  1.0/(double)data_ji.n_ji*beta_s_zji*beta_s_zji ) - 
										2.0/(double)(data_ji.n_ji - 1)*double_product;
				}
			}
			else
				throw std::runtime_error("Error in FC_beta_mv. è possibile che arrivi qui? Ho mv_data[j][i].n_ji == 0 ");
			
		}
	}
}