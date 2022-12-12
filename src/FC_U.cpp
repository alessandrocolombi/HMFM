//
// Created by dario on 15/12/2021.
//
#include <Rcpp.h>
#include <RcppEigen.h>

#include "FC_U.h"

void FC_U::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) {

    const unsigned int& d = gs_data.d;
    //const double& nu = gs_data.nu;
    const std::vector<unsigned int>& n_j = gs_data.n_j;
    const GDFMM_Traits::MatRow& S = gs_data.S;
    // T_j is computed for each group (T_j = sum of S_jm over m for each group j)

    Eigen::VectorXd T = S.rowwise().sum(); 
    // Sampler for new U
    sample::rgamma Gamma;
    //Rcpp::Rcout<<T(0)<<std::endl;

    // Rcpp::Rcout << "T = ";
    for (unsigned j=0; j<d; j++) { // for loop per livelli
        gs_data.U[j]= Gamma(gs_engine, n_j[j], 1.0/T(j)); //This is U'
        /*
        Old Order
        gs_data.U[j] *= nu; // compute U', U'|rest is gamma(n_j, T'_j), U' = nu * U
        */
        // Rcpp::Rcout << T(j) << " ";
    }
    // Rcpp::Rcout << std::endl;
    gs_data.update_log_sum();
    // Rcpp::Rcout<< "New log_sum : " << gs_data.log_sum <<std::endl;
    // Rcpp::Rcout<< "First value of U : " << gs_data.U[0];
}
