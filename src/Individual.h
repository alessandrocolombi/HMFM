#ifndef GDFMM_INDIVIDUAL_H
#define GDFMM_INDIVIDUAL_H

// [[Rcpp::depends(RcppEigen)]]
#include "include_headers.h"
#include "recurrent_traits.h"
#include <Rcpp.h>
#include <RcppEigen.h>

class Individual{
public :
    Individual(const std::string& _id, unsigned int _n, double _mean, double _var, const std::vector<double>& _obs_ji ) : ID(_id), n_ji(_n), mean_ji(_mean), var_ji(_var), obs_ji(_obs_ji){
        Ybar_star_ji = mean_ji;
        Vstar_ji    = var_ji;
        z_ji = GDFMM_Traits::VecCol::Zero(1);
        //Rcpp::Rcout<<"Caso senza covariate, avevo un bug qua. è stato sistemato?"<<std::endl;
        
    };
    Individual(const std::string& _id, unsigned int _n, double _mean, double _var, 
                const std::vector<double>& _obs_ji, const Rcpp::NumericMatrix& _X_ji ) : ID(_id), n_ji(_n), mean_ji(_mean), var_ji(_var), obs_ji(_obs_ji){
        
        Eigen::Map<Eigen::MatrixXd> X_eig(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>( _X_ji ));
        int num_obs = X_eig.rows();
        int num_cov = X_eig.cols();
        X_ji = GDFMM_Traits::MatCol::Zero(num_cov, num_obs);
        X_ji = X_eig.transpose();
                //Rcpp::Rcout<<"X_ji:"<<std::endl<<X_ji<<std::endl;
        X_jiX_jiT = X_ji*X_ji.transpose();
        z_ji = X_ji.rowwise().sum();

                //Rcpp::Rcout<<"X_jiX_jiT:"<<std::endl<<X_jiX_jiT<<std::endl;
                //Rcpp::Rcout<<"z_ji:"<<std::endl<<z_ji<<std::endl;
        Ybar_star_ji = mean_ji;
        Vstar_ji    = var_ji;
                //Rcpp::Rcout<<"Ybar_star_ji:"<<std::endl<<Ybar_star_ji<<std::endl;
                //Rcpp::Rcout<<"Vstar_ji:"<<std::endl<<Vstar_ji<<std::endl;

    };
    Individual() = default;

    std::string ID; // individual ID
    unsigned int n_ji; // number of observations taken for individual i in level j
    double mean_ji; // mean of observations taken for individual i in level j
    double var_ji; // variance of observations taken for individual i in level j
    std::vector<double> obs_ji; // vector containing all observations for individual i in level j
    GDFMM_Traits::MatCol X_ji; // (r x n_ji) matrix containing the design matrix for individual i in levele j. r is the number of covariates
    GDFMM_Traits::MatCol X_jiX_jiT; //(r x r) matrix that contains X_ji %*% t(X_ji)
    GDFMM_Traits::VecCol z_ji; // vector of length r that contains X_ji %*% ones(n_ji)
    double Vstar_ji;
    double Ybar_star_ji;

    virtual ~Individual() = default;
};
#endif


/*
Idea: definisco per ogni individuo altre due variabili, V_ji e Ybar_ji. Le inizializzo uguali a media e varianza. 
Poi, nel caso senza covariate le lascio sempre uguali, non le aggiorno mai. Nel caso con le covariate, le aggiorno
dopo aver aggiornato i beta. 
In Partition_mv::log_dmvnorm, invece di usare mean_ji e var_ji uso queste nuove variabili che nel caso con covariate 
rappresentano le statistiche sufficienti delle osservazioni Ystar_ji = Y_ji - X^T_ji*beta_j.
Per aggiornare le componenti della mistura, devo solo aggiornate FC_tau_mv::compute_cluster_summaries(), anche qua devo
usare le nuove variabili e devo aggiornare anche il modo in cui viene calcolata W
*/