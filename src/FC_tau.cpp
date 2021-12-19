
#include "FC_tau.h"



void FC_tau::update(GS_data& gs_data, sample::GSL_RNG gs_engine) {
    unsigned int Mna = gs_data.Mstar; // number of unallocated components
    unsigned int M = gs_data.M;
    unsigned int d = gs_data.d; // number of groups
    unsigned int K=gs_data.K; //number of clusters
    std::vector<unsigned int> n_j= gs_data.n_j; // number of observations per group
    std::vector< std::vector<unsigned int>> C=p.C; // C matrix of partition
    std::vector<unsigned int> clust_out= p.clust_out; // Vector of clusters
    GDFMM_Traits::MatUnsCol N=gs_data.N; // Matrix of observation oper cluster per group
    std::vector<unsigned int> ind_i; // i index of C elements
    std::vector<unsigned int> ind_j;// j index of C elements
    std::vector< std::vector<unsigned int>> Ctilde; // matrix of partition
    std::vector<std::vector<double>> *data=gs_data.data; //matrix of data we don't copy it since data can be big but we use a pointer
    std::string c=gs_data.prior;
    if (c == 'normal-inverse-gamma') {


        sample::rgamma Gamma;
        sample::rnorm rnorm;

        //tau-nonallocate

        for (int m = K; m < M; ++m){
             float sigma2_na=1 / Gamma(gs_engine, nu_0/2, (nu_0)*((sigma_0)/2)); // Non allocated Components' variance
             float mu_na=rnorm( mu_0, std::sqrt(sigma2_na / k_0)); // Non allocated Components' mean
             gs_data.mu[m]=mu_na;
             gs_data.sigma[m]=sigma2_na;
        }


         //tau-allocate

        for (int m = 0; m < K; ++m) {
            for (int j = 0; j <d ; ++j) {
                for (int i = 0; i <n_j[j] ; ++i) {
                    if(C[j][i]==clust_out[m]){
                        N(j,m)++;
                        ind_i.push_back(i);
                        ind_j.push_back(j);
                        Ctilde[j][m]=m;
                    }
                }
//set Ctilde in partition
            double nu_n_clust = nu_0 + N_k[m]; //manca da definire N_k[]
            double lpk = k_0 + N_k[m];
            double y_bar_clust= mean(ind_i, ind_j);
            double s2_clust= var(y_bar_clust, ind_i, ind_j);
//if (is.na(s2_clust[k])){s2_clust[k] <- 0}
            mu_n_clust = (k_0 * mu_0 + N_k[m] * y_bar_clust) / lpk;
            sigma2_n_clust = (nu_0 * (sigma_0 * sigma_0) + (N_k[m] - 1) * s2_clust+ k_0 * N_k[m] * (y_bar_clust - mu_0) * (y_bar_clust - mu_0) / (lpk));

            //Campionamento
            sigma2_a = 1 / rgamma(gs_engine, nu_n_clust/ 2, sigma2_n_clust / 2);
            mu_a=rnorm(gs_engine, mu_n_clust, sqrt(sigma2_a / lpk));
            gs_data.mu[m]=mu_a;
            gs_data.sigma[m]=sigma2_a;
        }



}
}

double mean(std::vector<unsigned int> ind_i, std::vector<unsigned int> ind_j){
    int count=0;
    double sum=0;
    for (int m = 0; m <ind_i.size(); ++m) {
        sum+=data->at(ind_j[m])->at(ind_i[m]);
        count++;

    }
    return sum/count;
}
 double var(double mean,std::vector<unsigned int> ind_i,std::vector<unsigned int>ind_j){
    double var;
    int count=0;
     for (int m = 0; m <ind_i.size() ; ++m) {
         var+=(data->at(ind_j[m])->at(ind_i[m]) - mean) * (data->at(ind_j[m])->at()ind_i[m]) - mean)) ;
         count++;
     }
     return var/count;
 }


// mean e var vanno bene qui e poi si possono chiamare tranquillamente da dentro senza definire il namespace?
