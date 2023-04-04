#include "FullConditional.h"
bool FullConditional::binary_decision(double p1, const sample::GSL_RNG& engine) const{
    // Sample u from a Uniform(0,1) and verify if is less than p
    sample::runif unif;
    double u = unif(engine);
    return u<p1;
}

void FullConditional::print() const{
  std::cout<<name<<"\n";
}

double FullConditional::log_raising_factorial(const unsigned int& n, const double& a)const
{

  if(n==0)
    return 0.0;
  if(a<0)
    throw std::runtime_error("Error in my_log_raising_factorial, can not compute the raising factorial of a negative number in log scale"); 
  else if(a==0.0){
    return -std::numeric_limits<double>::infinity();
  }
  else{
    
    double val_max{std::log(a+n-1)};
    double res{1.0};
    if (n==1)
      return val_max;
    for(std::size_t i = 0; i <= n-2; ++i){
      res += std::log(a + (double)i) / val_max;
    }
    return val_max*res;

  }
}




std::tuple<double,double,double,unsigned int>
FullConditional::compute_cluster_summaries( const GDFMM_Traits::uset_uiui& current_indicies,
                                       const std::vector<std::vector<Individual>>& data,
                                       const GDFMM_Traits::MatRow& beta )const
{
    // initialize return objects
    double W{0.0}; // W = 1/(sum(pi)) * sum_{i=1}^{N_m}(pi * Xbari)
    double mean_of_vars{0.0}; // sum_{i=1}^{N_m}( (pi-1)*Vi )
    double sum_piX2{0.0}; // sum_{i=1}^{N_m}( pi*Xbari^2 )
    unsigned int sum_pi{0}; // sum(pi)

    // for each element in the cluster
    std::for_each(current_indicies.cbegin(), current_indicies.cend(), 
        [&W, &mean_of_vars, &sum_piX2, &sum_pi, &data, &beta ](const std::pair<unsigned int, unsigned int>& ind){
            const Individual& data_ji = data.at(ind.first).at(ind.second);
            sum_pi += data_ji.n_ji;
            mean_of_vars += (double)(data_ji.n_ji - 1)*data_ji.Vstar_ji;
            sum_piX2 += data_ji.n_ji * data_ji.Ybar_star_ji * data_ji.Ybar_star_ji;
            W += data_ji.n_ji*data_ji.mean_ji - data_ji.z_ji.dot( beta.row(ind.first) );
        }
        );
    W = W/(double)sum_pi;

    return( std::tie(W,mean_of_vars,sum_piX2,sum_pi) );
}

std::tuple<double,double,double,unsigned int>
FullConditional::compute_cluster_summaries( const GDFMM_Traits::uset_uiui& current_indicies,
                                            const std::vector<std::vector<Individual>>& data )const
{
    // initialize return objects
    double W{0.0}; // W = 1/(sum(pi)) * sum_{i=1}^{N_m}(pi * Xbari)
    double mean_of_vars{0.0}; // sum_{i=1}^{N_m}( (pi-1)*Vi )
    double sum_piX2{0.0}; // sum_{i=1}^{N_m}( pi*Xbari^2 )
    unsigned int sum_pi{0}; // sum(pi)

    // for each element in the cluster
    std::for_each(current_indicies.cbegin(), current_indicies.cend(), 
        [&W, &mean_of_vars, &sum_piX2, &sum_pi, &data ](const std::pair<unsigned int, unsigned int>& ind){
            const Individual& data_ji = data.at(ind.first).at(ind.second);
            sum_pi += data_ji.n_ji;
            mean_of_vars += (double)(data_ji.n_ji - 1)*data_ji.var_ji;
            sum_piX2 += data_ji.n_ji * data_ji.mean_ji * data_ji.mean_ji;
            W += data_ji.n_ji*data_ji.mean_ji;
        }
        );
    W = W/(double)sum_pi;

    return( std::tie(W,mean_of_vars,sum_piX2,sum_pi) );
}



/*Old implementation
std::tuple<double,double,double,unsigned int>
FullConditional::compute_cluster_summaries(  const std::vector<unsigned int>& ind_i,
                                       const std::vector<unsigned int>& ind_j,
                                       const std::vector<std::vector<Individual>>& data,
                                       const GDFMM_Traits::MatRow& beta )const
{
    // get number of elements in the cluster
    //unsigned int N_m{ind_i.size()};
    //Rcpp::Rcout<<"Dentro FullConditional::compute_cluster_summaries.cpp - con COVARIATE"<<std::endl;
    // initialize return objects
    double W{0.0}; // W = 1/(sum(pi)) * sum_{i=1}^{N_m}(pi * Xbari)
    double mean_of_vars{0.0}; // sum_{i=1}^{N_m}( (pi-1)*Vi )
    double sum_piX2{0.0}; // sum_{i=1}^{N_m}( pi*Xbari^2 )
    unsigned int sum_pi{0}; // sum(pi)

    // for each element in the cluster

    //  questa funzione non va bene se NON ho le covariate. il problema è l'ultimo termine, quello con beta e z.
    //  devo differenziare il calcolo di z in base a r=0 e r>0
    for(size_t ii = 0; ii <ind_i.size() ; ii++){

        const Individual& data_ji = data.at(ind_j[ii]).at(ind_i[ii]);

        //Rcpp::Rcout<<"data_ji.n_ji:"<<std::endl<<data_ji.n_ji<<std::endl;
        //Rcpp::Rcout<<"data_ji.Ybar_star_ji:"<<std::endl<<data_ji.Ybar_star_ji<<std::endl;
        //Rcpp::Rcout<<"data_ji.Vstar_ji:"<<std::endl<<data_ji.Vstar_ji<<std::endl;
        //Rcpp::Rcout<<"data_ji.mean_ji:"<<std::endl<<data_ji.mean_ji<<std::endl;
        //Rcpp::Rcout<<"data_ji.z_ji:"<<std::endl<<data_ji.z_ji<<std::endl;

        sum_pi += data_ji.n_ji;
        mean_of_vars += (double)(data_ji.n_ji - 1)*data_ji.Vstar_ji;
        sum_piX2 += data_ji.n_ji * data_ji.Ybar_star_ji * data_ji.Ybar_star_ji;
        //Rcpp::Rcout<<"Voglio uno scalare"<<std::endl;
        //Rcpp::Rcout<<"data_ji.z_ji:"<<std::endl<<data_ji.z_ji<<std::endl;
        //Rcpp::Rcout<<"beta.row(ind_j[ii]):"<<std::endl<<beta.row(ind_j[ii])<<std::endl;
        //Rcpp::Rcout<<"data_ji.z_ji.dot( beta.row(ind_j[ii]) ):"<<std::endl<<data_ji.z_ji.dot( beta.row(ind_j[ii]) )<<std::endl;
        W += data_ji.n_ji*data_ji.mean_ji - data_ji.z_ji.dot( beta.row(ind_j[ii]) );

    }
    W = W/(double)sum_pi;
    return( std::tie(_W,_mean_of_vars,_sum_piX2,_sum_pi) );
}

std::tuple<double,double,double,unsigned int>
FullConditional::compute_cluster_summaries(  const std::vector<unsigned int>& ind_i,
                                       const std::vector<unsigned int>& ind_j,
                                       const std::vector<std::vector<Individual>>& data )const
{
    // get number of elements in the cluster
    //unsigned int N_m{ind_i.size()};
    //Rcpp::Rcout<<"Dentro FullConditional::compute_cluster_summaries.cpp - SENZA convariate"<<std::endl;
    // initialize return objects
    double W{0.0}; // W = 1/(sum(pi)) * sum_{i=1}^{N_m}(pi * Xbari)
    double mean_of_vars{0.0}; // sum_{i=1}^{N_m}( (pi-1)*Vi )
    double sum_piX2{0.0}; // sum_{i=1}^{N_m}( pi*Xbari^2 )
    unsigned int sum_pi{0}; // sum(pi)

    // for each element in the cluster

    //  questa funzione non va bene se NON ho le covariate. il problema è l'ultimo termine, quello con beta e z.
    //  devo differenziare il calcolo di z in base a r=0 e r>0
    for(size_t ii = 0; ii <ind_i.size() ; ii++){

        const Individual& data_ji = data.at(ind_j[ii]).at(ind_i[ii]);

            //Rcpp::Rcout<<"Dentro compute_cluster_summaries"<<std::endl;
            //Rcpp::Rcout<<"data_ji.n_ji:"<<std::endl<<data_ji.n_ji<<std::endl;
            //Rcpp::Rcout<<"data_ji.Ybar_star_ji:"<<std::endl<<data_ji.Ybar_star_ji<<std::endl;
            //Rcpp::Rcout<<"data_ji.mean_ji:"<<std::endl<<data_ji.mean_ji<<std::endl;
        sum_pi += data_ji.n_ji;
        mean_of_vars += (double)(data_ji.n_ji - 1)*data_ji.var_ji;
        sum_piX2 += data_ji.n_ji * data_ji.mean_ji * data_ji.mean_ji;
        W += data_ji.n_ji * data_ji.mean_ji;

    }
    W = W/(double)sum_pi;

    //Rcpp::Rcout<<"W:"<<std::endl<<W<<std::endl;
    //Rcpp::Rcout<<"mean_of_vars:"<<std::endl<<mean_of_vars<<std::endl;
    //Rcpp::Rcout<<"sum_piX2:"<<std::endl<<sum_piX2<<std::endl;
    //Rcpp::Rcout<<"sum_pi:"<<std::endl<<sum_pi<<std::endl;
    return( std::tie(_W,_mean_of_vars,_sum_piX2,_sum_pi) );
}

*/