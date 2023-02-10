#include "FC_Mstar.h"

FC_Mstar::FC_Mstar( std::string na, unsigned int _proposalMstar, 
                    bool FixPart, bool _keepfixed) : FullConditional(na,_keepfixed), proposal(_proposalMstar), Partition_fixed(FixPart){

    support_proposal.reserve(2*proposal);
    for(int i = -(int)proposal; i <= (int)proposal; i++){
        if (i != 0) 
          support_proposal.push_back(i);
    }
};

void FC_Mstar::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
    //if(Partition_fixed){
         //Rcpp::Rcout << "Mstar is not updated because partition is fixed "<< std::endl;
    //}
    //else{
        // From gs_data all needed variable are retrived


        const unsigned int& k = gs_data.K;
        const unsigned int& d = gs_data.d;
        const double& lambda = gs_data.lambda;
        const std::vector<double>& gamma = gs_data.gamma; //gamma variables
        const double& log_sum = gs_data.log_sum;
        const GDFMM_Traits::MatUnsCol& N = gs_data.N; // dxK matrix; 
        
        if(d == 0){


            // Random sampler for the Poisson Distribution is created
            sample::rpoisson Poisson;

            // UPDATE ROUTINE
            double p0 = lambda/( lambda + k*exp(log_sum));
            // Select, via extraction from a uniform, which distribution sample from
            bool select_p0 = binary_decision(p0, gs_engine);
             //Rcpp::Rcout<<"elogsum:"<<exp(-log_sum)<<std::endl;
            // Rcpp::Rcout<<"p0:"<<p0<<std::endl;
            
            if(select_p0)
                gs_data.Mstar = Poisson(gs_engine, lambda*exp(-log_sum)) + 1;
            else
                gs_data.Mstar = Poisson(gs_engine, lambda*exp(-log_sum));

                    //Rcpp::Rcout<<" Sampled Mstar = "<<gs_data.Mstar<<std::endl;
            

            // Update M in the Gibbs Sampler
            gs_data.M = k + gs_data.Mstar;
        }
        else{
            /*
            sample::sample_index sample_index;
            sample::runif runif;

            GDFMM_Traits::VecRow log_probs(11);
            GDFMM_Traits::VecRow probs(11);

            for(size_t m = 0; m <= 10; m++){
                log_probs(m) = log_full_cond_Mstar(m, k, lambda, gamma, N);
            }
            double probs_max = log_probs.maxCoeff();

            // scale values of probs_vec
            for(size_t m = 0; m <= 10; m++){
                probs(m) = exp(log_probs(m) - probs_max);
             //Rcpp::Rcout<<" p:"<<probs_vec(m)<<" ";
                if(std::isnan(probs(m)))
                    throw std::runtime_error("Error in FC_Mstar.cpp, get a nan in log_probs ");
            }
             //Rcpp::Rcout<<std::endl;
            // Assign y_ji to a component sampling from a multinomial
            gs_data.Mstar = sample_index(gs_engine, probs);
            gs_data.M = k + gs_data.Mstar;
            */
            // ORIGNAL VERSION
            //get sample index from GSL wrappers
            sample::sample_index sample_index;
            sample::runif runif;

            // This function maps a positive integer number in a relative number (integer with sign)
            auto map_to_Z = [](const unsigned int& n){
                if(n % 2 == 0)
                    return ( (int)(n/2) );
                else
                    return ( -(int)(n+1)/2 );
            };
            // This function maps a relative number (integer with sign) in a positive integer number 
            auto map_to_N = [](const int& z){
                if( z < 0 )
                    return( (unsigned int)(2*std::abs(z) - 1) );
                else
                    return( (unsigned int)(2*std::abs(z)) );
            };
            int z = map_to_Z(gs_data.Mstar); // map the current value of Mstar in Z
            int Uproposal = support_proposal[ sample_index( gs_engine, 2*proposal ) ]; //draw uniform proposal from -proposal to proposal expect 0
            int zproposal = Uproposal; // initialize proposed value for z
            binary_decision(0.5, gs_engine) ? (zproposal += z) : (zproposal -= z) ; // complete proposal step
            // compute log probability of accepting the move
            //double log_alpha =  log_full_cond_Mstar( (unsigned int)map_to_N(zproposal), k, lambda, gamma, N  ) - 
                                //log_full_cond_Mstar( (unsigned int)map_to_N(z),         k, lambda, gamma, N  );
            //Rcpp::Rcout<<"Mstar = "<<gs_data.Mstar<<std::endl;
            //Rcpp::Rcout<<"Z = "<<z<<std::endl;
            //Rcpp::Rcout<<"Uproposal = "<<Uproposal<<std::endl;
            //Rcpp::Rcout<<"zproposal = "<<zproposal<<std::endl;
            //Rcpp::Rcout<<"Mstar_proposal = "<<map_to_N(zproposal)<<std::endl;
            double log_alpha = log_prob_MH( gs_data.Mstar, 
                                            (unsigned int)map_to_N(zproposal), 
                                            k, lambda, gamma, N );

            //Rcpp::Rcout<<"alpha = "<<std::min(1.0,std::exp(log_alpha))<<std::endl;
            // Acceptance-rejection move
            double u = runif(gs_engine);
            if( u  < std::min(1.0, std::exp(log_alpha))  )
                gs_data.Mstar = (unsigned int)map_to_N(zproposal);

            // Update M in the Gibbs Sampler
            gs_data.M = k + gs_data.Mstar;
            
        }    


        // Initialize S according to new M -> non va bene farlo qui
        //gs_data.allocate_S(gs_data.M);
        // Initialize tau according to new M
        //gs_data.allocate_tau(gs_data.M);

        // Quando fanno questa allocate, distruggono il contenuto di tau e S. 
        // Sembra vada bene perché poi aggiornano gamma, S, tau che sembrano non aver bisogno del contenuto di gs_data.S o gs_data.mu/sigma
    //}
}

// This function computes the full conditional distribution for Mstar when one does not condition on U_1,...,U_d
double FC_Mstar::log_full_cond_Mstar( const unsigned int& m, const unsigned int& K, const double& Lambda,
                                      const std::vector<double>& Gamma, const GDFMM_Traits::MatUnsCol& N)const
{
    double log_prod{0.0};
    unsigned int Mstar_plus_K{m+K}; // compute Mstar+K
    unsigned int d{Gamma.size()};   // get number of levels

    // This for computes: log( prod_{j=1}^{d}( 1/(gamma_j*(Mstar+K))_{n_j} ) ) = - sum_{j=1}^{d}( log( gamma_j*(Mstar+K))_{n_j} ) )
    for(size_t j = 0; j < d; j++ ){
        log_prod -=  log_raising_factorial(N.row(j).sum(), Gamma[j]*(double)Mstar_plus_K);
    }
    return( std::log(Mstar_plus_K) + (double)m*std::log(Lambda) - gsl_sf_lnfact(m) + log_prod );
}

// This function computes the log acceptance probability for the MH move
double FC_Mstar::log_prob_MH( const unsigned int& m, const unsigned int& m_new, 
                              const unsigned int& K, const double& Lambda,
                              const std::vector<double>& Gamma, const GDFMM_Traits::MatUnsCol& N)const
{

    //Rcpp::Rcout<<std::endl;
    //Rcpp::Rcout<<"------------------------"<<std::endl;
    //Rcpp::Rcout<<std::endl;
    //Rcpp::Rcout<<"Mstar = "<<m<<std::endl;
    //Rcpp::Rcout<<"Mstar_proposed = "<<m_new<<std::endl;
    //Rcpp::Rcout<<"K = "<<K<<std::endl;
    //Rcpp::Rcout<<"Lambda = "<<Lambda<<std::endl;

    unsigned int d{Gamma.size()};   // get number of levels
    int q = (int)m_new - (int)m;    // get q, the difference with sign between the currest value of Mstar and the proposed value
    
    //double log_prod{0.0};
    //double quarto_termine{0.0};
    // first term: log( (m_new+k)/(m+k) )
    double res{ std::log(  (double)(m_new + K)  ) - std::log(  (double)(m+K) ) };
    //Rcpp::Rcout<<"first term = "<<std::log(  (double)(m_new + K)  ) - std::log(  (double)(m+K) )<<std::endl;
    // second term
    res += (double)q*std::log(Lambda);
    //Rcpp::Rcout<<"second term = "<<(double)q*std::log(Lambda)<<std::endl;
    // This for computes: log( prod_{j=1}^{d}( (gamma_j*(m+K))_{n_j}/(gamma_j*(m_new+K))_{n_j} ) ) 
    for(size_t j = 0; j < d; j++ ){
        res +=  log_raising_factorial(N.row(j).sum(), Gamma[j]*(double)(m+K)) - 
                log_raising_factorial(N.row(j).sum(), Gamma[j]*(double)(m_new+K)) ;
        //log_prod +=  log_raising_factorial(N.row(j).sum(), Gamma[j]*(double)(m+K)) - 
                     //log_raising_factorial(N.row(j).sum(), Gamma[j]*(double)(m_new+K)) ;
    }
    //Rcpp::Rcout<<"third term = "<<log_prod <<std::endl;
    // last term: ratio of log factorials
    (q > 0) ? (  res -= log_raising_factorial(std::abs(q), (double)(m+1))  ) : (  res += log_raising_factorial(std::abs(q), (double)( m+1-std::abs(q)) )  ) ; 
    //(q > 0) ? (  quarto_termine = -log_raising_factorial(std::abs(q), (double)(m+1))  ) : (  quarto_termine = log_raising_factorial(std::abs(q), (double)( m+1-std::abs(q)) )  ) ; 
    //Rcpp::Rcout<<"Fourth term = "<<quarto_termine<<std::endl;
    return res;
}