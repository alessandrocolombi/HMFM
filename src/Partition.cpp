#include "Partition.h"

void Partition::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
  // From gs_data all needed variable are retrived
  unsigned int k = gs_data.K; // number of cluster
  unsigned int d = gs_data.d; // number of group
  unsigned int M = gs_data.M; // number of components
  const GDFMM_Traits::MatRow& S = gs_data.S; // Matrix of weights
  std::vector<unsigned int> n_j = gs_data.n_j;// number of observation per group
  const std::vector<double>& mu = gs_data.mu; // Vector of means
  const std::vector<double>& sigma = gs_data.sigma; // Vector of standard deviations

  std::set<unsigned int> s;// A set which will be useful for cluster
  GDFMM_Traits::VecRow probsvec(M);


  // Define data taken from gs_data
  const std::vector<std::vector<double>>& data = gs_data.data;
  // Initialization of probs_max
  double probs_max;
  //get sample index from GSL wrappers
  sample::sample_index sample_index;
  //clear C in order to reuse it
  C.clear();
  // Generate matrix of "weights" for each observation
  for(unsigned j=0; j<d; j++){
    std::vector<std::vector<double>> probs;
    std::vector<GDFMM_Traits::VecRow> probsmat;
    // matrix of probability -> we need one of them for every j
    for(unsigned i=0; i<n_j[j]; i++){
      std::vector<double> vec(M);
      for(unsigned m=0; m<M; m++){
        vec[m]=log(S(j,m)) + log_norm(data[j][i], mu[m], sigma[m]);
      }
      probs.push_back(vec);
      //probs is a matrix with a variable number of rows while fixed (=M) number of columns
      //get the maximum element of a row
      probs_max=*max_element(probs[i].begin(), probs[i].end());
      //probs è una matrice che ha numero di righe variabile ma sempre M colonne
      for(unsigned m=0; m<M; m++){
        probsvec(m)=exp(probs[i][m] - probs_max);
       // Rcpp::Rcout<<" p:"<<probsvec(m)<<" ";
      }
      probsmat.push_back(probsvec);
    }

    // Asegno tramite il sample su probs a ogni cluster un'etichetta
    //If M==1 populate C matrix with ones
    if (M == 1){
      std::vector<unsigned int> v(n_j[j], 1);
      C.push_back(v);
    }
    else{
      std::vector<unsigned int> dis;

      for (unsigned i=0; i<n_j[j]; i++) {
        dis.push_back(sample_index(gs_engine, probsmat[i]));
      }
      C.push_back(dis);
      /* C is a matrix with the same dimension as the matrix of data where at each position
      we have the label of each element*/
    }

  }

  //create vector of allocated components
  /*for(unsigned int j=0; j<d; j++){
   for(unsigned int i=0; i<n_j[j]; i++){
   std::cout<<C[j][i];
   }
   std::cout<<std::endl;
  }
   */

  // empty clust_out vector and set in order to reuse it
  clust_out.clear() ;
  s.clear();

  //Assign to each value of clust_out
  for(unsigned int j=0; j<d; j++){
    for(unsigned int i=0; i<n_j[j]; i++){
      s.insert(C[j][i]); //insert every label inside a set
      clust_out.assign(s.begin(),s.end()); //get the vector of the label sorted and newly labeled e.g (0-1-2-3)
    }
  }

  //Print allocated components
  Rcpp::Rcout << "Allocated components:";
  for (auto it = clust_out.begin(); it !=clust_out.end(); ++it)
    Rcpp::Rcout << ' ' << *it;
  Rcpp::Rcout << std::endl;


  k = clust_out.size(); //Set K=the size of clust out
  gs_data.K = k; // updating K in the struct gs_data
  gs_data.allocate_N(k); // initialize N according to new K
  gs_data.update_Ctilde(C, clust_out);
  for(unsigned m=0; m<gs_data.K; m++){
      Rcpp::Rcout<<gs_data.N_k[m]<< " ";
  }
}


double Partition::log_norm(double x, double u, double s) const {
  const double ONE_OVER_SQRT_2PI = 0.39894228040143267793994605993438;
  return log((ONE_OVER_SQRT_2PI/std::sqrt(s))*std::exp(-0.5*(x-u)*(x-u)/s));
}
