#ifndef __RAF_STIRLING__
#define __RAF_STIRLING__

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>
#include "include_headers.h"
#include "recurrent_traits.h"
#include "GSL_wrappers.h"

//' Falling Factorial - Raf
//'
//' @export
// [[Rcpp::export]]
double falling_factorial(double x, int n)
{
  // if n=0 x_0=1 quindi log(x_0)=0
  if(n==0){return(0);
    }
  double out=0;
  int sign=+1; //sign of 
  double factor;
  for(int i=0;i<n;i++){
     factor = x-i;
      sign = sign*pow(-1,1*std:: signbit(factor)); 
      //Rcpp::Rcout<<"sign="<<sign<<"\n";
       // signbit returns 1 if the argument is negative 0 otherwise
      out += std::log(std::abs(factor)); 
       }
  
  return(sign*exp(out));
}


//' Compute Non Central C number - Direct formula Raf
//'
//' Compute the Charalambides non central C numbers using the direct formula, i.e the (2.60) one in the "Combinatorial methods in discrete distributions" book.
//' This version works for r and gamma real numbers, not only positive. However, it is unstable.
//' It computes the whole sequence, C(n,k,gamma,r) for k=0,...,n
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calcola_stirling(int n, double gamma, double r)
{
  
  
  Rcpp::NumericVector out(n+1); // The output vector
  
  double app; //just for  my convenience in programming
  double lgammak; // again just for my convenience
  Rcpp::NumericVector falling(n+1);   // again just for my convenience
  
  for(int k=0;k<=n;k++){
    falling[k]=falling_factorial(gamma*k+r,n);
 // log_ratio_gamma[k]=std::lgamma(k*gamma+r+1)-std::lgamma(k*gamma+r+1-n);
    //Sopra ho usato il fatto che (x)_n=n!(x\choose n) (guarda wiky) 
   // Rcpp::Rcout<<"falling="<<falling[k]<<"\n";
    //Rcpp::Rcout<<"k="<<k<<" n="<<n<<" gamma="<<gamma<<" r="<<r<<" \n";
    lgammak=std::lgamma(k+1); //log di k!
    
    out[k]=0;
    for(int j=0;j<=k;j++){
      app = R::lchoose(k,j)-lgammak;
      //Rcpp::Rcout<<"k="<<k<<" j="<<j<<"\n";
      out[k] += std::pow(-1,k-j)*std::exp(app)*falling[j];
    }
    out[k]= out[k];
  }
  
  return(out);
}



//' Compute Non Central C number - Recursive formula Raf
//'
//' Compute the Charalambides non central C numbers using the recursive formula, i.e the (2.67) one in the "Combinatorial methods in discrete distributions" book.
//' This version works for r and gamma real numbers, not only positive. However, it is not in log scale, therefore its usage is limited to small n cases.
//' It computes the whole sequence, C(n,k,gamma,r) for k=0,...,n.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calcola_stirling_ricor(unsigned int n, double gamma, double r)
{
 //  gamma=-gamma;
  Rcpp::NumericVector row_j(n+1,0.0); // The output vector initialize all the element to zero
  row_j[0]=1; /// 
  
  
  Rcpp::NumericVector row_jp1(n+1,0.0); // The output vector initialize all the element to zero
  row_jp1[0]=falling_factorial(r,1);
  row_jp1[1]=gamma; /// Row j+1=1
  
  
  
  for(int j=1;j<n;j++){
    
    
    std::copy(row_jp1.begin(),row_jp1.end(),row_j.begin()); // row j+1 becomes row j!
    //Rcpp::Rcout<<"j="<<j<<" row_j="<<row_j<<"\n";
    row_jp1[0]=falling_factorial(r,j+1);
    
    for(int k=1;k<=(j+1);k++){
      row_jp1[k]=(gamma*k+r-j)*row_j[k]+gamma*row_j[k-1];
    }
    
  }
  
  return(row_jp1[Rcpp::Range(0,n)]);  
}



//' Compute Central C number - Recursive formula Raf
//'
//' Same as \code{\link{calcola_stirling_ricor}} but with r=0, i.e for central numbers only.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calcola_stirling_ricor_centr(double gamma, unsigned int n)
{
  Rcpp::NumericVector row_j(n+1,0.0); // The output vector initialize all the element to zero
  row_j[0]=1; /// Row j=0
  
  Rcpp::NumericVector row_jp1(n+1,0.0); // The output vector initialize all the element to zero
  row_jp1[1]=gamma; /// Row j+1=1
  
  
  
  for(int j=1;j<n;j++){
    
    std::copy(row_jp1.begin(),row_jp1.end(),row_j.begin()); // row j+1 becomes row j!
    //Rcpp::Rcout<<"j="<<j<<" row_j="<<row_j<<"\n";
    
    
    for(int k=1;k<=(j+1);k++){
      row_jp1[k]=(gamma*k-j)*row_j[k]+gamma*row_j[k-1];
    }
    
  }
  
  return(row_jp1[Rcpp::Range(0,n)]);
}


//' Compute log of absolute values of non Central C number - Recursive formula for Raf
//'
//' \loadmathjax This is the main function in the computation of C numbers. It uses the (2.69) formula in the "Combinatorial methods in discrete distributions" book.
//' It computes \mjseqn{|C(n,k;-\gamma,-r)|} for each k=0,...,n. \mjseqn{\gamma} and r have to be positive. 
//' Attenzione nota importante. Questa in posizione (n,0) ha sempre lo stesso valore ed è pari al fattoriale decrescente (r)1. Ci sono due cose che non tornano,
//' prima di tutto il fatto che sia sempre uguale e poi che sia il fattoriale decrescente. In base alla relazione 
//' \mjsdeqn{(a + n -1)_{n} = (a)^{n}}, quello dovrebbe essere un raising factorial.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calcola_stirling_ricor_log( unsigned int n, double gamma, double r)
{
  
  double infinito = std::numeric_limits<double>::infinity();
  
  Rcpp::NumericVector lrow_j(n+1,-infinito); // The output vector initialize all the element to zero
  lrow_j[0]=0; /// lRow j=0
  
  Rcpp::NumericVector lrow_jp1(n+1,-infinito); // The output vector initialize all the element to zero
  lrow_jp1[0]=std::log(falling_factorial(r,1));
  lrow_jp1[1]=std::log(gamma); /// lRow j+1=1
  
  
  for(int j=1;j<n;j++){
    
    std::copy(lrow_jp1.begin(),lrow_jp1.end(),lrow_j.begin()); // lrow j+1 becomes lrow j!
    //Rcpp::Rcout<<"j="<<j<<" lrow_j="<<lrow_j<<"\n";
    
    //Rcpp::Rcout<<"exp(-inf)="<< std::exp(lrow_j[0]-lrow_j[1])<<"\n";

    // NB così è ancora sbagliata. qua non ci va il falling ma il rising factorial!
    lrow_jp1[0]=std::log(falling_factorial(r,j)); //il primo elemento non è mai aggiornato. aggiungo riga per confronto tempi. 
    for(int k=1;k<=(j);k++){
      lrow_jp1[k]= std::log(gamma*k+r+j)+lrow_j[k]+ std::log(1+gamma/(gamma*k+r+j)*std::exp(lrow_j[k-1]-lrow_j[k]));
    }
    
    lrow_jp1[j+1]= (j+1)*std::log(gamma);
    
  }
  
  return(lrow_jp1[Rcpp::Range(0,n)]);
}

//' Compute log of absolute values of Central C number - Recursive formula for Raf
//'
//' Same as \code{\link{calcola_stirling_ricor_log}} but with r=0, i.e for central numbers only.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calcola_stirling_ricor_log_centrali(unsigned int n, double gamma)
{
  
  double infinito = std::numeric_limits<double>::infinity();
  
  Rcpp::NumericVector lrow_j(n+1,-infinito); // The output vector initialize all the element to zero
  lrow_j[0]=0; /// lRow j=0
  
  Rcpp::NumericVector lrow_jp1(n+1,-infinito); // The output vector initialize all the element to zero
  lrow_jp1[1]=std::log(gamma); /// lRow j+1=1
  
  
  
  for(int j=1;j<n;j++){
    
    std::copy(lrow_jp1.begin(),lrow_jp1.end(),lrow_j.begin()); // lrow j+1 becomes lrow j!
    //Rcpp::Rcout<<"j="<<j<<" lrow_j="<<lrow_j<<"\n";
    
    //Rcpp::Rcout<<"exp(-inf)="<< std::exp(lrow_j[0]-lrow_j[1])<<"\n";
    
    for(int k=1;k<=(j);k++){
      lrow_jp1[k]= std::log(gamma*k+j)+lrow_j[k]+ std::log(1+gamma/(gamma*k+j)*std::exp(lrow_j[k-1]-lrow_j[k]));
    }
    
    lrow_jp1[j+1]= (j+1)*std::log(gamma);
    
  }
  
  return(lrow_jp1[Rcpp::Range(0,n)]);  
}


#endif