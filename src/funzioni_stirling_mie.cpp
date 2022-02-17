#ifndef	__STIRLING__
#define __STIRLING__

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>
#include "include_headers.h"
#include "recurrent_traits.h"
#include "GSL_wrappers.h"

/*
	// [[Rcpp::depends(RcppArmadillo)]]
	#include <RcppArmadillo.h>
	#include <RcppArmadilloExtensions/sample.h>


	#include <limits>

	#include <exception>      	// to terminate a function by the command throw
								// throw std::runtime_error("message")
								//also includes std::exception, std::terminate
*/

//' Calcola Stirling
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calcola_stirling_old(double gamma, int n){


	Rcpp::NumericVector out(n); // The output vector

	double app; //just for  my convenience in programming
	double lgammak; // again just for my convenience
	Rcpp::NumericVector Lgammamjg_over_Lgammajg(n); // again just for my convenience

	for(int k=1;k<=n;k++){
		Lgammamjg_over_Lgammajg[k-1]=std::lgamma(n+k*gamma)-std::lgamma(k*gamma);
		lgammak=std::lgammaf(k+1);

		out[k-1]=0;
		for(int j=1;j<=k;j++){
			app = R::lchoose(k,j)+Lgammamjg_over_Lgammajg[j-1]-lgammak;
			out[k-1] += std::pow(-1,j-k)*std::exp(app);
		}
	}

	return(out);

}


//' Calcola Stirling Ricorsivo
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calcola_stirling_ricor_old(double gamma, unsigned int n){
	gamma=-gamma;
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

	return(row_jp1[Rcpp::Range(1,n)]);

}



//' Calcola Stirling Ricorsivo Abs
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calcola_stirling_ricor_abs_old(double gamma, unsigned int n){

	Rcpp::NumericVector row_j(n+1,0.0); // The output vector initialize all the element to zero
	row_j[0]=1; /// Row j=0

	Rcpp::NumericVector row_jp1(n+1,0.0); // The output vector initialize all the element to zero

	row_jp1[0]=0; /// Row j+1=1
	row_jp1[1]=gamma; /// Row j+1=1



	for(int j=1;j<n;j++){

		std::copy(row_jp1.begin(),row_jp1.end(),row_j.begin()); // row j+1 becomes row j!
		//Rcpp::Rcout<<"j="<<j<<" row_j="<<row_j<<"\n";

		//row_jp1[0]=std::exp(std::lgamma(j-1+1)-std::lgamma(j-1-j+1));
		for(int k=1;k<=(j+1);k++){
			row_jp1[k]=(gamma*k+j)*row_j[k]+gamma*row_j[k-1];
		}

	}

	return(row_jp1[Rcpp::Range(1,n)]);

}








//' Calcola Stirling Ricorsivo Log
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calcola_stirling_ricor_log_old(double gamma, unsigned int n){

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

	return(lrow_jp1[Rcpp::Range(1,n)]);

}







//' Calcola Fattoriale Generalizzato Ricorsivo Log
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calcola_generalized_factorial_ricor_log_old(double gamma, unsigned int n){

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
			lrow_jp1[k]= std::log(j-gamma*k)+lrow_j[k]+ std::log(1+gamma/(j-gamma*k)*std::exp(lrow_j[k-1]-lrow_j[k]));
		}

		lrow_jp1[j+1]= (j+1)*std::log(gamma);

	}

	return(lrow_jp1[Rcpp::Range(1,n)]);

}







/// The function here give as output a matrix
//' Calcola Fattoriale Generalizzato Ricorsivo Log Matrice
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calcola_generalized_factorial_ricor_log_matrice_old(double gamma, unsigned int n){

	double infinito = std::numeric_limits<double>::infinity();

	Rcpp:: NumericMatrix out(n+1,n+1); // The output matrix initialize all the element to -inf
	
	// Same as above, using STL fill
        std::fill(out.begin(), out.end(), -infinito);
	
	out(0,0)=0;  /// lRow j=0 SERVE?
	out(1,1)=std::log(gamma); /// lRow j+1=1


	for(int j=2;j<=n;j++){


		//Rcpp::Rcout<<"exp(-inf)="<< std::exp(lrow_j[0]-lrow_j[1])<<"\n";

		for(int k=1;k<j;k++){
			out(j,k) = std::log(j-1-gamma*k)+out(j-1,k)+ std::log(1+gamma/(j-1-gamma*k)*std::exp(out(j-1,k-1)-out(j-1,k)));

		}

		out(j,j)= (j)*std::log(gamma);

	}

	return(out(Rcpp::Range(1,n),Rcpp::Range(1,n)));

}





#endif
