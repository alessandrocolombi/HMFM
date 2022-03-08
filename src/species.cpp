#ifndef __SPECIES_HPP__
#define __SPECIES_HPP__

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>
#include "include_headers.h"
#include "recurrent_traits.h"
#include "GSL_wrappers.h"
#include "utils.h"
#include <gsl/gsl_sf.h>
#include "ComponentPrior_factory.h"

//' Title Rcpp function
//'
//' @export
// [[Rcpp::export]]
int try_rcpp(int x){
	Rcpp::Rcout<<"Inside test function"<<std::endl;
	return x+10;
}

//' Raising Factorial
//'
//' \loadmathjax This function computes the rising factorial \mjseqn{(a)^{n}} using the gsl code for the Pochhammer symbol, i.e
//' \mjsdeqn{(a)^{n} = \frac{\Gamma(a+n)}{\Gamma(a)}}.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
//' @export
// [[Rcpp::export]]
double raising_factorial(const unsigned int& n, const double& a)
{
	return( gsl_sf_poch(a, (double)n ) );
}

//' log Raising Factorial
//'
//' \loadmathjax This function computes the logarithm of the rising factorial \mjseqn{(a)^{n}} using the gsl code for the log of Pochhammer symbol.
//' See \code{\link{raising_factorial}} and \code{\link{compute_Pochhammer}} for details. 
//' @export
// [[Rcpp::export]]
double log_raising_factorial(const unsigned int& n, const double& a) //secondo me troppo generale, può essere semplificata
{
	return( gsl_sf_lnpoch(a, (double)n ) );
}

//' Falling Factorial
//'
//' \loadmathjax This function computes the falling factorial \mjseqn{ a_{n} }. See \code{\link{raising_factorial}} for details.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
//' @export
// [[Rcpp::export]]
double my_falling_factorial(const unsigned int& n, const double& a)
{
	if(n%2 == 0) //n is even
		return( gsl_sf_poch(-a, (double)n ) );
	else //n is odd, change sign
		return( -1*gsl_sf_poch(-a, (double)n ) );
}

//' log Falling Factorial
//'
//' \loadmathjax This function computes the logarithm of the falling factorial \mjseqn{ a_{n} } using the gsl code for the log of Pochhammer symbol.
//' See \code{\link{my_falling_factorial}} and \code{\link{compute_Pochhammer}} for details. 
//' @export
// [[Rcpp::export]]
double my_log_falling_factorial(const unsigned int& n, const double& a)
{
	if(n%2 == 0) //n is even
		return( gsl_sf_lnpoch(-a, (double)n ) );
	else //n is odd, change sign
		return( -1*gsl_sf_lnpoch(-a, (double)n ) );
}

//' Pochhammer Symbol
//'
//' \loadmathjax This function computes the Pochhammer symbol, 
//' \mjsdeqn{(a)^{x} = \frac{\Gamma(a+x)}{\Gamma(a)}}.
//' Where \code{x} is a real number. When x is an integer, such a function coincides with the rising factorial defined in \code{\link{raising_factorial}}.
//' The raising (here denote with the upper apex) and the falling factorial (here denote with the lower apex) are related by the following relationship
//' \mjsdeqn{(a)_{n} = (-1)^{n}(a)^{n}}.
//' @export
// [[Rcpp::export]]
double compute_Pochhammer(const unsigned int& x, const double& a)
{
	return( gsl_sf_poch(a,x ) );
}

//' Pochhammer log Symbol
//'
//' \loadmathjax This function computes the Pochhammer symbol in log form. See \code{\link{compute_Pochhammer}} for details. 
//' @export
// [[Rcpp::export]]
double compute_log_Pochhammer(const unsigned int& x, const double& a)
{
	return( gsl_sf_lnpoch(a,x ) );
}

//' Build matrix of logC numbers
//'
//' This is the recursive function called by \code{\link{my_logC}} that builds the matrix containing all the log(|C(n,k)|) numbers.
//' It gets as input the element (n,k) to build, the scale s and location r (here defined as positive and non-negative numbers) and the
//' matrix res to be build. The matrix is constructed diagonal by diagonal, starting from the bottom.
//' Important remark, note that log(|C(n,0)|) = log(raising factorial (n,r)). 
//'
// [[Rcpp::export]]
void build_logC_matrix(const unsigned int& n, const unsigned int& k, const double& s, const double& r, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>& res){

	//Rcpp::Rcout<<"Costruisco n = "<<n<<", k = "<<k<<std::endl;
	if(k>n)
		throw std::runtime_error("Can not call build_logC_matrix if k > n");

	// Boundary conditions
	if(n == 0 && k == 0)
		res(n,k) = 0;
	else if(k == 0) 
		res(n,k) = std::log(raising_factorial(n,r));
	else if(n == k){
		build_logC_matrix(n-1, k-1, s, r, res);
		res(n,k) = std::log(s) + res(n-1,k-1);
	}
	else{
		double coef(s*k + r + n - 1);
		build_logC_matrix(n-1,k-1,s,r,res); // the matrix is always constructed diagonal by diagonal
		res(n,k) = std::log(coef) + res(n-1,k) + std::log( 1 + s/coef * std::exp( res(n-1,k-1) - res(n-1,k) ) );
	}
}

//' My logC
//'
//' This function is the recursive formula 2.69 in "Combinatorial methods in discrete distributions" book by Charalambides. 
//' It returns an (n+1 x n+1) matrix containing all the log(|C(nn,k)|) numbers, for nn = 0,...,n+1 and k = 0,...,nn. 
//' scale and location must be negative and non-positive, respectively.
//' As a consequence, it is memory expensive.
//' @export
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>	
my_logC(const unsigned int& n, const double& scale, const double& location){

	if(!( (scale<0) & (location<=0) ) )
		throw std::runtime_error("Error in my_logC. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative and location in non positive");
	
	const double& s = -scale; //s is strictly positive
	const double& r = -location; //r is non-negative

	MatCol res(MatCol::Constant(n+1,n+1,0.0)); //nxn matrix, n is in standard notation, i.e it starts from 1
	for(int k = n; k>=0; k--){
		//Rcpp::Rcout<<"Nel ciclo for, n = "<<n<<" , k = "<<k<<std::endl;
		build_logC_matrix(n,k,s,r,res);
		//Rcpp::Rcout<<"res finito k = "<<std::endl<<res<<std::endl;
	}
	return(res);
}


//' Compute log of absolute values of non Central C number 
//'
//' \loadmathjax This is the main function in the computation of C numbers. It uses the (2.69) formula in the "Combinatorial methods in discrete distributions" book.
//' It computes \mjseqn{log(|C(n,k; scale, location)|)} for each k=0,...,n. 
//' scale and location must be negative and non-positive, respectively.
//' It uses eigen objects, apparetly it is slower than using Rcpp vectors.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd my_logC2(const unsigned int& n, const double& scale, const double& location){

	if(!( (scale<0) & (location<=0) ) )
		throw std::runtime_error("Error in my_logC. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative and location in non positive");
	
	const double& s = -scale; //s is strictly positive
	const double& r = -location; //r is non-negative

	
	VecCol LogC_old(VecCol::Constant(n+1,0.0));

	if(n == 0)
		return LogC_old; // nothing to do in this case

	// Compute the first row
	LogC_old(0) = std::log(raising_factorial(1,r));
	LogC_old(1) = std::log(s);

	VecCol LogC_update(VecCol::Constant(n+1,0.0));
	double coef(0.0);
	for(std::size_t nn = 2; nn <= n; ++nn ){ //for each row
		LogC_update(0) = std::log(raising_factorial(nn,r));
		for(std::size_t k = 1; k < nn; ++k){ //for each column but the first and the last one
			coef = s*k + r + nn - 1;
			LogC_update(k) = std::log(coef) + LogC_old(k) + std::log( 1 + s/coef*std::exp( LogC_old(k-1) - LogC_old(k) ) );
		}
		LogC_update(nn) = nn*std::log(s); //update last element
		LogC_old.swap(LogC_update); //avoid copy, LogC_update has to be overwritten but in this way no copy is performed to update LogC_old.
	} 
	return (LogC_old);
}


//' My logC2 - central
//'
//' This function is the recursive formula 2.69 in "Combinatorial methods in discrete distributions" book by Charalambides for central numbers.
//' This is the specialization of \code{\link{my_logC2}} for central C numbers.
//' scale must be negative.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd my_logC2_central(const unsigned int& n, const double& scale){

	if(!( scale<0 ) )
		throw std::runtime_error("Error in my_logC2. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative.");
	
	const double& s = -scale; //s is strictly positive
	const double inf = std::numeric_limits<double>::infinity();

	VecCol LogC_old(VecCol::Constant(n+1,0.0));

	if(n == 0)
		return LogC_old; // nothing to do in this case

	// Compute the first row
	LogC_old(0) = -inf;
	LogC_old(1) = std::log(s);

	VecCol LogC_update(VecCol::Constant(n+1,0.0));
	LogC_update(0) = -inf;
	double coef(0.0);
	for(std::size_t nn = 2; nn <= n; ++nn ){ //for each row
		for(std::size_t k = 1; k < nn; ++k){ //for each column but the first and the last one
			coef = s*k + nn - 1;
			LogC_update(k) = std::log(coef) + LogC_old(k) + std::log( 1 + s/coef*std::exp( LogC_old(k-1) - LogC_old(k) ) );
		}
		LogC_update(nn) = nn*std::log(s); //update last element
		LogC_old.swap(LogC_update); //avoid copy, LogC_update has to be overwritten but in this way no copy is performed to update LogC_old.
	} 
	return (LogC_old);
}



//' compute_logC - Compute log of absolute values of non Central C number 
//'
//' \loadmathjax This is the main function in the computation of C numbers. It uses the (2.69) formula in the "Combinatorial methods in discrete distributions" book.
//' It computes \mjseqn{log(|C(n,k; scale, location)|)} for each k=0,...,n. 
//' This implementation uses Rcpp vectors.
//' @param scale must be strictly negative.
//' @param locatio must be non positive. Set to 0 for central C numbers.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector compute_logC(const unsigned int& n, const double& scale, const double& location){

	if(!( (scale<0) & (location<=0) ) )
		throw std::runtime_error("Error in my_logC. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative and location in non positive");
	
	const double& s = -scale; //s is strictly positive
	const double& r = -location; //r is non-negative

	
	Rcpp::NumericVector LogC_old(n+1, 0.0);

	if(n == 0)
		return LogC_old; // nothing to do in this case

	// Compute the first row
	LogC_old[0] = std::log(raising_factorial(1,r));
	LogC_old[1] = std::log(s);

	//Rcpp::NumericVector LogC_update(LogC_old);
	Rcpp::NumericVector LogC_update(n+1, 0.0);
	double coef(0.0);
	for(std::size_t nn = 2; nn <= n; ++nn ){ //for each row
		LogC_update[0] = std::log(raising_factorial(nn,r));
		for(std::size_t k = 1; k < nn; ++k){ //for each column but the first and the last one
			coef = s*k + r + nn - 1;
			LogC_update[k] = std::log(coef) + LogC_old[k] + std::log( 1 + s/coef*std::exp( LogC_old[k-1] - LogC_old[k] ) );
		}
		LogC_update[nn] = nn*std::log(s); //update last element
		
		//std::copy(LogC_update.begin(),LogC_update.end(),LogC_old.begin());
		std::swap(LogC_old, LogC_update); //avoid copy, LogC_update has to be overwritten but in this way no copy is performed to update LogC_old.
	} 
	return (LogC_old);
}


// Da Testare
double compute_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max = 100 ){
	
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Vprior, the length of n_i (group sizes) and gamma has to be equal");
	double res(0.0);
	for(std::size_t Mstar=0; Mstar <= M_max; ++Mstar){
		res += raising_factorial(k, (double)(Mstar+1) ) * qM.eval_prob(Mstar + k) * 
		       std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 1.0, std::multiplies<>(), 
		       					   [&Mstar, &k](const double& nj, const double& gamma_j){return 1/compute_Pochhammer(nj, gamma_j*(Mstar + k));} );		       
	}
	return res;
}


// Tutta da testare
double compute_log_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max = 100 ){
	
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Vprior, the length of n_i (group sizes) and gamma has to be equal");

	// Initialize vector of results
	std::vector<double> log_vect_res(M_max+1, -std::numeric_limits<double>::infinity() );
	// Initialize quantities to find the maximum
	unsigned int idx_max(0); 
	double val_max(log_vect_res[idx_max]);

	// Start the loop, let us compute all the elements
	for(std::size_t Mstar=0; Mstar <= M_max; ++Mstar){

		// Formula implementation
		log_vect_res[Mstar] = log_raising_factorial(k,(double)(Mstar+1) ) + 
							  qM.log_eval_prob(Mstar + k) - 
							  std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(), 
			       					   			  [&Mstar, &k](const double& nj, const double& gamma_j){return compute_log_Pochhammer(nj, gamma_j*(Mstar + k));} 
			       					   			);	
		// Check if it is the new maximum			       					   			   
        if(log_vect_res[Mstar]>val_max){ 
        	idx_max = Mstar;
        	val_max = log_vect_res[Mstar];
        } 
        	       					   			    
	}

	// Formula to compute the log of all the sums in a stable way
	return (val_max + 
			std::log(1 + 
				    std::accumulate(   log_vect_res.cbegin(), log_vect_res.cbegin()+idx_max, 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )  +
				    std::accumulate(   log_vect_res.cbegin()+idx_max+1, log_vect_res.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )  
		            ) 
		   ); 
}


//questa è sola per 2 gruppi
double compute_Kprior_unnormalized(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma){
	
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Kprior, the length of n_i (group sizes) and gamma has to be equal");
	if(n_i.size() != 2)
		throw std::runtime_error("Error in compute_Kprior, the length of n_i (group sizes) must be equal to 2");
	if(k == 0)
		return 0.0;

	double inf = std::numeric_limits<double>::infinity();

	std::vector<double> log_a(k+1,-inf);    // This vector contains all the quantities that depend only on r1 

	// Initialize quantities to find the maximum of log_a
	unsigned int idx_max1(0); 
	double val_max1(log_a[idx_max1]);
	
	// Compute all C numbers required
	Rcpp::NumericVector absC1 = compute_logC(n_i[0], -gamma[0], 0.0); //absC1[i] = |C(n1,i,-gamma1)| for i = 0,...,n1
	Rcpp::NumericVector absC2 = compute_logC(n_i[1], -gamma[1], 0.0); //absC2[i] = |C(n1,i,-gamma1)| for i = 0,...,n2

	// Start for loop
	for(std::size_t r1=0; r1 <= k; ++r1){

		// Compute a_r1 using its definition
		log_a[r1] = gsl_sf_lnchoose(k,r1) - my_log_falling_factorial(r1,k) + absC1[k-r1];

		// Prepare for computing the second term

		// Initialize vector of results
		std::vector<double> log_vect_res(k-r1+1, -inf );
		// Initialize quantities to find the maximum of log_vect_res
		unsigned int idx_max2(0); 
		double val_max2(log_vect_res[idx_max2]);
		
		// Inner loop on r2
		for(std::size_t r2=0; r2<= k-r1; ++r2){
			// Compute b_r2*c_r1r2
			log_vect_res[r2] = gsl_sf_lnchoose(k-r1,r2) - std::lgamma(k-r2+1) + absC2[k-r2];

			// Check if it is the new maximum of log_vect_res		       					   			   
        	if(log_vect_res[r2]>val_max2){ 
        		idx_max2 = r2;
        		val_max2 = log_vect_res[r2];
        	} 
		}

		// Update log_a:  log(a_i*alfa_i) = log(a_i) + log(alfa_i)
		log_a[r1] += val_max2 + 
					 std::log(1 + 
				              std::accumulate(   log_vect_res.cbegin(), log_vect_res.cbegin()+idx_max2, 0.0, [&val_max2](double& acc, const double& x){return acc + exp(x - val_max2 );}   )  +
				              std::accumulate(   log_vect_res.cbegin()+idx_max2+1, log_vect_res.cend(), 0.0, [&val_max2](double& acc, const double& x){return acc + exp(x - val_max2 );}   )  
		                     );
		// Check if it is the new maximum of log_a		       					   			   
       	if(log_a[r1]>val_max1){ 
       		idx_max1 = r1;
       		val_max1 = log_a[r1];
       	}			 
	}

	// Complete the sum over all elements in log_a
	return (val_max1 + 
			std::log(1 + 
				    std::accumulate(   log_a.cbegin(), log_a.cbegin()+idx_max1, 0.0, [&val_max1](double& acc, const double& x){return acc + exp(x - val_max1 );}   )  +
				    std::accumulate(   log_a.cbegin()+idx_max1+1, log_a.cend(), 0.0, [&val_max1](double& acc, const double& x){return acc + exp(x - val_max1 );}   )  
		            ) 
		   ); 
}




//' Test ComponentPrior
//' @export
// [[Rcpp::export]]
void Test_Prior(){
	Rcpp::String form = "Poisson";
	double lambda(5.0);
	ComponentPrior_Parameters param;
	param.Lambda = lambda;
	auto qM_ptr = Select_ComponentPrior(form, param);
	ComponentPrior& qM(*qM_ptr);
	std::string i_am = qM_ptr->showMe();
	std::string i_am2 = qM.showMe();
	Rcpp::Rcout<<"I am: "<<i_am<<std::endl;
	Rcpp::Rcout<<"I am - II - : "<<i_am2<<std::endl;
	form = "NegativeBinomial";
	param.p = 0.45;
	param.n_succ = 2;
	auto qM_ptr2 = Select_ComponentPrior(form, param);
	i_am = qM_ptr2->showMe();
	Rcpp::Rcout<<"I am: "<<i_am<<std::endl;
	// Devo testare se le densità sono giuste!
	Rcpp::Rcout<<"===================================="<<std::endl;
	Rcpp::Rcout<<"============== Densità ============="<<std::endl;
	Rcpp::Rcout<<"===================================="<<std::endl;
	Rcpp::Rcout<<"Poisson1(0)="<<qM_ptr->eval_prob(0)<<std::endl;
	Rcpp::Rcout<<"Poisson1(1)="<<qM_ptr->eval_prob(1)<<std::endl;
	Rcpp::Rcout<<"Poisson1(2)="<<qM_ptr->eval_prob(2)<<std::endl;
	Rcpp::Rcout<<"Poisson1(3)="<<qM_ptr->eval_prob(3)<<std::endl;
	Rcpp::Rcout<<"Poisson1(4)="<<qM_ptr->eval_prob(4)<<std::endl;
	Rcpp::Rcout<<"NegativeBinomial1(0)="<<qM_ptr2->eval_prob(0)<<std::endl;
	Rcpp::Rcout<<"NegativeBinomial1(1)="<<qM_ptr2->eval_prob(1)<<std::endl;
	Rcpp::Rcout<<"NegativeBinomial1(2)="<<qM_ptr2->eval_prob(2)<<std::endl;
	Rcpp::Rcout<<"NegativeBinomial(3)="<<qM_ptr2->eval_prob(3)<<std::endl;
	Rcpp::Rcout<<"NegativeBinomial(4)="<<qM_ptr2->eval_prob(4)<<std::endl;
	// Devo testare se le densità sono giuste!
	Rcpp::Rcout<<"===================================="<<std::endl;
	Rcpp::Rcout<<"============ Log Densità ==========="<<std::endl;
	Rcpp::Rcout<<"===================================="<<std::endl;
	Rcpp::Rcout<<"log_Poisson1(0)="<<qM_ptr->log_eval_prob(0)<<std::endl;
	Rcpp::Rcout<<"log_Poisson1(1)="<<qM_ptr->log_eval_prob(1)<<std::endl;
	Rcpp::Rcout<<"log_Poisson1(2)="<<qM_ptr->log_eval_prob(2)<<std::endl;
	Rcpp::Rcout<<"log_Poisson1(3)="<<qM_ptr->log_eval_prob(3)<<std::endl;
	Rcpp::Rcout<<"log_Poisson1(4)="<<qM_ptr->log_eval_prob(4)<<std::endl;
	Rcpp::Rcout<<"log_NegativeBinomial1(0)="<<qM_ptr2->log_eval_prob(0)<<std::endl;
	Rcpp::Rcout<<"log_NegativeBinomial1(1)="<<qM_ptr2->log_eval_prob(1)<<std::endl;
	Rcpp::Rcout<<"log_NegativeBinomial1(2)="<<qM_ptr2->log_eval_prob(2)<<std::endl;
	Rcpp::Rcout<<"log_NegativeBinomial1(3)="<<qM_ptr2->log_eval_prob(3)<<std::endl;
	Rcpp::Rcout<<"log_NegativeBinomial1(4)="<<qM_ptr2->log_eval_prob(4)<<std::endl;
	Rcpp::Rcout<<"===================================="<<std::endl;
	Rcpp::Rcout<<"================ Moda =============="<<std::endl;
	Rcpp::Rcout<<"===================================="<<std::endl;
	Rcpp::Rcout<<"La moda Poisson1 quando lambda = "<<lambda<<" è pari a "<<qM_ptr->get_mode()<<std::endl;
	Rcpp::Rcout<<"La moda NegativeBinomial1 quando p = "<<lambda<<" è pari a "<<param.p<<" e n è "<<param.n_succ<<" è pari a "<<qM_ptr2->get_mode()<<std::endl;
}


#endif
