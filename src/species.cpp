#include "species.h"

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	log_stable_sum
//------------------------------------------------------------------------------------------------------------------------------------------------------

// These functions computes log(sum_i(a_i)) using a stable formula for log values. Let us denote a* to the the maximum value of vector a which is attained when i = i*.
// ---> log(sum_i(a_i)) = log(a*) + log[ 1 + sum_{i not i*}(exp{log(a_i) - log(a*)}) ]
// See that only the logarithm of the elements of a are needed. Hence, it is likely that one has already computed them in log scale. If so, set is_log = T
// In this versione of the function, the maximum and its position are passed as parameters. No check is performed to be sure that such information were true or not.
// It is assumed that the value in val_max is on the same scale of the values of a, i.e it is in log scale if is_log is set to true.
double log_stable_sum(const std::vector<double>& a, const bool is_log, const double& val_max, const unsigned int& idx_max){
	
	double inf = std::numeric_limits<double>::infinity();
	
	if(a.size() == 0)
		return 0.0;

	if(is_log==TRUE){ // a contains values in log scale
		
		if(val_max == -inf) // need to handle the case where all values are -inf
			return -inf;

		return (val_max +
				std::log(1 +
					    std::accumulate(   a.cbegin(), a.cbegin()+idx_max, 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )  +
					    std::accumulate(   a.cbegin()+idx_max+1, a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )
				        )
			   );
	}
	else{

		if(val_max < 0)
			throw std::runtime_error("log_stable_sum, if is_log is FALSE, the maximum value can not be negative ");
		if(val_max == 0)
			return 0.0;
		// Do not checks if values of a are strictly positive
		return ( std::log(val_max) +
				 std::log(1 +
					      std::accumulate(   a.cbegin(), a.cbegin()+idx_max, 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   )  +
					      std::accumulate(   a.cbegin()+idx_max+1, a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   )
			             )
			   );	
	}
}

// As before but gets an iterator poining to the maximum value
double log_stable_sum(const std::vector<double>& a, const bool is_log, const GDFMM_Traits::vector_d_citerator& it_max){

	double inf = std::numeric_limits<double>::infinity();

	if(a.size() == 0)
		return 0.0;
	double val_max{*it_max};
	if(is_log){ // a contains values in log scale
		
		if(val_max == -inf) // need to handle the case where all values are -inf
			return -inf;

		return (    val_max +
					std::log(1 +
						    std::accumulate(   a.cbegin(), it_max, 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )  +
						    std::accumulate(   it_max+1, a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )
				            )
			   );
	}
	else{

		if(val_max < 0)
			throw std::runtime_error("log_stable_sum, if is_log is FALSE, the maximum value can not be negative ");
		if(val_max == 0)
			return 0.0;

		// Do not checks if values of a are strictly positive
		return ( std::log(val_max) +
				 std::log(1 +
					      std::accumulate(   a.cbegin(), it_max, 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   )  +
					      std::accumulate(   it_max+1, a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   )
			             )
			   );	
	}
}

// In this specialized version of the function, the position of the max value is not provided. Hence, one additional operation is done. 
// Since exp( log(a*) - log(a*)) = 1, the previous formula becomes 
// ---> log(sum_i(a_i)) = log(a*) + log[ sum_{i}(exp{log(a_i) - log(a*)}) ]
double log_stable_sum(const std::vector<double>& a, const bool is_log, const double& val_max){

	double inf = std::numeric_limits<double>::infinity();

	if(a.size() == 0)
		return 0.0;

	// Do not checks if it is really the max value
	if(is_log){ // a contains values in log scale

		if(val_max == -inf) // need to handle the case where all values are -inf
			return -inf;

		return (val_max +
					std::log( std::accumulate(   a.cbegin(), a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )   )
			   );
	}
	else{

		if(val_max < 0)
			throw std::runtime_error("log_stable_sum, if is_log is FALSE, the maximum value can not be negative ");
		if(val_max == 0)
			return 0.0;

		return ( std::log(val_max) +
				 std::log( std::accumulate(   a.cbegin(), a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   ) )
			   );	
	}

}

// In this version of the formula, the maximum value is computed
double log_stable_sum(const std::vector<double>& a, const bool is_log){
	if(a.size() == 0)
		return 0.0;

	// Computes maximum value
	GDFMM_Traits::vector_d_citerator it_max{std::max_element(a.cbegin(), a.cend())};
	//double val_max{*it_max};
	// Calls the specialized version
	return log_stable_sum(a,is_log,it_max);
}


int try_rcpp(int x){
	Rcpp::Rcout<<"Inside test function"<<std::endl;
	return x+10;
}

double raising_factorial_poch(const unsigned int& n, const double& a)
{
	return( gsl_sf_poch(a, (double)n ) );
}

double log_raising_factorial_poch(const unsigned int& n, const double& a) //secondo me troppo generale, può essere semplificata
{	
	if(a<=0)
		throw std::runtime_error("Error in log_raising_factorial, can not compute the raising factorial of a negative number in log scale");

	return( gsl_sf_lnpoch(a, (double)n ) );
}


double log_raising_factorial(const unsigned int& n, const double& a){

	if(n==0)
		return 0.0;
	if(a<0)
		throw std::runtime_error("Error in my_log_raising_factorial, can not compute the raising factorial of a negative number in log scale"); // mi si blocca se a è 0 ma questo caso è da gestire nel caso del cacolo dei numeri centrali!
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
			if(i>100)
				return -1.0;
		}
		return val_max*res;

	}
}

double raising_factorial(const unsigned int& n, const double& a){
	if(n==0)
		return 1.0;
	if(n==1)
		return a;
	if(a<=0){
		double res{1.0};
		for(std::size_t i = 0; i <= n-1; ++i){
			res *= ( a + (double)i ) ;
		}
		return res;
	}
	else{
		return std::exp(log_raising_factorial(n,a));
	}

}


double my_falling_factorial_old(const unsigned int& n, const double& a)
{
	if(n%2 == 0) //n is even
		return( gsl_sf_poch(-a, (double)n ) );
	else //n is odd, change sign
		return( -1*gsl_sf_poch(-a, (double)n ) );
}

double my_log_falling_factorial_old(const unsigned int& n, const double& a) //questo non va per a negativi!
{
	if(n%2 == 0) //n is even
		return( gsl_sf_lnpoch(-a, (double)n ) );
	else //n is odd, change sign
		return( -1*gsl_sf_lnpoch(-a, (double)n ) );
}

double my_log_falling_factorial(const unsigned int& n, const double& a) 
{
	if(n==0)
		return 0.0;
	if(a<0)
		throw std::runtime_error("Error in my_log_falling_factorial, can not compute the falling factorial of a negative number in log scale");
	else if(a==0.0){
		return -std::numeric_limits<double>::infinity();
	}
	if(a-n+1<=0)
		throw std::runtime_error("Error in my_log_falling_factorial, can not compute the falling factorial (a)_n in log scale if a <= n-1");
	else{
		double val_max{std::log(a)};
		double res{1.0};
		if (n==1)
			return val_max;
		for(std::size_t i = 1; i <= n-1; ++i){
			res += std::log(a - (double)i) / val_max;
		}
		return val_max*res;
	}
}

double my_falling_factorial(const unsigned int& n, const double& a)
{
	if(n==0)
		return 1.0;
	if(n==1)
		return a;
	if(a<=0){
		double res{1.0};
		for(std::size_t i = 0; i <= n-1; ++i){
			res *= ( a + (double)i ) ;
		}
		return res;
	}
	else{
		return std::exp(my_log_falling_factorial(n,a));
	}
}



double compute_Pochhammer(const unsigned int& x, const double& a)
{
	return( gsl_sf_poch(a,x ) );
}

double compute_log_Pochhammer(const unsigned int& x, const double& a)
{
	return( gsl_sf_lnpoch(a,x ) );
}

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


Rcpp::NumericVector compute_logC(const unsigned int& n, const double& scale, const double& location){

	if(!( (scale<0) & (location<=0) ) )
		throw std::runtime_error("Error in my_logC. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative and location in non positive");

	const double& s = -scale; //s is strictly positive
	const double& r = -location; //r is non-negative


	Rcpp::NumericVector LogC_old(n+1, 0.0);

	if(n == 0)
		return LogC_old; // nothing to do in this case

	// Compute the first row
	LogC_old[0] = log_raising_factorial(1,r);
	LogC_old[1] = std::log(s);

	//Rcpp::NumericVector LogC_update(LogC_old);
	Rcpp::NumericVector LogC_update(n+1, 0.0);
	double coef(0.0);
	for(std::size_t nn = 2; nn <= n; ++nn ){ //for each row
		LogC_update[0] = log_raising_factorial(nn,r);
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
double compute_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max ){

	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Vprior, the length of n_i (group sizes) and gamma has to be equal");
	double res(0.0);
	for(std::size_t Mstar=0; Mstar <= M_max; ++Mstar){
		res += raising_factorial(k, (double)(Mstar+1) ) * qM.eval_prob(Mstar + k) *
		       std::inner_product( n_i.cbegin(), n_i.cend(),gamma.cbegin(), 1.0, std::multiplies<>(),
		       					   [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return 1/raising_factorial( nj, gamma_j*((double)Mstar + (double)k));} );
		        // nj is an integer, this is just a raising factorial, not a pochammer
	}
	return res;
}


// Sembra funzionare bene anche nel caso d=2
double compute_log_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max ){

	//Rcpp::Rcout<<"Dentro compute_log_Vprior"<<std::endl;
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) must be positive");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Vprior, the length of n_i (group sizes) and gamma has to be equal");

	// Initialize vector of results
	std::vector<double> log_vect_res(M_max+1, -std::numeric_limits<double>::infinity() );
	// Initialize quantities to find the maximum
	unsigned int idx_max{0};
	double val_max(log_vect_res[idx_max]);

	// Start the loop, let us compute all the elements
	for(std::size_t Mstar=0; Mstar <= M_max; ++Mstar){
		//Rcpp::Rcout<<"Mstar = "<<Mstar<<std::endl;
		// Formula implementation
		
		log_vect_res[Mstar] = log_raising_factorial(k,Mstar+1 ) +
							  qM.log_eval_prob(Mstar + k) -
							  std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
			       					   			  [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return log_raising_factorial( nj, gamma_j*(Mstar + k) );}
			       					   			);
		/*
		double _a = log_raising_factorial(k, Mstar+1 );	
		Rcpp::Rcout<<"log_raising_factorial = "<<_a<<std::endl;		       					   			
		double _b = qM.log_eval_prob(Mstar + k);			
		Rcpp::Rcout<<"qM.log_eval_prob = "<<_b<<std::endl;		       					   			
		Rcpp::Rcout<<" ORA PARTE IL SECONDO log_raising_factorial !"<<std::endl;		       					   			
		double _c = std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
			       					   			  [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return log_raising_factorial( nj, gamma_j*(Mstar + k) );}
			       					   );
		Rcpp::Rcout<<"inner_product = "<<_c<<std::endl;	
		log_vect_res[Mstar] = _a +_b -_c;		
		*/
		//Rcpp::Rcout<<" ---> calcolato log_vect_res[Mstar] = "<<log_vect_res[Mstar]<<std::endl;	 
		      					   			       					   		
		// Check if it is the new maximum
        if(log_vect_res[Mstar]>val_max){
        	idx_max = Mstar;
        	val_max = log_vect_res[Mstar];
        }

	}

	// Formula to compute the log of all the sums in a stable way
	return log_stable_sum(log_vect_res, TRUE, val_max, idx_max);
	/*
	return (val_max +
			std::log(1 +
				    std::accumulate(   log_vect_res.cbegin(), log_vect_res.cbegin()+idx_max, 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )  +
				    std::accumulate(   log_vect_res.cbegin()+idx_max+1, log_vect_res.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )
		            )
		   );
	*/	   
}


//questa è sola per 1 o 2 gruppi
double compute_Kprior_unnormalized(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma){

	//Rcpp::Rcout<<"Dentro a compute_Kprior_unnormalized con n_i: ";
	//for(auto __v : n_i)
		//Rcpp::Rcout<<__v<<", ";
	//Rcpp::Rcout<<std::endl;


	double inf = std::numeric_limits<double>::infinity();

	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Kprior, the length of n_i (group sizes) and gamma has to be equal");
	if(n_i.size() > 2 || n_i.size() == 0)
		throw std::runtime_error("Error in compute_Kprior, the length of n_i (group sizes) must be equal to 1 or 2");
	if(k == 0)
		return -inf;
	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  )
		return -inf;

	if(n_i.size()==1){ // one group only 
		Rcpp::NumericVector absC = compute_logC(n_i[0], -gamma[0], 0.0); //absC[i] = |C(n1,i,-gamma1)| for i = 0,...,n1
		return absC[k];
	}
	else{ // two groups case

		// Compute all C numbers required
		Rcpp::NumericVector absC1 = compute_logC(n_i[0], -gamma[0], 0.0); //absC1[i] = |C(n1,i,-gamma1)| for i = 0,...,n1
		Rcpp::NumericVector absC2 = compute_logC(n_i[1], -gamma[1], 0.0); //absC2[i] = |C(n2,i,-gamma2)| for i = 0,...,n2

		const int start1 = std::max( {0, (int)k - (int)n_i[0]} );
		const int start2 = std::max( {0, (int)k - (int)n_i[1]} );
		const int end1   = std::min( {(int)k, (int)n_i[1]} );
					//Rcpp::Rcout<<"start1 = "<<start1<<std::endl; 
					//Rcpp::Rcout<<"start2 = "<<start2<<std::endl; 
					//Rcpp::Rcout<<"end1   = "<<end1<<std::endl; 
  
					//std::vector<double> log_a(k+1, -inf);    // old version
		std::vector<double> log_a(end1-start1+1, -inf);    // This vector contains all the quantities that depend only on r1
		// Initialize quantities to find the maximum of log_a
		unsigned int idx_max1(0);
		double val_max1(log_a[idx_max1]);

		// Start for loop
		unsigned int outer_indx{0};
		for(std::size_t r1=start1; r1 <= end1; ++r1){

			// Compute a_r1 using its definition
			log_a[outer_indx] = gsl_sf_lnchoose(k,r1) - my_log_falling_factorial(r1,(double)k) + absC1[k-r1];
			
			// Prepare for computing the second term
		
			// Initialize vector of results
					//std::vector<double> log_vect_res(k-r1+1, -inf ); // old version
			std::vector<double> log_vect_res(k-r1-start2+1, -inf );

			// Initialize quantities to find the maximum of log_vect_res
			unsigned int idx_max2(0);
			double val_max2(log_vect_res[idx_max2]);
		
			// Inner loop on r2
			unsigned int inner_indx{0};
			for(std::size_t r2=start2; r2<= k-r1; ++r2){

				// Compute b_r2*c_r1r2
				log_vect_res[inner_indx] = gsl_sf_lnchoose(k-r1,r2) + std::lgamma(k-r2+1) +  absC2[k-r2];
				
				// Check if it is the new maximum of log_vect_res
	        	if(log_vect_res[inner_indx]>val_max2){
	        		idx_max2 = inner_indx;
	        		val_max2 = log_vect_res[inner_indx];
	        	}
	        			//Rcpp::Rcout<<"Computing for r1 = "<<r1<<" and r2 = "<<r2<<std::endl;
				inner_indx++;
			}
		
		 
			// Update log_a:  log(a_i*alfa_i) = log(a_i) + log(alfa_i)
			log_a[outer_indx] += log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2);
			      
			// Check if it is the new maximum of log_a
	       	if(log_a[outer_indx]>val_max1){
	       		idx_max1 = outer_indx;
	       		val_max1 = log_a[outer_indx];
	       	}
	       	outer_indx++;
		}

		// Complete the sum over all elements in log_a
		return log_stable_sum(log_a, TRUE, val_max1, idx_max1);

	}
}

double compute_Kprior_unnormalized_recursive(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma){

	double inf = std::numeric_limits<double>::infinity();

	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Kprior_unnormalized_recursive, the length of n_i (group sizes) and gamma has to be equal");
	if( n_i.size() == 0)
		throw std::runtime_error("Error in compute_Kprior_unnormalized_recursive, the length of n_i (group sizes) must be greater than 0");
	if( n_i.size() <= 2)
		return compute_Kprior_unnormalized(k,n_i,gamma);
	if(k == 0)
		return -inf;
	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  )
		return -inf;

	
	// if here, n_i.size()>2 
	std::vector<double> log_a(k+1, -inf);    
	// Initialize quantities to find the maximum of log_a
	unsigned int idx_max1(0);
	double val_max1(log_a[idx_max1]);

	// Start for loop
	for(std::size_t k1=0; k1 <= k; ++k1){

		Rcpp::Rcout<<"----> Dentro a k1 = "<<k1<<std::endl;
		// Compute recursive formula
		log_a[k1] = compute_Kprior_unnormalized_recursive(k1, {n_i.cbegin(), n_i.cend()-1} , {gamma.cbegin(), gamma.cend()-1} ); 
		Rcpp::Rcout<<"log_a["<<k1<<"] = "<<"NNP(K^(2) = "<<k1<<") = "<<log_a[k1]<<std::endl;
		// Prepare for computing the second term
		
		// Initialize vector of results
		std::vector<double> log_vect_res(k1+1, -inf ); 

		// Initialize quantities to find the maximum of log_vect_res
		unsigned int idx_max2(0);
		double val_max2(log_vect_res[idx_max2]);
		
		// Inner loop on r2
		unsigned int inner_indx{0};
		for(std::size_t k2=k-k1; k2<= k; ++k2){
			Rcpp::Rcout<<"***** Dentro a k2 = "<<k2<<std::endl;
			
			// Compute last element
			Rcpp::Rcout<<"n_i[n_i.size()-1] = "<<n_i[n_i.size()-1]<<std::endl;
			Rcpp::Rcout<<"gamma[gamma.size()-1] = "<<gamma[gamma.size()-1]<<std::endl;
			Rcpp::Rcout<<"NNP(K^(1) = "<<k2<<") = "<< compute_Kprior_unnormalized(  k2, {n_i[n_i.size()-1]}, {gamma[gamma.size()-1]}  ) <<std::endl;
			Rcpp::Rcout<<"coef = "<<std::exp( gsl_sf_lnchoose(k2,k-k1) + my_log_falling_factorial(k1+k2-k, k1) )<<std::endl;
			Rcpp::Rcout<<"coef giusto = "<< ( gsl_sf_fact(k1)*gsl_sf_fact(k2)  )/(gsl_sf_fact(k-k1)*gsl_sf_fact(k-k2)*gsl_sf_fact(k1+k2-k)  ) <<std::endl;
			log_vect_res[inner_indx] = gsl_sf_lnchoose(k2,k-k1) + my_log_falling_factorial(k1+k2-k, k1) + compute_Kprior_unnormalized(  k2, {n_i[n_i.size()-1]}, {gamma[gamma.size()-1]}  );
				
		 	// Check if it is the new maximum of log_vect_res
	       	if(log_vect_res[inner_indx]>val_max2){
	       		idx_max2 = inner_indx;
	       		val_max2 = log_vect_res[inner_indx];
	       	}
	   		inner_indx++;
		}
		// Update log_a:  log(a_i*alfa_i) = log(a_i) + log(alfa_i)
					Rcpp::Rcout<<"Vettore interno"<<std::endl;
					for(auto __v : log_vect_res)
						Rcpp::Rcout<<__v<<", ";
					Rcpp::Rcout<<std::endl;
		Rcpp::Rcout<<"Adding effect: log_stable_sum(log_vect_res) = "<<log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2)<<std::endl;			
		log_a[k1] += log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2);
		Rcpp::Rcout<<"Updated log_a["<<k1<<"] = "<<log_a[k1]<<std::endl;	      
		// Check if it is the new maximum of log_a
	   	if(log_a[k1]>val_max1){
	    	idx_max1 = k1;
	      	val_max1 = log_a[k1];
	    }   
	}

			Rcpp::Rcout<<"Vettore finale"<<std::endl;
			for(auto __v : log_a)
				Rcpp::Rcout<<__v<<", ";
			Rcpp::Rcout<<std::endl;
	// Complete the sum over all elements in log_a
	return log_stable_sum(log_a, TRUE, val_max1, idx_max1);
}

double p_distinct_prior_c_old(const unsigned int& k, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max  ){

	// Component prior preliminary operations
	ComponentPrior_Parameters qM_params;
	if(prior == "Poisson"){
		qM_params.Lambda = prior_param["lambda"];
	}
	else if(prior == "NegativeBinomial"){
		qM_params.p = prior_param["p"];
		qM_params.n_succ = prior_param["r"];
	}
	else{
		throw std::runtime_error("Error in p_distinct_prior_c, not implemented prior requested by R function");
	}

	//Rcpp::Rcout<<"Print ComponentPrior_Parameters: qM_params.Lambda = "<<qM_params.Lambda<<"; qM_params.p = "<<qM_params.p<<"; qM_params.n_succ = "<<qM_params.n_succ<<std::endl;
	auto qM_ptr = Select_ComponentPrior(prior, qM_params);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// Convert Rcpp vector
	std::vector<unsigned int> n_i   = Rcpp::as< std::vector<unsigned int> >(n_j);
	std::vector<double> gamma       = Rcpp::as< std::vector<double> >(gamma_j);

	// Compute normalization constant
			//Rcpp::Rcout<<"Calcolo log_V:"<<std::endl;
	//double V{ compute_Vprior(k, n_i, gamma, qM, M_max) }; 
	double log_V{ compute_log_Vprior(k, n_i, gamma, qM, M_max) };
			//Rcpp::Rcout<<"log_V = "<<log_V<<std::endl;

	// Compute unnormalized probability
			//Rcpp::Rcout<<"Calcolo log_K:"<<std::endl;
	double log_K{compute_Kprior_unnormalized(k, n_i, gamma)};
			//Rcpp::Rcout<<"log_K = "<<log_K<<std::endl;

	//return 
	return std::exp(log_V + log_K);
}

double p_distinct_prior_c(const unsigned int& k, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max  ){

	// Component prior preliminary operations
	ComponentPrior_Parameters qM_params;
	if(prior == "Poisson"){
		qM_params.Lambda = prior_param["lambda"];
	}
	else if(prior == "NegativeBinomial"){
		qM_params.p = prior_param["p"];
		qM_params.n_succ = prior_param["r"];
	}
	else{
		throw std::runtime_error("Error in p_distinct_prior_c, not implemented prior requested by R function");
	}

	//Rcpp::Rcout<<"Print ComponentPrior_Parameters: qM_params.Lambda = "<<qM_params.Lambda<<"; qM_params.p = "<<qM_params.p<<"; qM_params.n_succ = "<<qM_params.n_succ<<std::endl;
	auto qM_ptr = Select_ComponentPrior(prior, qM_params);
	ComponentPrior& qM(*qM_ptr);
	//Rcpp::Rcout<<"Selected prior is --> "<<qM.showMe()<<std::endl;

	// Convert Rcpp vector
	std::vector<unsigned int> n_i   = Rcpp::as< std::vector<unsigned int> >(n_j);
	std::vector<double> gamma       = Rcpp::as< std::vector<double> >(gamma_j);

	// Compute normalization constant
	Rcpp::Rcout<<"Calcolo log_V:"<<std::endl;
	double log_V{ compute_log_Vprior(k, n_i, gamma, qM, M_max) };
	Rcpp::Rcout<<"log_V = "<<log_V<<std::endl;

	// Compute unnormalized probability
	Rcpp::Rcout<<"Calcolo log_K:"<<std::endl;
	double log_K{compute_Kprior_unnormalized_recursive(k, n_i, gamma)};
	Rcpp::Rcout<<"log_K = "<<log_K<<std::endl;

	//return 
	return std::exp(log_V + log_K);
}


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


void Test_prod_sum(){

	double inf = std::numeric_limits<double>::infinity();

	std::vector<unsigned int> n{1,1,1};
	std::vector<double> gamma{1.0,2.0,3.0};
	unsigned int M{2};
	unsigned int k{0};

	std::vector<double> log_a{1.0, 2.0, 1.0};
	std::vector<double> a{ std::exp(1.0), std::exp(2.0), std::exp(1.0)};
	unsigned int idx{1};
	double max_log{log_a[idx]};
	double max_nolog{a[idx]};

	double res1 = combined_product(n,  gamma,  M,  k);
	double res2 = combined_sum(n,  gamma,  M,  k);

	// Test log_stable_sum
	double res3 = log_stable_sum(log_a, TRUE, max_log, idx);
	double res4 = log_stable_sum(a, FALSE, max_nolog, idx);

	double res5 = log_stable_sum(log_a, TRUE);
	double res6 = log_stable_sum(a, FALSE);

	Rcpp::Rcout<<"res combined_product = "<<res1<<std::endl;
	Rcpp::Rcout<<"res combined_sum = "<<res2<<std::endl;
	Rcpp::Rcout<<"res log_stable_sum (log and max) = "<<res3<<std::endl;
	Rcpp::Rcout<<"res log_stable_sum (no log and max) = "<<res4<<std::endl;
	Rcpp::Rcout<<"res log_stable_sum (log) = "<<res5<<std::endl;
	Rcpp::Rcout<<"res log_stable_sum ( no log) = "<<res6<<std::endl;

	Rcpp::Rcout<<"--------------------------------------------"<<std::endl;
	std::vector<double> a_zero{ 0.0, 0.0, 0.0 };
	std::vector<double> log_a_inf{ -inf, -inf, -inf };
	Rcpp::Rcout<<"res log_stable_sum( a_zero ) = "<<log_stable_sum(a_zero, FALSE, 0.0, 0)<<std::endl;
	Rcpp::Rcout<<"res log_stable_sum( log_a_inf ) = "<<log_stable_sum(log_a_inf, TRUE, -inf, 0)<<std::endl;
}