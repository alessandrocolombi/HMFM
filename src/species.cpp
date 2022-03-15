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

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Factorials and Pochammer
//------------------------------------------------------------------------------------------------------------------------------------------------------


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


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	C numbers
//------------------------------------------------------------------------------------------------------------------------------------------------------


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

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	A priori functions
//------------------------------------------------------------------------------------------------------------------------------------------------------


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


// Works for d>2
double compute_log_Vprior(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max ){

	if(k == 0)
		throw std::runtime_error("Error in compute_log_Vprior. It does not make any sense to evaluate this function when k=0. The behaviuor is indefined");

	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) must be positive");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_log_Vprior, the length of n_i (group sizes) and gamma has to be equal");

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) 
		throw std::runtime_error("Error in compute_log_Vprior. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

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
		// Checks
		double _a = log_raising_factorial(k, Mstar+1 );	
		Rcpp::Rcout<<"log_raising_factorial = ("<< Mstar+1 <<")_("<< k <<") = "<<_a<<std::endl;		       					   			
		double _b = qM.log_eval_prob(Mstar + k);			
		Rcpp::Rcout<<"qM.log_eval_prob( "<< Mstar + k <<") = "<<_b<<std::endl;	

		Rcpp::Rcout<<"Ora calcolo inner_product:"<<std::endl;	       					   		
		Rcpp::Rcout<<"Stampo n_i: ";		
			for(auto __v : n_i)
				Rcpp::Rcout<<__v<<", ";
			Rcpp::Rcout<<std::endl;	
		Rcpp::Rcout<<"Stampo gamma: ";		
			for(auto __v : gamma)
				Rcpp::Rcout<<__v<<", ";
			Rcpp::Rcout<<std::endl;	
		Rcpp::Rcout<<"log[ (  "<<gamma[0]<<"("<< Mstar<<"+"<<k<<")  )^("<<n_i[0]<<") ] = "<<log_raising_factorial( n_i[0], gamma[0]*(Mstar + k) )<<std::endl;	
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

	//Rcpp::Rcout<<"Dentro a compute_Kprior_unnormalized: ";
	//for(auto __v : n_i)
		//Rcpp::Rcout<<__v<<", ";
	//Rcpp::Rcout<<std::endl;

	double inf = std::numeric_limits<double>::infinity();

	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Kprior, the length of n_i (group sizes) and gamma has to be equal");
	if(n_i.size() > 2 || n_i.size() == 0)
		throw std::runtime_error("Error in compute_Kprior, the length of n_i (group sizes) must be equal to 1 or 2");
	
	if(k == 0){
		//Need to handle degenerate cases. The probability of k=0 is 1 when n_j = 0. It is needed to have coherence when some elements of n_i are equal to 0.
		if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //if(n_i.size()==1 & n_i[0]==0)  // old version
			return 0.0;
		}	
		else
			return -inf;
	}
	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;



	if(n_i.size()==1){ // one group only 
		Rcpp::NumericVector absC = compute_logC(n_i[0], -gamma[0], 0.0); //absC[i] = |C(n1,i,-gamma1)| for i = 0,...,n1
		return absC[k];
	}
	else{ // two groups case

		// Compute all C numbers required
		Rcpp::NumericVector absC1 = compute_logC(n_i[0], -gamma[0], 0.0); //absC1[i] = |C(n1,i,-gamma1)| for i = 0,...,n1
		Rcpp::NumericVector absC2 = compute_logC(n_i[1], -gamma[1], 0.0); //absC2[i] = |C(n2,i,-gamma2)| for i = 0,...,n2

		const unsigned int start1 = std::max( 0, (int)k - (int)n_i[0] ); //max(0, k-n1)
		const unsigned int start2 = std::max( 0, (int)k - (int)n_i[1] ); //max(0, k-n2)
		const unsigned int end1   = std::min( (int)k, (int)n_i[1] );     //min(k,n2)
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
							//log_a[outer_indx] = gsl_sf_lnchoose(k,r1) - my_log_falling_factorial(r1,(double)k) + absC1[k-r1];  //<---- old version
			log_a[outer_indx] =  absC1[k-r1];
			
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
							//log_vect_res[inner_indx] = gsl_sf_lnchoose(k-r1,r2) + std::lgamma(k-r2+1) +  absC2[k-r2];   //<---- old version
				log_vect_res[inner_indx] = gsl_sf_lnchoose(k-r2,r1) + my_log_falling_factorial(k-r1-r2,(double)(k-r1)) +  absC2[k-r2];
				
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

//questa va bene per qualunque d, se d<=2 chiama compute_Kprior_unnormalized
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

				//Rcpp::Rcout<<"----> Dentro a k1 = "<<k1<<std::endl;
		// Compute recursive formula
		log_a[k1] = compute_Kprior_unnormalized_recursive(k1, {n_i.cbegin(), n_i.cend()-1} , {gamma.cbegin(), gamma.cend()-1} ); 
				//Rcpp::Rcout<<"log_a["<<k1<<"] = "<<"NNP(K^(2) = "<<k1<<") = "<<log_a[k1]<<std::endl;
		// Prepare for computing the second term
		
		// Initialize vector of results
		std::vector<double> log_vect_res(k1+1, -inf ); 

		// Initialize quantities to find the maximum of log_vect_res
		unsigned int idx_max2(0);
		double val_max2(log_vect_res[idx_max2]);
		
		// Inner loop on r2
		unsigned int inner_indx{0};
		for(std::size_t k2=k-k1; k2<= k; ++k2){
					//Rcpp::Rcout<<"***** Dentro a k2 = "<<k2<<std::endl;
			
			// Compute last element
					//Rcpp::Rcout<<"n_i[n_i.size()-1] = "<<n_i[n_i.size()-1]<<std::endl;
					//Rcpp::Rcout<<"gamma[gamma.size()-1] = "<<gamma[gamma.size()-1]<<std::endl;
					//Rcpp::Rcout<<"NNP(K^(1) = "<<k2<<") = "<< compute_Kprior_unnormalized(  k2, {n_i[n_i.size()-1]}, {gamma[gamma.size()-1]}  ) <<std::endl;
					//Rcpp::Rcout<<"coef = "<<std::exp( gsl_sf_lnchoose(k2,k-k1) + my_log_falling_factorial(k1+k2-k, k1) )<<std::endl;
					//Rcpp::Rcout<<"coef giusto = "<< ( gsl_sf_fact(k1)*gsl_sf_fact(k2)  )/(gsl_sf_fact(k-k1)*gsl_sf_fact(k-k2)*gsl_sf_fact(k1+k2-k)  ) <<std::endl;
			log_vect_res[inner_indx] = gsl_sf_lnchoose(k2,k-k1) + my_log_falling_factorial(k1+k2-k, k1) + compute_Kprior_unnormalized(  k2, {n_i[n_i.size()-1]}, {gamma[gamma.size()-1]}  );
				
		 	// Check if it is the new maximum of log_vect_res
	       	if(log_vect_res[inner_indx]>val_max2){
	       		idx_max2 = inner_indx;
	       		val_max2 = log_vect_res[inner_indx];
	       	}
	   		inner_indx++;
		}
		// Update log_a:  log(a_i*alfa_i) = log(a_i) + log(alfa_i)
					//Rcpp::Rcout<<"Vettore interno"<<std::endl;
					//for(auto __v : log_vect_res)
						//Rcpp::Rcout<<__v<<", ";
					//Rcpp::Rcout<<std::endl;
					//Rcpp::Rcout<<"Adding effect: log_stable_sum(log_vect_res) = "<<log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2)<<std::endl;			
		log_a[k1] += log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2);
					//Rcpp::Rcout<<"Updated log_a["<<k1<<"] = "<<log_a[k1]<<std::endl;	      
		// Check if it is the new maximum of log_a
	   	if(log_a[k1]>val_max1){
	    	idx_max1 = k1;
	      	val_max1 = log_a[k1];
	    }   
	}

			//Rcpp::Rcout<<"Vettore finale"<<std::endl;
			//for(auto __v : log_a)
				//Rcpp::Rcout<<__v<<", ";
			//Rcpp::Rcout<<std::endl;
	// Complete the sum over all elements in log_a
	return log_stable_sum(log_a, TRUE, val_max1, idx_max1);
}

//Direct formula per d=2
double compute_SK_prior_unnormalized(const unsigned int& k, const unsigned int& s, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma)
{

			//Rcpp::Rcout<<"Dentro compute_SK_prior_unnormalized"<<std::endl;
	double inf = std::numeric_limits<double>::infinity();

	//Checks
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_SK_prior_unnormalized, the length of n_i (group sizes) and gamma has to be equal");
	if(n_i.size() > 2 || n_i.size() == 0)
		throw std::runtime_error("Error in compute_SK_prior_unnormalized, the length of n_i (group sizes) must be equal to 1 or 2");

	// d=1
	if(n_i.size() == 1)
		throw std::runtime_error("Error in compute_SK_prior_unnormalized, the length of n_i (group sizes) can not be zero");
	// k=0 case
	if(k==0){
		if( s==0 & *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //if(n_i.size()==1 & n_i[0]==0)  // old version
			throw std::runtime_error("Error in compute_SK_prior_unnormalized, k and s can not be zero even if all n_j are zero");
			return 0.0;
		}	
		else
			return -inf;
	}
	// k<s case, if here k>0
	if(k<s)
		return -inf;
	// k>n1+n2 case
	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;
	// s > min(n1,n2) case
	if( s > *std::min_element(n_i.cbegin(),n_i.cend()) )
		return -inf;

	//If here, need to compute the formula
	// Compute all C numbers required
	Rcpp::NumericVector absC1 = compute_logC(n_i[0], -gamma[0], 0.0); //absC1[i] = |C(n1,i,-gamma1)| for i = 0,...,n1
	Rcpp::NumericVector absC2 = compute_logC(n_i[1], -gamma[1], 0.0); //absC2[i] = |C(n2,i,-gamma2)| for i = 0,...,n2

	//Compute and check the range of r1
	const unsigned int start1 = std::max( 0, (int)k - (int)n_i[1] ); //max(0,k-n2)
	const unsigned int end1   = std::min( (int)(k-s), (int)n_i[0] - (int)s ); //min(k-s,n1-s)
	const int range_size = (int)end1-(int)start1+1;
				//Rcpp::Rcout<<"start1 = "<<start1<<std::endl; 
				//Rcpp::Rcout<<"end1   = "<<end1<<std::endl;
				//Rcpp::Rcout<<"Range size = "<<range_size<<std::endl; 

	if(range_size <= 0) //All Cnumbers would be 0
		return -inf;			 
	
	std::vector<double> log_a(range_size, -inf);    // This vector contains all the quantities that depend only on r1
	// Initialize quantities to find the maximum of log_a
	unsigned int idx_max1(0);
	double val_max1(log_a[idx_max1]);
	// Start for loop
	unsigned int outer_indx{0};
	for(std::size_t r1=start1; r1 <= end1; ++r1){
		// Compute a_r1 using its definition
		log_a[outer_indx] = gsl_sf_lnchoose(k-r1,s) + my_log_falling_factorial(s,(double)(s+r1)) + absC1[s+r1] + absC2[k-r1] ;			

		// Check if it is the new maximum of log_a
	       if(log_a[outer_indx]>val_max1){
	       	idx_max1 = outer_indx;
	       	val_max1 = log_a[outer_indx];
	       }
	       outer_indx++;
	}
	// Complete the sum over all elements in log_a
	return log_stable_sum(log_a, TRUE, val_max1, idx_max1);
	
	//return -1.0;

}

//Recursive formula for d>2
double compute_SK_prior_unnormalized_recursive(const unsigned int& k, const unsigned int& s, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma){

			//Rcpp::Rcout<<"Dentro compute_SK_prior_unnormalized_recursive"<<std::endl;
	double inf = std::numeric_limits<double>::infinity();

	//Checks
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_SK_prior_unnormalized_recursive, the length of n_i (group sizes) and gamma has to be equal");
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_SK_prior_unnormalized_recursive, the length of n_i (group sizes) must be equal to 1 or 2");
	// d=1
	if(n_i.size() == 1){
		throw std::runtime_error("Error in compute_SK_prior_unnormalized_recursive, the length of n_i (group sizes) can not be zero");
		//Se finisco qua, una possibilità è chiamare compute_Kprior_unnormalized
	}
	// d=2 case
	if(n_i.size() == 2)
		return compute_SK_prior_unnormalized(k,s,n_i,gamma);
	// k=0 case
	if(k==0){
		if( s==0 & *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ 
			Rcpp::Rcout<<"Sono in compute_SK_prior_unnormalized_recursive e sono finito nel caso P(K=0,S=0) con tutti i termini di n_i nulli"<<std::endl;
			return 0.0;
		}	
		else
			return -inf;
	}
	// k<s case, if here k>0
	if(k<s)
		return -inf;
	// k>n1+...+nd case
	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;
	// s > min(n1,...,nd) case
	if( s > *std::min_element(n_i.cbegin(),n_i.cend()) )
		return -inf;

	
	// if here, n_i.size()>2 in a non trivial case
	
	// Initialize quantities for sum on s1
	std::vector<double> log_s1(k-s+1, -inf);    
	// Initialize quantities to find the maximum of log_s1
	unsigned int idx_max0(0);
	double val_max0(log_s1[idx_max0]);

	// Start loop on s1
	for(std::size_t s1 = 0; s1<=k-s; ++s1){
				//Rcpp::Rcout<<"^^^^ Dentro s1 = "<<s1<<std::endl;
		// Compute coeff that depends only on s1
		log_s1[s1] = my_log_falling_factorial(s,(double)(s+s1));
		
		// Initialize quantities for inner sums
		const unsigned int k_s_s1 = k-s-s1;
		// Compute the log of the inner sum on h and t
		std::vector<double> log_a(k_s_s1+1, -inf);    
		// Initialize quantities to find the maximum of log_a
		unsigned int idx_max1(0);
		double val_max1(log_a[idx_max1]);

		// Start for loop on h
		for(std::size_t h=0; h <= k_s_s1; ++h){

					//Rcpp::Rcout<<"----> Dentro a h = "<<h<<std::endl;
			// Compute recursive formula
			log_a[h] = compute_SK_prior_unnormalized_recursive(k-h, s+s1, {n_i.cbegin(), n_i.cend()-1} , {gamma.cbegin(), gamma.cend()-1} ); 
					//Rcpp::Rcout<<"log_a["<<h<<"] = "<<"NNP(S^(2) = "<<s+s1<<", K^(2) = "<<k-h<<") = "<<log_a[h]<<std::endl;


			// Prepare for computing the second term
			const unsigned int k_s_s1_h = k-s-s1-h;
			// Initialize vector of results
			std::vector<double> log_vect_res(k_s_s1_h+1, -inf ); 

			// Initialize quantities to find the maximum of log_vect_res
			unsigned int idx_max2(0);
			double val_max2(log_vect_res[idx_max2]);
			
			// Inner loop on t
			unsigned int inner_indx{0};
			for(std::size_t t=0; t<= k_s_s1_h; ++t){
						//Rcpp::Rcout<<"***** Dentro a t = "<<t<<std::endl;
				
				// Compute last element
							//Rcpp::Rcout<<"********************"<<std::endl;
							//Rcpp::Rcout<<"Coef = "<<std::exp(my_log_falling_factorial(s,(double)(s+s1)) + gsl_sf_lnchoose(k_s_s1_h,t) + gsl_sf_lnchoose(s+h+t,s) + my_log_falling_factorial(t,(double)(h+t)))<<std::endl;
							//Rcpp::Rcout<<"Coef giusto = "<< ( gsl_sf_fact(s+s1)*gsl_sf_fact(k-s-s1-h)*gsl_sf_fact(s+h+t) )/(gsl_sf_fact(s)*gsl_sf_fact(s1)*gsl_sf_fact(h)*gsl_sf_fact(t)*gsl_sf_fact(k-s-s1-h-t)  ) <<std::endl;
							//Rcpp::Rcout<<"--------------------"<<std::endl;
				log_vect_res[inner_indx] = gsl_sf_lnchoose(k_s_s1_h,t) + gsl_sf_lnchoose(s+h+t,s) + my_log_falling_factorial(t,(double)(h+t)) +
										   compute_Kprior_unnormalized( s+h+t, {n_i[n_i.size()-1]}, {gamma[gamma.size()-1]}  );
					
			 	// Check if it is the new maximum of log_vect_res
		       	if(log_vect_res[inner_indx]>val_max2){
		       		idx_max2 = inner_indx;
		       		val_max2 = log_vect_res[inner_indx];
		       	}
		   		inner_indx++;
			}
			// Update outer vector
			log_a[h] += log_stable_sum(log_vect_res, TRUE, val_max2, idx_max2);
						//Rcpp::Rcout<<"Updated log_a["<<h<<"] = "<<log_a[h]<<std::endl;	      
			// Check if it is the new maximum of log_a
		   	if(log_a[h]>val_max1){
		    	idx_max1 = h;
		      	val_max1 = log_a[h];
		    }   
		}

		// End of loop on h. Update the sum on s1
		log_s1[s1] += log_stable_sum(log_a, TRUE, val_max1, idx_max1);

		// Check if it is the new maximum of log_s1
	   	if(log_s1[s1]>val_max0){
	    	idx_max0 = s1;
	      	val_max0 = log_s1[s1];
	    }   

	}


	// Complete the sum over all elements in log_s1
	return log_stable_sum(log_s1, TRUE, val_max0, idx_max0);

}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	A posteriori functions
//------------------------------------------------------------------------------------------------------------------------------------------------------

std::vector<double> build_log_qM_post(const unsigned int& k, const std::vector<unsigned int>& n_i, const std::vector<double>& gamma, 
									  const ComponentPrior& qM, unsigned int M_max )
{
	// checks !!
	if(n_i.size() == 0)
		throw std::runtime_error("Error in build_log_qM_post, the length of n_i (group sizes) must be positive");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in build_log_qM_post, the length of n_i (group sizes) and gamma has to be equal");


	if(k == 0){
		// The behaviour of such cases is indefined beacuse i am not sure if the need to be handled for coherence or not
		if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //if(n_i.size()==1 & n_i[0]==0)  // old version
			throw std::runtime_error("Error in build_log_qM_post. The behaviuor of this function when k=0 and n_i is exactly zero is indefined.");
		}	
		else
			throw std::runtime_error("Error in build_log_qM_post. It does not make any sense to evaluate this function when k=0 and n_i is not exactly zero. The behaviuor is indefined");
	}

	// Compute all terms from 0 to M_max+k
	const unsigned int end{M_max+k};
	// Initialize vector of results
	std::vector<double> log_vect_res(end+1, -std::numeric_limits<double>::infinity() );
	// Initialize quantities to find the maximum
	unsigned int idx_max{0};
	double val_max(log_vect_res[idx_max]);

	// Start the loop, let us compute all the elements
	for(std::size_t Mstar=0; Mstar <= end; ++Mstar){

		// Formula implementation
		log_vect_res[Mstar] = log_raising_factorial(k,Mstar+1 ) +
							  qM.log_eval_prob(Mstar + k) -
							  std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
			       					   			  [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return log_raising_factorial( nj, gamma_j*(Mstar + k) );}
			       					   			); 
			      					   			       					   		
		// Check if it is the new maximum
	    if(log_vect_res[Mstar]>val_max){
	       	idx_max = Mstar;
	       	val_max = log_vect_res[Mstar];
	    }

	}

	// Formula to compute the log of all the sums in a stable way
	double log_norm_const{log_stable_sum(log_vect_res, TRUE, val_max, idx_max)};
			//Rcpp::Rcout<<"Stampo log_vect_res: ";		
			//for(auto __v : log_vect_res)
				//Rcpp::Rcout<<__v<<", ";
			//Rcpp::Rcout<<std::endl;
		//
			//Rcpp::Rcout<<"log_norm_const:"<<std::endl<<log_norm_const<<std::endl;

	// Normalize all terms
	std::transform(log_vect_res.begin(), log_vect_res.end(), log_vect_res.begin(), [&log_norm_const](double& x){return x - log_norm_const;}); // MA NON HA MODIFICATO NIENTE!!

			//Rcpp::Rcout<<"POST NORMALIZZAZIONE"<<std::endl;
				//Rcpp::Rcout<<"Stampo log_vect_res: ";		
			//for(auto __v : log_vect_res)
				//Rcpp::Rcout<<__v<<", ";
			//Rcpp::Rcout<<std::endl;
	//return
	return log_vect_res;
}

double compute_log_Vpost(const unsigned int& r, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i, 
						 const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max )
{

	//Checks
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vpost, the length of n_i (group sizes) must be positive");
	if(n_i.size() != m_i.size())
		throw std::runtime_error("Error in compute_log_Vpost, the length of n_i (old group sizes) and m_i (new group sizes) has to be equal");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_log_Vpost, the length of n_i (group sizes) and gamma has to be equal");

	// Special cases
	if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ 
		Rcpp::Rcout<<"The compute_log_Vpost has been called but vector n_i of previous observations is made of all zeros. Call the compute_log_Vprior function instead."<<std::endl;
		return compute_log_Vprior(r, m_i, gamma, qM, M_max );
	}
	if(k == 0)
		throw std::runtime_error("Error in compute_log_Vpost. If here, n_i is not exactly zero but k is set to 0. This makes no sense, such behaviuor is indefined");

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) 
		throw std::runtime_error("Error in compute_log_Vpost. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	//if(r == 0) --> r may be 0 in posterior calculations
		
	std::vector<double> log_qM_post{ build_log_qM_post(k, n_i, gamma, qM, M_max) }; //should be a vector of size M_max+k+1
	//Rcpp::Rcout<<"Stampo log_qM_post: ";		
	//for(auto __v : log_qM_post)
		//Rcpp::Rcout<<__v<<", ";
	//Rcpp::Rcout<<std::endl;

	const unsigned int start = r;
	const unsigned int end   = M_max;

	// Initialize vector of results
	std::vector<double> log_vect_res(end-start+1, -std::numeric_limits<double>::infinity() );
	// Initialize quantities to find the maximum
	unsigned int idx_max{0};
	double val_max(log_vect_res[idx_max]);

	// Can not use inner_product with 3 inputs. Need to compute in advance the vector of the arguments of log_raising_factorial
	std::vector<double> fact_argument(n_i.size(), 0.0);
	// Start the loop, let us compute all the elements
	unsigned int indx{0};
	for(std::size_t Mstar=0; Mstar <= M_max; ++Mstar){
		//Rcpp::Rcout<<"Mstar = "<<Mstar<<std::endl;

		// Fill the elements of fact_argument
		std::transform(n_i.cbegin(),n_i.cend(),gamma.cbegin(),fact_argument.begin(),
					   [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return ( (double)nj + (double)(Mstar+k)*gamma_j );} );

		// Formula implementation
		if(Mstar >= r){
			
			log_vect_res[indx] = my_log_falling_factorial(r, (double)Mstar ) +
								  log_qM_post[Mstar + k] -
								  std::inner_product( m_i.cbegin(),m_i.cend(),fact_argument.cbegin(), 0.0, std::plus<>(),
				       					   			  [](const unsigned int& mj, const double& xj){return log_raising_factorial( mj, xj );}
				       					   			);

			// Check if it is the new maximum
			if(log_vect_res[indx]>val_max){
			    idx_max = indx;
			    val_max = log_vect_res[indx];
			}	    

			indx++;

		}


	}
			//Rcpp::Rcout<<"Stampo log_vect_res: ";		
			//for(auto __v : log_vect_res)
				//Rcpp::Rcout<<__v<<", ";
			//Rcpp::Rcout<<std::endl;

	// Formula to compute the log of all the sums in a stable way
	return log_stable_sum(log_vect_res, TRUE, val_max, idx_max);   
}	

double compute_log_Vpost_naive(const unsigned int& r, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i, 
						 	   const std::vector<double>& gamma, const ComponentPrior& qM, unsigned int M_max )
{
	//Checks
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_log_Vpost, the length of n_i (group sizes) must be positive");
	if(n_i.size() != m_i.size())
		throw std::runtime_error("Error in compute_log_Vpost, the length of n_i (old group sizes) and m_i (new group sizes) has to be equal");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_log_Vpost, the length of n_i (group sizes) and gamma has to be equal");

	// Special cases
	if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ 
		Rcpp::Rcout<<"The compute_log_Vpost has been called but vector n_i of previous observations is made of all zeros. Call the compute_log_Vprior function instead."<<std::endl;
		return compute_log_Vprior(r, m_i, gamma, qM, M_max );
	}
	if(k == 0)
		throw std::runtime_error("Error in compute_log_Vpost. If here, n_i is not exactly zero but k is set to 0. This makes no sense, such behaviuor is indefined");

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) 
		throw std::runtime_error("Error in compute_log_Vpost. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");

	//if(r == 0) --> r may be 0 in posterior calculations
		
	double log_Vprior{compute_log_Vprior(k, n_i, gamma, qM, 1000)};
	Rcpp::Rcout<<"log_Vprior:"<<std::endl<<log_Vprior<<std::endl; 
		// Initialize vector of results
	const unsigned int start = r;
	const unsigned int end   = M_max;
		
		std::vector<double> log_vect_res(end-start+1, -std::numeric_limits<double>::infinity() );
		// Initialize quantities to find the maximum
		unsigned int idx_max{0};
		double val_max(log_vect_res[idx_max]);

		// Start the loop, let us compute all the elements
		unsigned int indx{0};
		for(std::size_t Mstar=start; Mstar <= end; ++Mstar){
			//Rcpp::Rcout<<"Mstar = "<<Mstar<<std::endl;
			// Formula implementation
			double inner_prod_2{0.0};
			for(std::size_t j=0; j<n_i.size(); ++j){
				inner_prod_2 += log_raising_factorial( m_i[j], gamma[j]*(Mstar + k)+n_i[j] );
			}
			
			log_vect_res[indx] =  gsl_sf_lnfact(Mstar+k) - gsl_sf_lnfact(Mstar-r) +
								  qM.log_eval_prob(Mstar + k) -
								  inner_prod_2 - 
								  std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
				       					   			  [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return log_raising_factorial( nj, gamma_j*(Mstar + k) );}
				       					   			);
			
			// Checks
			/*
			double _a = gsl_sf_lnfact(Mstar+k) - gsl_sf_lnfact(Mstar-r);
			Rcpp::Rcout<<"log_raising_factorial = ("<< Mstar+1 <<")_("<< k <<") = "<<_a<<std::endl;		       					   			
			double _b = qM.log_eval_prob(Mstar + k);			
			Rcpp::Rcout<<"qM.log_eval_prob( "<< Mstar + k <<") = "<<_b<<std::endl;	

			Rcpp::Rcout<<"Ora calcolo inner_product:"<<std::endl;	       					   		
			Rcpp::Rcout<<"Stampo n_i: ";		
				for(auto __v : n_i)
					Rcpp::Rcout<<__v<<", ";
				Rcpp::Rcout<<std::endl;	
			Rcpp::Rcout<<"Stampo gamma: ";		
				for(auto __v : gamma)
					Rcpp::Rcout<<__v<<", ";
				Rcpp::Rcout<<std::endl;	
			Rcpp::Rcout<<"log[ (  "<<gamma[0]<<"("<< Mstar<<"+"<<k<<")  )^("<<n_i[0]<<") ] = "<<log_raising_factorial( n_i[0], gamma[0]*(Mstar + k) )<<std::endl;	
			double _c = std::inner_product( n_i.cbegin(),n_i.cend(),gamma.cbegin(), 0.0, std::plus<>(),
				       					   			  [&Mstar, &k](const unsigned int& nj, const double& gamma_j){return log_raising_factorial( nj, gamma_j*(Mstar + k) );}
				       					  );
			Rcpp::Rcout<<"inner_product = "<<_c<<std::endl;	
			log_vect_res[Mstar] = _a +_b -_c;		
			
			//Rcpp::Rcout<<" ---> calcolato log_vect_res[Mstar] = "<<log_vect_res[Mstar]<<std::endl;	 
			*/      					   			       					   		
			// Check if it is the new maximum
	        if(log_vect_res[indx]>val_max){
	        	idx_max = indx;
	        	val_max = log_vect_res[indx];
	        }
	        indx++;
		}
		// Formula to compute the log of all the sums in a stable way
		Rcpp::Rcout<<"log_stable_sum(log_vect_res, TRUE, val_max, idx_max):"<<std::endl<<log_stable_sum(log_vect_res, TRUE, val_max, idx_max)<<std::endl;
		return ( log_stable_sum(log_vect_res, TRUE, val_max, idx_max) - log_Vprior );
}	

//questa è sola per 1 o 2 gruppi
double compute_Kpost_unnormalized(const unsigned int& r, const unsigned int& k, const std::vector<unsigned int>& m_i, const std::vector<unsigned int>& n_i, 
						 		  const std::vector<double>& gamma)
{

	//Rcpp::Rcout<<"Dentro a compute_Kpost_unnormalized: ";
	//for(auto __v : n_i)
		//Rcpp::Rcout<<__v<<", ";
	//Rcpp::Rcout<<std::endl;

	double inf = std::numeric_limits<double>::infinity();

	//Checks
	if(n_i.size() == 0)
		throw std::runtime_error("Error in compute_Kpost_unnormalized, the length of n_i (group sizes) must be positive");
	if(n_i.size() != m_i.size())
		throw std::runtime_error("Error in compute_Kpost_unnormalized, the length of n_i (old group sizes) and m_i (new group sizes) has to be equal");
	if(n_i.size() != gamma.size())
		throw std::runtime_error("Error in compute_Kpost_unnormalized, the length of n_i (group sizes) and gamma has to be equal");
	

	// Special cases
	if( *std::max_element(n_i.cbegin(),n_i.cend()) == 0 ){ //check it is a prior case
		Rcpp::Rcout<<"The compute_Kpost_unnormalized has been called but vector n_i of previous observations is made of all zeros. Call the compute_Kprior_unnormalized function instead"<<std::endl;
		return compute_Kprior_unnormalized(r, m_i, gamma );
	}

	if( k > std::accumulate(n_i.cbegin(), n_i.cend(), 0.0)  ) // check that the known data are coherent
		throw std::runtime_error("Error in compute_Kpost_unnormalized. It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined");
	
	if( r > std::accumulate(m_i.cbegin(), m_i.cend(), 0.0)  ) //The probability of having more distinc values than observations must be zero
		return -inf;

	if( *std::max_element(m_i.cbegin(),m_i.cend()) == 0 ){
		    if( r == 0 ) 
		    	return 0.0;
		    else
		    	return -inf; //This case is handle just for completeness and to sake of clarity. However, the case n1=n2=0 and r>0 is handled in the previous if (r > n1+n2)
	}

	// if here, it is not possible that n1=n2=0. 
	if( (n_i.size()==1) || ( n_i.size()==2 && m_i[1] == 0 ) ){ // one group only or m2=0
		Rcpp::NumericVector absC = compute_logC(  m_i[0], -gamma[0], - ( (double)k*gamma[0] + (double)n_i[0] )  ); //absC[i] = |C(m1,i,-gamma1,-(k*gamma1 + n1))| for i = 0,...,m1
		return absC[r];
	}
	else if( n_i.size()==2 && m_i[0] == 0 ){ // m1=0
		Rcpp::NumericVector absC = compute_logC(  m_i[1], -gamma[1], - ( (double)k*gamma[1] + (double)n_i[1] )  ); //absC[i] = |C(m2,i,-gamma2,-(k*gamma2 + n2))| for i = 0,...,m2
		return absC[r];
	}
	else{ // two groups case. 

		// Compute all C numbers required
		Rcpp::NumericVector absC1 = compute_logC(  m_i[0], -gamma[0], - ( (double)k*gamma[0] + (double)n_i[0] )  ); //absC[i] = |C(m1,i,-gamma1,-(k*gamma1 + n1))| for i = 0,...,m1
		Rcpp::NumericVector absC2 = compute_logC(  m_i[1], -gamma[1], - ( (double)k*gamma[1] + (double)n_i[1] )  ); //absC[i] = |C(m2,i,-gamma2,-(k*gamma2 + n2))| for i = 0,...,m2

		const unsigned int start1 = std::max( 0, (int)r - (int)m_i[0] ); //max(0, r-m1)
		const unsigned int start2 = std::max( 0, (int)r - (int)m_i[1] ); //max(0, r-m2)
		const unsigned int end1   = std::min( (int)r, (int)m_i[1] );     //min(r,m2)
					//Rcpp::Rcout<<"start1 = "<<start1<<std::endl; 
					//Rcpp::Rcout<<"start2 = "<<start2<<std::endl; 
					//Rcpp::Rcout<<"end1   = "<<end1<<std::endl; 
  
		std::vector<double> log_a(end1-start1+1, -inf);    // This vector contains all the quantities that depend only on r1
		// Initialize quantities to find the maximum of log_a
		unsigned int idx_max1(0);
		double val_max1(log_a[idx_max1]);

		// Start for loop
		unsigned int outer_indx{0};
		for(std::size_t r1=start1; r1 <= end1; ++r1){

			// Compute a_r1 using its definition
			log_a[outer_indx] =  absC1[r-r1];
			
			// Prepare for computing the second term
			const unsigned int end2{r-r1};     								//r-r1, this is just for coherence of notation, no operation is needed
			// Initialize vector of results
			std::vector<double> log_vect_res(end2-start2+1, -inf );

			// Initialize quantities to find the maximum of log_vect_res
			unsigned int idx_max2(0);
			double val_max2(log_vect_res[idx_max2]);
		
			// Inner loop on r2
			unsigned int inner_indx{0};
			for(std::size_t r2=start2; r2<= end2; ++r2){

				// Compute 
				log_vect_res[inner_indx] = gsl_sf_lnchoose(r-r2,r1) + my_log_falling_factorial(r-r1-r2,(double)(r-r1)) +  absC2[r-r2];
				
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

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Rcpp call functions
//------------------------------------------------------------------------------------------------------------------------------------------------------


double p_distinct_prior_c_old(const unsigned int& k, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, 
							  const Rcpp::List& prior_param, unsigned int M_max  )
{

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

double p_distinct_prior_c(const unsigned int& k, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, 
						  const Rcpp::List& prior_param, unsigned int M_max  )
{

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

 
double p_shared_prior_c(const unsigned int& s, const Rcpp::NumericVector& n_j, const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, 
					 	const Rcpp::List& prior_param, unsigned int M_max  )
{

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

	// Need to compute the probability for each possible K and marginalize it out
	double res{0.0};
	const unsigned int Kmax{ std::accumulate(n_i.cbegin(), n_i.cend(), 0) }; 

	const unsigned int start{ std::max(1,(int)s) }; //max(1,s) 
	//start can not be 0! However, it is possible that s=0 and k should range form s to n1+n2
	//The problem is that when k=0, it does not make any sense to compute compute_log_Vprior(0), one gets none. Actually it rises an error at the moment.
	//Of course, P(K=0,S=0) is zero (at least, it is zero when n1+n2>0). This term is computed separately and should be 0 (or it can be omitted)
	if(s==0){ //can be omitted in this function 
		res += std::exp(compute_SK_prior_unnormalized_recursive(0, s, n_i, gamma));
	}
	for(unsigned int k=start; k<=Kmax; k++){
				//Rcpp::Rcout<<"++++ k = "<<k<<std::endl;
		// Compute normalization constant
				//Rcpp::Rcout<<"Calcolo log_V:"<<std::endl;
		double log_V{ compute_log_Vprior(k, n_i, gamma, qM, M_max) };
				//Rcpp::Rcout<<"log_V = "<<log_V<<std::endl;

		// Compute unnormalized probability
				//Rcpp::Rcout<<"Calcolo log_SK:"<<std::endl;
		double log_SK{compute_SK_prior_unnormalized_recursive(k, s, n_i, gamma)};
				//Rcpp::Rcout<<"log_SK = "<<log_SK<<std::endl;

		// Compute P(K=k,S=s) and cumulate 
		res += std::exp(log_V + log_SK);
				//Rcpp::Rcout<<"P(K="<<k<<", S="<<s<<") = "<<std::exp(log_V + log_SK)<<std::endl;
				//Rcpp::Rcout<<"res:"<<std::endl<<res<<std::endl;
	}

	//return 
	return res;
}

double p_distinct_posterior_c(const unsigned int& r, const unsigned int& k, const Rcpp::NumericVector& m_j, const Rcpp::NumericVector& n_j, 
						      const Rcpp::NumericVector& gamma_j, const Rcpp::String& prior, const Rcpp::List& prior_param, unsigned int M_max )
{
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
	std::vector<unsigned int> m_i   = Rcpp::as< std::vector<unsigned int> >(m_j);
	std::vector<unsigned int> n_i   = Rcpp::as< std::vector<unsigned int> >(n_j);
	std::vector<double> gamma       = Rcpp::as< std::vector<double> >(gamma_j);

	// Compute normalization constant
			Rcpp::Rcout<<"Calcolo log_V posteriori:"<<std::endl;
	double log_Vpost{ compute_log_Vpost(r, k, m_i, n_i, gamma, qM, M_max ) };
			Rcpp::Rcout<<"log_Vpost = "<<log_Vpost<<std::endl;

			Rcpp::Rcout<<"Calcolo NAIVE log_V posteriori:"<<std::endl;
	double log_Vpost_naive{ compute_log_Vpost_naive(r, k, m_i, n_i, gamma, qM, M_max ) };
			Rcpp::Rcout<<"log_Vpost_naive = "<<log_Vpost_naive<<std::endl;		
	// Compute unnormalized probability
			Rcpp::Rcout<<"Calcolo log_K posteriori:"<<std::endl;
	double log_Kpost{compute_Kpost_unnormalized(r, k, m_i, n_i, gamma)};
			Rcpp::Rcout<<"log_Kpost = "<<log_Kpost<<std::endl;

	//return 
	//return std::exp(log_Vpost + log_Kpost);
	return std::exp(log_Vpost_naive + log_Kpost);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Tests
//------------------------------------------------------------------------------------------------------------------------------------------------------


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