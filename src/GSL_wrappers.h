#ifndef __GSLWRAPPERS_H__
#define __GSLWRAPPERS_H__

#include "include_headers.h"
#include "recurrent_traits.h"

//GSL
#include <gsl/gsl_rng.h>     //For random number generators
#include <gsl/gsl_randist.h> //For random variates and probability density functions
#include <gsl/gsl_cdf.h> 	 //For cumulative density functions
#include <gsl/gsl_bspline.h> //For spline operations
#include <gsl/gsl_linalg.h> //For cholesky decomposition

//Load GDFMM traits
using namespace GDFMM_Traits;

/*
	Reference for random number generators:    https://www.gnu.org/software/gsl/doc/html/rng.html
	Reference for random number distributions: https://www.gnu.org/software/gsl/doc/html/randist.html#
*/

namespace sample{ //use the sample:: namespace to aviod clashes with R or oterh packages

	/*--------------------------------------------------------------
		Random number generator wrapper
	----------------------------------------------------------------*/

	//This class simply wraps in c++ code the construction and desctruction of a gsl_rng object.
	//I had to remove std::random_device because there is a bug when compiling in window (returs always the same value).
	// reference for bug -> https://en.cppreference.com/w/cpp/numeric/random/random_device, https://sourceforge.net/p/mingw-w64/bugs/338/
	class GSL_RNG{ 
		public:

			//constructor 1: takes one seed. If seed is 0, generates random seed
			GSL_RNG(unsigned int const & _seed){ 
				if(_seed == 0){
					seed = static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count());
					std::seed_seq seq = {seed}; //seed provived here has to be random. Than std::seed_seq adds entropy becasuse steady_clock is not sufficientyl widespread
					std::vector<unsigned int> seeds(1);
					seq.generate(seeds.begin(), seeds.end());
					seed = seeds[0];
					//std::cout<<"seed = "<<seed<<std::endl;
				}
				else{
					seed = _seed;
				}
				gsl_rng_env_setup();
				r = gsl_rng_alloc(gsl_rng_default);
				gsl_rng_set(r,seed);	
			}

			//constructor 0: default constructor. It is equvalent to the previous one using seed=0.
			GSL_RNG(){
				gsl_rng_env_setup();
				r = gsl_rng_alloc(gsl_rng_default);
				seed = static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count());
				std::seed_seq seq = {seed}; //seed provived here has to be random. Than std::seed_seq adds entropy becasuse steady_clock is not sufficientyl widespread
				std::vector<unsigned int> seeds(1);
				seq.generate(seeds.begin(), seeds.end());
				seed = seeds[0];
				gsl_rng_set(r,seed);
			}

			//descructor, can not be removed because gsl_rng_alloc() is allocating memory
			~GSL_RNG(){
				gsl_rng_free(r);
			}

			//print some information
			void print_info()const{
				printf ("generator type: %s\n", gsl_rng_name(r));
				std::cout<<"seed = "<<seed<<std::endl;
			}

			//call operator, return the engine ready to generate a number
			gsl_rng* operator()()const{
				return r;
			}

			//set seed, not tested much
			inline void set_seed(unsigned int const & s){
				seed = s;
				gsl_rng_set(r,seed);
			}

			//getter
			inline unsigned int get_seed() const{
				return seed;
			}
		private:
			gsl_rng * r; //the random number generator
			unsigned int seed; //the seed
	};

	/*--------------------------------------------------------------
		Random number distribution wrappers
	----------------------------------------------------------------*/
		
	//Callable object to draw a sample from sampling from Unif([0,1])
	struct runif
	{
		double operator()(GSL_RNG const & engine)const{
			return gsl_rng_uniform(engine()); //gsl_rng_uniform is a function, nothing has to be de-allocated
		}
		double operator()()const{
			return runif()(GSL_RNG ()); 
			/*GSL_RNG () creates a GSL_RNG obj calling the default constructor and than calls it call operator. 
			In other words, it generates a random number generator, generates a number and destroys the generator.
			It is equivalent to return gsl_rng_uniform( GSL_RNG ()() ).
			In this second case, one () is for the default constructor, the second () is for the call operator.  */
		}
	};

	//Callable object to draw a sample from Unif({0,...,N-1})
	struct runif_int 
	{
		unsigned int operator()(GSL_RNG const & engine, unsigned int const & N)const{
			return gsl_rng_uniform_int(engine(), N); //gsl_rng_uniform_int is a function, nothing has to be de-allocated.
		}
		unsigned int operator()(unsigned int const & N)const{
			return runif_int()(GSL_RNG (), N);
			/*GSL_RNG () creates a GSL_RNG obj calling the default constructor and than calls it call operator. 
			In other words, it generates a random number generator, generates a number and destroys the generator*/
		}
	};

	//Callable object to draw a sample from N(mean,sd). 
	// --> NB  it takes the standard deviation as input! <--
	struct rnorm
	{
		//Gets the engine
		//N(mean,sd)
		double operator()(GSL_RNG const & engine, double const & mean, double const & sd)const{
			return gsl_ran_gaussian_ziggurat(engine(),sd) + mean;
		}

		//Gets the engine
		//N(0,1)
		double operator()(GSL_RNG const & engine)const{
			return gsl_ran_gaussian_ziggurat(engine(), 1.0);
		}

		//Engine defaulted
		//N(mean,sd)
		double operator()(double const & mean, double const & sd)const{
			return rnorm()(GSL_RNG (), mean, sd);
		}

		//Engine defaulted
		//N(0,1)
		double operator()()const{
			return gsl_ran_gaussian_ziggurat(GSL_RNG ()(),1.0); //the first () is for the constructor, the second il for the call operator
		}
	};
	
	//Callable object to draw a sample from Gamma(shape,scale). 
	// --> NB  Watch out notation! gsl uses the scale parameter, in Europe we are used to the rate parameter (rate = 1/scale) <--	
	struct rgamma{

		//Gets the engine
		//Gamma(shape,scale)		
		double operator()(GSL_RNG const & engine, double const & shape, double const & scale)const{
			return gsl_ran_gamma(engine(),shape,scale);
		}

		//Engine defaulted
		//Gamma(shape,scale)		
		double operator()(double const & shape, double const & scale)const{
			return gsl_ran_gamma(GSL_RNG ()(),shape,scale);
		}
	};

	//Callable object to draw a sample from Chi-squared(k). 
	// --> NB  k is the degree of freedom. Chi-squared(k) = Gamma(shape = k/2, scale = 2) <--	
	struct rchisq{

		//Gets the engine
		//Chi-squared(k)		
		double operator()(GSL_RNG const & engine, double const & k)const{
			return gsl_ran_chisq(engine(),k);
		}

		//Engine defaulted
		//Chi-squared(k)
		double operator()(double const & k)const{
			return gsl_ran_chisq(GSL_RNG ()(), k);
		}
	};


	//Callable object to draw a sample from Dirichlet(alpha[0],...,alpha[K-1]). 
	//Both input and output Type are template parameters. Return type of a template function can not be template(not sure). For sure, it can be if the template parameter RetType describes the class
	//as in this case.
	// --> NB: Keep record of tested types! It is very challenging to check the type in the code. For example, alpha can not be a std::vector but I you do not get any compiler error if
	//		you insert a matrix. For sure, VecCol and VecRow works well, other types may be dangerous (not tested). <--
	template<typename RetType = VecCol> 
	struct rdirichlet{
		
		//Default constructor, used to check that Return
		rdirichlet(){
			//Check for Return Type
			static_assert( std::is_same_v<RetType, VecCol> || 
						   std::is_same_v<RetType, VecRow> ||
						   std::is_same_v<RetType, std::vector<double> >  ,
						  "______ ERROR, invalid Return Type requested in rdirichlet. Can handle only VecRow, VecCol and std::vector<double> _____");			
		}
		//Gets the engine
		//Dirichlet(alpha[0],...,alpha[K-1])
		template<typename Derived>
		RetType operator()(GSL_RNG const & engine, Eigen::MatrixBase<Derived> const & alpha)const{

			unsigned int dim = alpha.size();
			if(dim < 1){
				throw std::runtime_error("length of alpha (concentration parameter) in rdirichlet has to be positive");
			}

			//Declare return objec

			RetType return_obj;
			if constexpr(std::is_same_v<RetType, VecCol> || std::is_same_v<RetType, VecRow> ){

				return_obj = RetType::Zero(dim);				//Eigen obj that has to be returned
			}
			else if constexpr(std::is_same_v<RetType, std::vector<double> >){ //static_assert avoids any other possibility
				
				return_obj.resize(dim);
				return_obj = std::vector<double>(dim, 0.0);		//std::vector<double> that has to be returned
			}
			gsl_ran_dirichlet(engine(), dim, alpha.derived().data(), return_obj.data());
			//Return 
			return(return_obj);
		}

		//Engine defaulted
		//Dirichlet(alpha[0],...,alpha[K-1])		
		template<typename Derived>
		RetType operator()(Eigen::MatrixBase<Derived> const & alpha)const{
			return rdirichlet<RetType>()(GSL_RNG (), alpha);
		}

	};
}

#endif