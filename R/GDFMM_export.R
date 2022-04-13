#### GDFMM FUNCTIONS ####

# Here are collected all functions needed/usefull to run GDFMM Gibbs Sampler
# and analyze its results

#' sum_check : check if elements in a vector sum to an expected value. Elements of the vector
#'             are supposed to be computed as expected_sum*p. If expectations are not respected
#'             return a corrected version of vec, according to p
#'
#' @param vec vector to integer value to be checked
#' @param expected_sum expected value for sum(vec)
#' @param p theoretical weights applied to divide expected_sum in vec's elements
#' @return corrected version of vec (if needed)
#' @export
sum_check <- function(vec, expected_sum, p){
  if( !dplyr::near( sum(p), 1) )
    stop("weights must sum to 1!")

  round_diff = vec - expected_sum*p # large if value in vec is OVER-estimate

  if(sum(vec) > expected_sum ){
    i_max = which.max(round_diff)
    vec[i_max] = vec[i_max] - 1
    vec = sum_check(vec, expected_sum, p)
  }
  if(sum(vec) < expected_sum){
    i_min = which.min(round_diff)
    vec[i_min] = vec[i_min] + 1
    vec = sum_check(vec, expected_sum, p)
  }

  return(vec)
}


#' point_density : compute the point estimate with credible intervals for density in
#'                 in a point for a gaussian distribution and
#'
#' @param x point where density is computed
#' @param mu_vec vector of mean for the gaussian distribution
#'               extracted from posterior distribution of mu
#' @param sigma_vec vector of standard deviation for the gaussian distribution
#' @param alpha level of confidence
#' @return c(inf, mean, sup) of point estimate of the density in x
#' @export
point_density <- function(x, mu_vec, sigma_vec, alpha = 0.05) {
  estimated_density = mean(dnorm(x, mu_vec, sigma_vec))
  err = sd(dnorm(x, mu_vec, sigma_vec))*qnorm(1-alpha/2)
  ret = c(estimated_density - err, estimated_density, estimated_density + err)
  return(ret)
}

#' dnorm_est : compute the point-wise estimate of the density for points in the given grid.
#'             mu_vec e sigma_vec are respectively the vector of means and standard deviation
#'             estimate for the gaussian distribution
#'
#' @param grid vector of point where I want to evaluate the density
#' @param mu_vec vector of mean for the gaussian distribution
#' @param sigma_vec vector of standard deviation for the gaussian distribution
#' @param alpha level of confidence
#' @return named matrix with (Inf, estimate, Sup)
#' @export
dnorm_est <- function(grid, mu_vec, sigma_vec, alpha = 0.05){
  n = length(grid)
  CI_vec = data.frame( "Inf." = numeric(n), "Est." = numeric(n), "Sup." = numeric(n))

  for(i in 1:n){
    CI_vec[i,] = point_density(grid[i], mu_vec, sigma_vec, alpha)
  }

  return(CI_vec)
}

#' dmix : density of a mixture of gaussian distributions
#'
#' @param x point or vector of values where I want to evaluate the density
#' @param w vector of weights of the
#' @param mu_vec vector of mean for the gaussian distribution
#' @param sigma_vec vector of standard deviation for the gaussian distribution
#' @return named matrix with (Inf, estimate, Sup)
#' @export
dmix <- function(x, w_j, mu_vec, sigma_vec){
  if( !dplyr::near(sum(w_j), 1) )
    stop("weigths don't sum to one")

  if(length(w_j) != length(mu_vec) || length(w_j) != length(mu_vec) )
    stop("length of w_j, mu_vec and sigma_vec differs")

  K = length(w_j)

  n = length(x)
  mix_density = numeric(n)

  for(k in 1:K){
    mix_density = mix_density + w_j[k]*dnorm(x, mu_vec[k], sigma_vec[k])
  }

  return(mix_density)
}


#' rmix : function to extract a random sample of dimension n from a Gaussian mixture
#'        defined by the 3 vectors p, mu, sigma (that must have same length!)
#'
#' @param n dimension of the random sample that has to be extracted
#' @param p vector of weights for the components of the mixture
#' @param mu vector of means for the components
#' @param sigma vector of standard deviations for the components
#' @return random sample derived my the specified mixture
#' @export
rmix <- function(n, p, mu, sigma){
  if( !dplyr::near(sum(p),1) )
    stop("weights for the components must sum to 1")

  if( length(p) != length(mu) || length(p) != length(sigma) )
    stop("number of weights for the components must agree with the dimensions of
         mu's vector and sigma's vector AND *mu* and *sigma* must have
         length")
  # rnormm cannot take 0 weighted components ==> eliminate them
  ind_0 = which( p == 0 )

  if(length(ind_0) == 0){
    sample_p = p
    sample_mu = mu
    sample_sigma = sigma
  }else{
    sample_p = p[-ind_0]
    sample_mu = mu[-ind_0]
    sample_sigma = sigma[-ind_0]
  }

  n_p = sapply(sample_p, function(x){round(x*n)})
  n_p = sum_check(n_p, n, sample_p)
  n_p = as.vector(n_p, mode = "integer")

  mix_sample = numeric(n)
  i = 1
  M = length(sample_p)

  for(m in 1:M){
    ind_m = seq(i, i+n_p[m]-1)
    mix_sample[ind_m] = rnorm(n_p[m], sample_mu[m], sample_sigma[m])
    i = i + n_p[m]
  }

  # shuffle sampled values
  mix_sample = mix_sample[sample(n)]

  return( mix_sample )
}

#' simulate_data : function to simulate data from the specified configuration of Gaussian
#'                 mixtures in different groups. Results for each one of the \code{n_simul}
#'                 datasets are saved in a directory create with the time-stamp of the simulation.
#' @param n_simul number of datasets to be simulated
#' @param group_dim vector of numerosity for each group in a dataset
#' @param p_mix matrix with each row that specifies weights for the mixture in the corrispondent
#'              group
#' @param mu vector of means for the Gaussian components of the mixture
#' @param sigma vector of standard deviations for the Gaussian components of the mixture
#' @param burnin number of iterations to be discarded for the GDFMM
#' @param n_iter number of iterations to be saved for the GDFMM run
#' @param thin thinning of the GDFMM run
#' @param seed seed for the GDFMM run (0 ==> random seed)
#' @param option list of option for the GDFMM model. See GDFMM_sampler help for the list of
#'               needed values
#' @param dir path to the directory where data about simlation have to be saved
#' @return some metrics to evaluate the goodness of GDFMM
#' @export
simulate_data <- function(n_simul, group_dim, p_mix, mu, sigma,
                          burnin = 2000, n_iter = 2000, thin = 3, seed = 1234,
                          option, dir = ".")
  {
      # dimension check
  if( length(group_dim) != dim(p_mix)[1] )
    stop("*group_dim* must be a vector with same number of elements as number of *p_mix* rows")

  if( dim(p_mix)[2] != length(mu) || dim(p_mix)[2] != length(sigma) )
    stop("*p_mix* must have number of columns equal to the length of mu and sigma (vectors of same length)")

  # create a directory to save results and set it as working directory
  original_dir = getwd()
  name_dir = paste("/simulation_GDFMM_", Sys.time(), sep = "")
  name_dir = stringr::str_replace_all(name_dir,":",".")
  name_dir = paste(dir, name_dir, sep = "")

  dir.create(name_dir)
  setwd(name_dir)

  # create a .txt file to store all informations of the simulation

  cat("SIMULATION options are : \n", file = "simulation_INFO.txt")
  cat("- Datasets simulated : ", n_simul, "\n", file = "simulation_INFO.txt", append = T)
  cat("- Numerosity of the groups : ", group_dim, "\n", file = "simulation_INFO.txt", append = T)
  cat("- Mixture specifics : mu -> (", mu, ")  sigma-> (", sigma, ") \n \n",
      file = "simulation_INFO.txt", append = TRUE)

  option_text = ""
  opt_names = names(option)
  for(name_i in 1:length(opt_names)){
    option_text = paste(option_text, "-", opt_names[name_i], "=",
                        option[[opt_names[name_i]]], "\n", sep = " ")
  }
  cat("GDFMM options are : \n", option_text, file = "simulation_INFO.txt", append = TRUE)
  cat("- n_iter = ", n_iter, " seed = ", seed, file = "simulation_INFO.txt", append = TRUE)

  # first set the seed
  set.seed(seed)

  # retrieve d: number of groups
  d = length(group_dim)
  # retrieve M: number of components
  M = dim(p_mix)[2]
  # retrieve max length of a group
  n_max = max(group_dim)

  # create vectors to store metrics of the simulation
  KLdiv_vec = numeric(n_simul)
  Mdiff_vec = numeric(n_simul)

  for(i in 1:n_simul){
    cat("\014")
    cat("Iteration ", i, " of ", n_simul, "\n")
      # groups creation

    # initialize data matrix
    data <- matrix(NA, nrow = d, ncol = n_max)
    group_list = list()
    # create the groups
    for(j in 1:d){
      n_j <- group_dim[j]
      group_list[[j]] = rmix(n_j, p_mix[j,], mu, sigma)
      data[j, 1:n_j] <- group_list[[j]]
    }

    # 1st run of GDFMM --> obj: estimate partition
    GS <- GDFMM_sampler(data, n_iter, burnin, thin, seed = seed, option = option)
    cat("First run of GDFMM done --> computing the partition \n")

    # compute best partion via binder-loss
    part_matrix <- GS$Partition
    sim_matrix <- psm(part_matrix)
    binder_dahl <- dlso(part_matrix, loss = 'binder', estimate = sim_matrix)

    est_part = as.vector(binder_dahl)
    cat("Partition computed \n")
    # compute the estimate number of clusters
    k_est = length( levels( as.factor(est_part) ) )

    option_fixed = c(option, list("partition" = est_part) )

    GS_fixed <- GDFMM_sampler(data, niter = n_iter, burnin = burnin, thin = thin, seed = seed,
                              FixPartition = T, option = option_fixed)

    # for each group compute mean of weights for the estimated components
    w_list = GS_fixed$w_jk
    w = matrix(nrow = d, ncol = k_est)

    for(j in 1:d){
      w[j, ] = rowMeans(w_list[[j]])
    }

    # compute mean of mu from GS
    mu_list = GS_fixed$mu
    mu_GS = c()
    for(k in 1:k_est){
      mu_GS = c(mu_GS, mean(mu_list[[k]]))
    }

    # # compute mean of sigma from GS
    sigma_list = GS_fixed$sigma
    sigma_GS = c()
    for(k in 1:k_est){
      sigma_GS = c(sigma_GS, sqrt( mean(sigma_list[[k]]) ) )
    }

    # create grid to evaluate density
    max_val = max(data, na.rm = T)
    min_val = min(data, na.rm = T)
    range = max_val - min_val
    grid = seq(min_val - 0.1*range, max_val + 0.1*range, length = 1000)

    # compute real mixture density and estimated one; then compare them with
    # Kullback-Leibler Divergence
    KLdiv = 0.0
    GS_mix = list()
    real_mix = list()
    for(j in 1:d){
      GS_mix[[j]] = dmix(grid, w[j,], mu_GS, sigma_GS)
      real_mix[[j]] = dmix(grid, p_mix[j,], mu, sigma)
      KLdiv = KLdiv + KLD(real_mix[[j]], GS_mix[[j]])$sum.KLD.px.py
    }


    # save results for the current iteration
    save(list = c("GS", "GS_fixed", "est_part", "group_list", "grid", "GS_mix", "real_mix"),
         file = paste("run", i,".rda", sep = "") )

    # store metrics for the currente iteratiom
    KLdiv_vec[i] = KLdiv/d
    Mdiff_vec[i] = abs(M - k_est)
  }
  cat("DONE")

  # save list of metrics vector
  metrics = list("KLD" = KLdiv_vec, "M_diff" = Mdiff_vec)
  save(metrics, file = "metrics_list.rda")
  save(p_mix, file = "p_mix.rda")

  # restore original working directory
  setwd(original_dir)

  # return divergence metrics obtained through the simulation
  return(list("KLD" = KLdiv_vec, "M_diff" = Mdiff_vec))
}


#' GDFMM Gibbs Sampler: function to run the GDFMM model. There is the possibility to fix
#'                      the partition, passing TRUE to FixPartition and specifying the
#'                      partion in the option. Default prior for P0 is an inverse gamma
#'
#' @param data input data
#' @param niter number of iterations
#' @param burnin burnin period
#' @param thin thinning value
#' @param seed seed for GSL random engine (0 ==> random seed)
#' @param P0.prior string with the prior to be used as P0
#' @param FixPartition TRUE if we want to fix the partition
#' @param option list with initial values, hyperparameters and other options.
#'               Value always needed are : "Mstar0", "Lambda0", "mu0", "nu0", "sigma0",
#'               "k0", "Adapt_MH_hyp1", "Adapt_MH_hyp2", "Adapt_MH_power_lim", "Adapt_MH_var0",
#'               "alpha_gamma", "beta_gamma", "alpha_lambda", "beta_lambda".
#'               If FixPartition = TRUE, also value for "partition" is needed
#' @return results of Gibbs Sampler
#' @export
GDFMM_sampler <- function(data, niter, burnin, thin, seed,
                            P0.prior = "Normal-InvGamma", FixPartition = F, option) {

  #Data check and pre-processing
  #--> handle here different types of input types. User you be able to pass the data the simplest possible ways. For example, this
  #    function should be able to handle both matrixes and data.frames (or others if needed).

  #Check number of iterations

  #Check P0.prior

  #Check options
      # if(FixPartition){ check che una partition sia passata}
  if( any(is.na(data)) )
    stop("There are nan in data")
  # This is just an example, of course you can save the c++ output and perform further operations in R
  return( GDFMM:::GDFMM_sampler_c(data, niter, burnin, thin, seed, P0.prior, FixPartition, option))
}



#' p_distinct_prior
#'
#' This function computes the a priori probability that the number of distinct species is equal to \code{k}.
#' @param k integer, the number of distinct species whose probability has to be computed.
#' @param n_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
p_distinct_prior = function(k,n_j, gamma, prior = "Poisson", ..., Max_iter = 100){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check trivial cases
  if(k<0)
    stop("Error, the number of distinct species k can not be negative")
  if(k==0)
    return (0)
  if(k > sum(n_j))
    return (0)

  # Compute non trivial cases
  return (  p_distinct_prior_c(k,n_j,gamma,prior,prior_params,Max_iter)  )
}


#' p_shared_prior
#'
#' This function computes the a priori probability that the number of shared species is equal to \code{s}.
#' @param s integer, the number of shared species whose probability has to be computed.
#' @param n_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
p_shared_prior = function(s,n_j, gamma, prior = "Poisson", ..., Max_iter = 100){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check trivial cases
  if(s<0)
    stop("Error, the number of shared species s can not be negative")
  if(s > min(n_j))
    return (0)



  # Compute non trivial cases
  return (  p_shared_prior_c(s,n_j,gamma,prior,prior_params,Max_iter)  )
}



#' p_distinct_posterior
#'
#' This function computes the probability a posteriori that the number of distinct species is equal to \code{r}.
#' @param r integer, the number of distinct species whose probability has to be computed.
#' @param k integer, the number of distinct species that have been observed in the sample of sizes given by \code{n_j}.
#' @param m_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param n_j a positive scalar or a vector of positive integers that defines the size of the groups that have been alread observed. It must be of the same size of \code{m_j}.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
p_distinct_posterior = function(r, k, m_j, n_j, gamma, prior = "Poisson", ..., Max_iter = 100){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if(length(n_j)!=length(m_j))
    stop("The length of m_j must be equal to the length of n_j")
  if( any(m_j<0) || any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check Special cases
  if(k < 0)
    stop("Error, the number of distinct species k that have been observed can not be negative")
  if(sum(n_j) == 0)
    cat("p_distinct_posterior has been called but vector n_i of previous observations is made of all zeros. Call the p_distinct_prior function instead. \n")
  if(k > sum(n_j))
    stop("It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined")
  if(k == 0)
    stop("The number of distinct species k that have been observed can not be 0. Such value makes no sense when sum(n_j)>0 and the behaviuor is left undefined when sum(n_j)=0. ")
  if(r > sum(m_j) )
    return (0)
  if( r == 0 && sum(m_j) == 0 ){
    return (1) # just for coherence
  }

  # Compute non trivial cases
  return (  p_distinct_posterior_c(r,k,m_j,n_j,gamma,prior,prior_params,Max_iter)  )
}

#' p_shared_posterior
#'
#' This function computes the posterior probability that the number of shared species is equal to \code{t}.
#' @param t integer, the number of shared species whose probability has to be computed.
#' @param k integer, the number of distinct species that have been observed in the sample of sizes given by \code{n_j}.
#' @param m_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param n_j a positive scalar or a vector of positive integers that defines the size of the groups that have been alread observed. It must be of the same size of \code{m_j}.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
p_shared_posterior = function(t, k, m_j, n_j, gamma, prior = "Poisson", ..., Max_iter = 100){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if(length(n_j)!=length(m_j))
    stop("The length of m_j must be equal to the length of n_j")
  if( any(m_j<0) || any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check Special cases
  if(k < 0)
    stop("Error, the number of distinct species k that have been observed can not be negative")
  if(k == 0)
    stop("The number of distinct species k that have been observed can not be 0. Such value makes no sense when sum(n_j)>0 and the behaviuor is left undefined when sum(n_j)=0. ")
  if(sum(n_j) == 0)
    cat("p_distinct_posterior has been called but vector n_i of previous observations is made of all zeros. Call the p_distinct_prior function instead. \n")
  if(k > sum(n_j))
    stop("It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined")
  if(t<0)
    stop("Error, the number of shared species t can not be negative")
  if(t > min(m_j))
    return (0)

  # Compute non trivial cases
  return (p_shared_posterior_c(t, k, m_j, n_j, gamma, prior, prior_params, Max_iter ) )
}



#' Expected_prior
#'
#' This function computes the expected value of the a priori probability of distinct or shared species.
#' @param n_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param type string, select \code{distinct} to compute the expected value of the number of distinct species or select \code{shared} to compute the expected value of the number of shared species.
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
Expected_prior = function(n_j, gamma, type, prior = "Poisson", ..., Max_iter = 100, tol = 10^-12){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")
  if(type != "distinct" & type != "shared")
    stop("type can only be distinct or shared")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check trivial cases
  if(sum(n_j) == 0)
    stop("All values in n_j are 0.")

  if(type == "shared" & length(n_j) < 2 )
    stop("At least two groups are needed to compute the expected values of shared components")


  # Compute non trivial cases
  return( Expected_prior_c(n_j, gamma, type, prior, prior_params, Max_iter, tol  )  )
}


#' Expected_posterior
#'
#' This function computes the expected value of the a posteriori probability of distinct or shared species.
#' @param k integer, the number of distinct species that have been observed in the sample of sizes given by \code{n_j}.
#' @param m_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param n_j an positive integer in the case of exchangeable data or a vector of size \code{d} in the case of partially exchangeable data.
#' @param gamma real valued, it must be of the same size of \code{n_j}
#' @param type string, select \code{distinct} to compute the expected value of the number of distinct species or select \code{shared} to compute the expected value of the number of shared species.
#' @param prior a string that indicates the type of prior to be used for the number of components. It can only be equal to \code{"Poisson"} or \code{"NegativeBinomial"}.
#' @param ... the addition parameters to be used in the prior. Use \code{lambda} for the "Poisson"case (must be strictly positive) and \code{r} (positive integer) and \code{p} (real in (0,1)) for the "NegativeBinomial" case.
#'
#' @export
Expected_posterior = function(k, m_j, n_j, gamma, type, prior = "Poisson", ..., Max_iter = 100, tol = 10^-12){
  l = list(...)
  L = length(l)

  #checks
  if(length(n_j)!=length(gamma))
    stop("The length of n_j must be equal to the length of gamma")
  if( any(n_j<0) || any(gamma<=0) )
    stop("The elements of n_j must the non negative and the elements of gamma must be strictly positive")
  if(Max_iter<=0)
    stop("The number of iterations must be strictly positive")
  if(type != "distinct" & type != "shared")
    stop("type can only be distinct or shared")

  # read prior parameters
  prior_params = list("lambda" = -1, "r" = -1, "p" = -1)
  if(prior == "Poisson"){
    if(L!=1)
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter expected ")
    if(! names(l)=="lambda")
      stop("Error when reading the prior parameters: when prior is Poisson, only one parameter named lambda is expected. The name must be passed explicitely ")

    prior_params$lambda = l$lambda
  }
  else if(prior == "NegativeBinomial"){
    if(L!=2)
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly two parameters expected ")
    if( !any(! names(l) %in% names(prior_params)) )  #check names
      stop("Error when reading the prior parameters: when prior is NegativeBinomial, exactly one parameters named r and p are expected. The names must be passed explicitely ")

    prior_params$r = l$r
    prior_params$p = l$p
  }
  else
    stop("prior can only be equal to Poisson or NegativeBinomial")

  # Check trivial cases
  if(sum(n_j) == 0){
    stop("All values in n_j are 0. Call the Expected_prior function instead")
  }
  if(sum(m_j) == 0)
    stop("All values in m_j are 0.")

  if(type == "shared" & length(m_j) < 2 )
    stop("At least two groups are needed to compute the expected values of shared components")

  if(k < 0)
    stop("Error, the number of distinct species k that have been observed can not be negative")
  if(k > sum(n_j))
    stop("It is not possible that k is higher than the sum of the elements of n_i. The behaviuor is indefined")
  if(k == 0)
    stop("The number of distinct species k that have been observed can not be 0. Such value makes no sense when sum(n_j)>0 and the behaviuor is left undefined when sum(n_j)=0. ")

  # Compute non trivial cases
  return( Expected_posterior_c(k, m_j, n_j, gamma, type, prior, prior_params, Max_iter, tol)  )
}




#' genera_mix_gas
#'
#' @export
genera_mix_gas <- function(n = 200, pro=c(0.5,0.5), means = c(-1,1), sds=sqrt(c(1,1))){
    p <- length(pro)
    if(length(means)!=p){stop("The number of component do not coincides with the number of means components")}
    if(length(sds)!=p){stop("The number of component do not coincides with the number of standard deviation components")}

    y <- vector(length=n)
    clu <- vector(length=n)

    for(i in 1:n){
      u <- runif(1)

      cmp <- sample(1:p,1,prob=pro)

      y[i]=rnorm(1,mean=means[cmp],sd=sds[cmp])
      clu[i] <- cmp
    }

    return(list(y=y,clu=clu))
}


#' pred_uninorm
#'
#' @export
pred_uninorm <- function(idx_group, grid, fit){

    n_iter <- length(fit$K) #number of iterations
    l_grid <- length(grid)
    MIX    <- matrix(0, nrow=n_iter, ncol=l_grid)

    # This loop computes the predictive
    for(it in 1:n_iter){

      # Get sampled values
      M_it <- fit$K[it] + fit$Mstar[it] # compute the number of components (allocated or not)
      S_it = fit$S[[it]][idx_group,]    # get (S_{j,1}^(it), ..., S_{j,M}^(it)), where j is idx_group and M is M_it
      T_it = sum(S_it)                  # compute the sum of the vector above
      mu_it   <- fit$mu[[it]]           # get the mean, (mu_{1}^{(it)}, ..., mu_{M}^{(it)})
      sig2_it <- fit$sigma[[it]]        # get the variances, (sigma^2_{1}^{(it)}, ..., sigma^2_{M}^{(it)})

      # XX is a l_grid x M_it matrix, it contains the Normal kernels evauated over the grid
      # XX[i,m] = Norm(grid[i] | mu_{m}^{(it)}, sigma^2_{m}^{(it)})
      XX = t(sapply(1:M_it, simplify = "matrix",
      function(m){
      dnorm( x = grid, mean=mu_it[m], sd=sqrt(sig2_it[m]) )
      }
      ))
      # XX <- matrix(ncol=l_grid,nrow=M_it)
      # for(m in 1:M_it){ XX[m,] <- dnorm(grid,mean=mu_it[m],sd=sqrt(sig2_it[m]))}

      # Compute predicted density at iteration it
      MIX[it,] <- (S_it/T_it) %*% XX
    }

    # Density estimation and credible bounds
    pred_est <- apply(MIX,2,quantile,prob=c(0.025,0.5,0.975))
    return(pred_est)
}

#' predictive
#'
#' Same as \code{\link{pred_uninorm}} before but using apply insted of for
#'
#' @param idx_group [integer] the index of the group
#' @param grid [vector] a grid where the normal kernel is evaluated
#' @param fit [list] the output of \code{\link{GDFMM_sampler}}
#'
#' @return [matrix] of size \code{n x length(grid)} containing the quantiles of level \code{0.025,0.5,0.975}
#' @export
predictive <- function(idx_group, grid, fit){

    n_iter <- length(fit$K) #number of iterations
    l_grid <- length(grid)
    #MIX    <- matrix(0, nrow=n_iter, ncol=l_grid)

    # MIX is a n_iter x l_grid matrix
    # This loop computes the predictive
    MIX = t(sapply(1:n_iter, simplify = "matrix",
                    function(it){
                      # Get sampled values
                      M_it <- fit$K[it] + fit$Mstar[it] # compute the number of components (allocated or not)
                      S_it = fit$S[[it]][idx_group,]    # get (S_{j,1}^(it), ..., S_{j,M}^(it)), where j is idx_group and M is M_it
                      T_it = sum(S_it)                  # compute the sum of the vector above
                      mu_it   <- fit$mu[[it]]           # get the mean, (mu_{1}^{(it)}, ..., mu_{M}^{(it)})
                      sig2_it <- fit$sigma[[it]]        # get the variances, (sigma^2_{1}^{(it)}, ..., sigma^2_{M}^{(it)})

                      # XX is a l_grid x M_it matrix, it contains the Normal kernels evauated over the grid
                      # XX[i,m] = Norm(grid[i] | mu_{m}^{(it)}, sigma^2_{m}^{(it)})
                      XX = t(sapply(1:M_it, simplify = "matrix",
                                    function(m){
                                      dnorm( x = grid, mean=mu_it[m], sd=sqrt(sig2_it[m]) ) # returns a vector of length equal to l_grid
                                    }
                                  ))
                      # Compute predicted density at iteration it
                      (S_it/T_it) %*% XX
                    }
                ))


    # Density estimation and credible bounds
    pred_est <- apply(MIX,2,quantile,prob=c(0.025,0.5,0.975))
    return(pred_est)
}

#' pred_uninorm
#'
#' This function can easily be parallelized!
#' @inheritParams predictive
#' @return [list] of length \code{d} where each element is the return object of \code{\link{predictive}}.
#' @export
predictive_all_groups <- function(grid, fit){
  d = nrow(fit$gamma)
  lapply(1:d, predictive, grid = grid, fit = fit)
}
