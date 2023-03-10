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

#' simulate_data (never checked)
#'
#' function to simulate data from the specified configuration of Gaussian
#' mixtures in different groups. Results for each one of the \code{n_simul}
#' datasets are saved in a directory create with the time-stamp of the simulation.
#' @param n_simul [integer] number of datasets to be simulated
#' @param group_dim [vector] of numerosity for each group in a dataset
#' @param p_mix [matrix] with each row that specifies weights for the mixture in the corrispondent
#'              group
#' @param mu [vector] of means for the Gaussian components of the mixture
#' @param sigma [vector] of standard deviations for the Gaussian components of the mixture
#' @param burnin [integer] number of iterations to be discarded for the GDFMM
#' @param n_iter [integer] number of iterations to be saved for the GDFMM run
#' @param thin thinning of the GDFMM run
#' @param seed [integer] seed for the GDFMM run (0 ==> random seed)
#' @param option [list] of option for the GDFMM model. See GDFMM_sampler help for the list of
#'               needed values
#' @param dir path to the directory where data about simlation have to be saved
#' @return some metrics to evaluate the goodness of GDFMM
#' @export
simulate_data_OLD_POLIMI <- function(n_simul, group_dim, p_mix, mu, sigma,
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



#' arrange_partition
#'
#' This function takes as input a partition and fix it according to the notation used to define partitions in the sampler.
#' @param partition [vector] the input partition.
#' @return [vector] containing the partition following the requirement of the sampler.
#' @export
arrange_partition = function(partition){

  num = sort(unique(partition)) #get passed indices
  idx = seq(0, length(num)-1, by=1) #create sequence with correct indices

  for(h in 1:length(num)){
    if(num[h]!=idx[h])
      partition[partition == num[h]] = idx[h]
  }

  # A possible way to do this operation using apply is
  #apply(as.array(partition), MARGIN = 1, FUN = function(x){
          #h = which(num == x)
          #if(num[h]!=idx[h]) return (idx[h])
          #else return (x)
  #})
  # but it looks more expensive to me
  return (partition)
}

#' set_options
#'
#' Use this function to set up the options for the conditional Gibbs sampler.
#' @param partition [vector] of length equal to the number of data. If \code{NULL}, all data points are put in the same cluster.
#' @param Mstar0 [integer] the initial value of non-allocated components
#' @param Lambda0 [double] the initial value for Lambda.
#' @param gamma0 [double] the initial value for the gamma parameters.
#' @param mu0 [double] the mean parameter in the prior of mu.
#' @param k0 [double] the parameter in the prior of mu.
#' @param sigma0 [double] the rate parameter in the prior of sigma.
#' @param nu0 [double] the shape parameter in the prior of sigma.
#' @param beta0 [double] the prior mean for the regression coefficients.
#' @param Sigma0 [double] the prior covariance matrix for the regression coefficients.
#' @param IncludeCovariates [bool] set \code{TRUE} if covariates are provided.
#' @param Adapt_MH_hyp1 [double] default is 0.7.
#' @param Adapt_MH_hyp2 [double] default is 0.234.
#' @param Adapt_MH_power_lim [double] default is 10.
#' @param Adapt_MH_var0 [double] default is 1.
#' @param proposal_Mstar [integer] must be strictly positive. The proposal distribution for Metropolis-Hasting update of Mstar is a discrete uniform
#'        distribution between \code{-proposal_Mstar} and \code{proposal_Mstar}. 0 value is escluded. 
#' @param alpha_gamma [double] the shape parameter in the prior of gamma.
#' @param beta_gamma [double] the rate parameter in the prior of gamma.
#' @param alpha_lambda [double] the shape parameter in the prior of lambda.
#' @param beta_lambda [double] the rate parameter in the prior of lambda.
#' @param UpdateU [bool] set \code{TRUE} if U must be updated. Set \code{FALSE} to fix it to a common value.
#' @param UpdateM [bool] set \code{TRUE} if Mstar must be updated. Set \code{FALSE} to fix it to a common value.
#' @param UpdateGamma [bool] set \code{TRUE} if gamma must be updated. Set \code{FALSE} to fix it to a common value.
#' @param UpdateS [bool] set \code{TRUE} if S must be updated. Set \code{FALSE} to fix it to a common value.
#' @param UpdateTau [bool] set \code{TRUE} if tau must be updated. Set \code{FALSE} to fix it to a common value.
#' @param UpdateLambda [bool] set \code{TRUE} if Lambda must be updated. Set \code{FALSE} to fix it to a common value.
#' @param UpdateBeta [bool] set \code{TRUE} if Beta must be updated. Set \code{FALSE} to fix it to a common value.
#'
#' @export
set_options = function( partition = NULL, Mstar0 = 2,
                        Lambda0 = 3, mu0 = 0, sigma0 = 1, gamma0 = 1,
                        beta0 = c(0,0), Sigma0 = 10*diag(2),
                        IncludeCovariates = FALSE,
                        Adapt_MH_hyp1 = 0.7,Adapt_MH_hyp2 = 0.234, Adapt_MH_power_lim = 10, Adapt_MH_var0=1,
                        proposal_Mstar = 1,
                        k0 = 1/10, nu0 = 10, alpha_gamma = 1, beta_gamma = 1, alpha_lambda = 1, beta_lambda = 1,
                        init_mean_cluster = NULL, init_var_cluster = NULL,
                        UpdateU = T, UpdateM = T, UpdateGamma = T, UpdateS = T, UpdateTau = T, UpdateLambda = T, UpdateBeta = F
)
{
  option<-list("Mstar0" = Mstar0, "Lambda0" = Lambda0, "mu0" = mu0,"sigma0"= sigma0, "gamma0" = gamma0,
                "beta0" = beta0, "Sigma0" = Sigma0, "IncludeCovariates" = IncludeCovariates,
               "Adapt_MH_hyp1"= Adapt_MH_hyp1,"Adapt_MH_hyp2"= Adapt_MH_hyp2, "Adapt_MH_power_lim"=Adapt_MH_power_lim, "Adapt_MH_var0"=Adapt_MH_var0,
               "proposal_Mstar" = proposal_Mstar,
               "k0"= k0, "nu0"=nu0, "alpha_gamma"=alpha_gamma,
               "beta_gamma"=beta_gamma, "alpha_lambda"=alpha_lambda, "beta_lambda"=beta_lambda,
               "init_mean_cluster" = init_mean_cluster, "init_var_cluster" = init_var_cluster,
               "UpdateU" = UpdateU, "UpdateM" = UpdateM, "UpdateGamma" = UpdateGamma, "UpdateS" = UpdateS,
               "UpdateTau" = UpdateTau, "UpdateLambda" = UpdateLambda, "UpdateBeta" = UpdateBeta, "partition" = partition
  )
  return (option)
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
#' @param option the output of \code{\link{set_options}} function
#' @return results of Gibbs Sampler
#' @export
GDFMM_sampler <- function(data, niter, burnin, thin, seed,
                            P0.prior = "Normal-InvGamma", FixPartition = F, option = NULL) {

  n = ncol(data)*nrow(data) - sum(is.na(data)) #get number of data points

  #Check option to be in the correct form
  option_temp = set_options(partition = NULL)
  if(is.null(option)) # no option, set default
    option = option_temp

  if(length(option) != length(option_temp))
    stop("option parameter is malformed. Its length is not the expected one. Use set_options() function to set it correctely.")
  if(!all(names(option) == names(option_temp) ))
    stop("option parameter is malformed. The names are not the expected ones. Use set_options() function to set them correctely.")

  #Check partiton
  if(is.null(option$partition)){ # set empty partition
    if(FixPartition)
        stop("If FixPartition is selected, a partition must be provided in option$partition")
    option$partition = rep(0,n)
  }else{
    cat("\n Check that provided partition is well formed. It must start from 0 and all values must be contiguous \n")
    option$partition = arrange_partition(option$partition)

    # check that partiton and data are coherent
    if(n != length(option$partition))
      stop("The number of points in the data is not coherent with the length of the partition. Are there missing values in the data? Such implementation is not able to deal with them")
  }

  # Check initial values for tau
  K_init = length(table(option$partition)) # compute initial number of clusters
  if(is.null(option$init_var_cluster))
    option$init_var_cluster = rgamma(n=K_init+option$Mstar0,
                                     shape = option$nu0/2,
                                     rate  = option$nu0*option$sigma0/2 )
  if(length(option$init_var_cluster)!=K_init+option$Mstar0)
    stop("The length of option$init_var_cluster must be equal to the initial number of clusters deduced from the initial partition plus Mstar0 ")
  if(is.null(option$init_mean_cluster))
    option$init_mean_cluster = rnorm(n=K_init+option$Mstar0,
                                     option$mu0, sqrt(option$init_var_cluster/option$k0))
  if(length(option$init_mean_cluster)!=K_init+option$Mstar0)
    stop("The length of option$init_mean_cluster must be equal to the initial number of clusters deduced from the initial partition plus Mstar0")

  # Check proposal for Mstar
  option$proposal_Mstar = floor(option$proposal_Mstar)  
  if(option$proposal_Mstar <= 0)
    stop("proposal_Mstar must be a strictly positive integer")


  #if( any(is.na(data)) )
    #stop("There are nan in data") --> per come sto passando i dati non posso fare questo controllo. malissimo in ottica missing data

 # //sigma[m] =  1 / Gamma(gs_engine, nu0/2, 2 / (nu0*sigma0));
 #  //mu[m] = rnorm(gs_engine, mu0, sqrt(sigma[m]));


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
#' This function computes the expected value and the variance of the a priori probability of distinct or shared species.
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
#' This function computes the expected value and the variance of the a posteriori probability of distinct or shared species.
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
#' This function generate a random sample from a mixture of uni-dimensional Gaussian distributed random variables.
#'
#' @param n [integer] the number of points to be generated.
#' @param pro [vector] the mixture proportions. Its length must be equal to the one of \code{means} and \code{sds}.
#' @param means [vector] the centers of the mixture componentes. Its length must be equal to the one of \code{pro} and \code{sds}.
#' @param sds [vector] the standard deviations of the mixture componentes. Its length must be equal to the one of \code{pro} and \code{means}.
#' @return [list] with a vector named \code{y} containing the sampled values and a vector named \code{clu} with the cluster membership of each data point.
#' @export
genera_mix_gas <- function(n = 200, pro=c(0.5,0.5), means = c(-1,1), sds=sqrt(c(1,1))){
    p <- length(pro)
    if(length(means)!=p){stop("The number of component do not coincides with the number of means components")}
    if(length(sds)!=p){stop("The number of component do not coincides with the number of standard deviation components")}

    y <- vector(length=n)
    clu <- vector(length=n)

    for(i in 1:n){
      cmp <- sample(1:p,1,prob=pro)

      y[i]=rnorm(1,mean=means[cmp],sd=sds[cmp])
      clu[i] <- cmp
    }

    return(list(y=y,clu=clu))
}


#' pred_uninorm: this is the old version of Raf
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
#' This function computes the predictive distribution for group \code{idx_group} generated from the \code{\link{GDFMM_sampler}} or \code{\link{ConditionalSampler}}
#' @param idx_group [integer] the index of the group of interest.
#' @param grid [vector] a grid where the normal kernel is evaluated.
#' @param fit [list] the output of a conditional sampler, \code{\link{GDFMM_sampler}} or \code{\link{ConditionalSampler}}
#' @param burnin [integer] the number of draws from \code{\link{GDFMM_sampler}} that must be discarded.
#'
#' @return [matrix] of size \code{n x length(grid)} containing the quantiles of level \code{0.025,0.5,0.975}.
#' @export
predictive <- function(idx_group, grid, fit, burnin = 1){
    n_iter <- length(fit$mu) #number of iterations
    l_grid <- length(grid)
                            #MIX    <- matrix(0, nrow=n_iter, ncol=l_grid)
                            # MIX is a n_iter x l_grid matrix

    # This loop computes the predictive
    MIX = t(sapply(burnin:n_iter, simplify = "matrix",
                    function(it){
                      # Get sampled values
                      M_it  = length(fit$mu[[it]]) # compute the number of components (allocated or not)
                      #M_it <- fit$K[it] + fit$Mstar[it] # compute the number of components (allocated or not)
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

#' predictive_all_groups
#'
#' This function computes the predictive distribution for all groups generated from the \code{\link{GDFMM_sampler}}.
#' @inheritParams predictive
#' @return [list] of length \code{d} where each element is the return object of \code{\link{predictive}}.
#' @export
predictive_all_groups <- function(grid, fit, burnin = 1){
  d = nrow(fit$gamma)
  lapply(1:d, predictive, grid = grid, fit = fit, burnin = burnin)
}



#' generate_data
#'
#' VECCHIA, NON è DA USARE
#' This function generate data to be used in the \code{\link{GDFMM_sampler}} or \code{\link{GDFMM_marginal_sampler}}.
#' A major drawback of this function is that the current formulation allows only to use equal mixture proportions.
#' This is indeed the old version, the newest function is \code{\link{generate_data_prob}}.
#' @param d [integer] the number of levels in the data.
#' @param K [integer] the total number of clusters.
#' @param mu [vector] of length \code{K} defining the means of the clusters.
#' @param sd [vector] of length \code{K} defining the standard deviations of the clusters.
#' @param n_j [vector] of length \code{d} defining the cardinality of each groups.
#' @param seed [integer] the random seed.
#' @return [list] with a matrix of size \code{d x max(n_j)} named \code{data} containing the data to be fed to \code{\link{GDFMM_sampler}}
#' and a vector named \code{real_partition} with the cluster membership of each data point.
#' @export
generate_data <- function(d, K = 3, mu= c(-20,0,20), sd = c(1,1,1), n_j = rep(200, d), seed = 124123 )
{
  warning("VECCHIA, NON è DA USARE")
  if(length(mu) != K || length(sd) != K ) stop("The length of mu and sd must be equal to K")
  if(length(n_j) != d ) stop("The length of n_j must be equal to d")

  n = sum(n_j) #total number of data points
  p = matrix(0, nrow = d, ncol = K) # matrix with components weights

  set.seed(seed)
  Kgruppo = c() # used to save the number of generated clusters in each level
  componenti_gruppo = NULL # used to state what components are used to generate data in each level

  data = matrix(NA, nrow = d, ncol = max(n_j))     # d x max(n_j) matrix
  cluster = matrix(NA, nrow = d, ncol = max(n_j))  # d x max(n_j) matrix
  real_partition = c()      # real_partition is a vector of length sum(n_j), it collects all the group membership.
  # values are collected level by level, so first all the values in level 1, the all values in level 2 and so on
  # cluster label must always start from 0!
  for(j in 1:d){
    Kgruppo[j] = sample(1:K,1) # number of clusters in each level
    componenti_gruppo[[j]] = sample(1:K,Kgruppo[j], replace = F) # choose the components
    p[j,1:Kgruppo[j]] = rep(1/Kgruppo[j], Kgruppo[j]) # set the weights all equals
    appoggio = genera_mix_gas(n = n_j[j], pro = p[j,1:Kgruppo[j]], means = mu[ componenti_gruppo[[j]] ],
                              sds = sd[ componenti_gruppo[[j]] ] )

    data[j, 1:n_j[j]] = appoggio$y
    #cluster[j, 1:n_j[j]] = appoggio$clu, #errore, genera_mix_gas usa sempre indici che partono da 1!
    cluster[j, 1:n_j[j]] = unlist(lapply(1:n_j[j], function(h){componenti_gruppo[[j]][appoggio$clu[h]]}))
    real_partition = c(real_partition, cluster[j, 1:n_j[j]])
  }
  # In real partition devo avere valore da 0 a K senza buchi. Per esempio, se ho solo 0 e 2 non va bene!
  # quella che viene modificata dentro il sampler è
  # partion_within_sampler = arrange_partition(real_partition)
  return( list(data = data, real_partition = real_partition) )
}

#' generate_data
#'
#' VECCHIA, NON è DA USARE
#' This function generate data to be used in the \code{\link{GDFMM_sampler}} or \code{\link{GDFMM_marginal_sampler}}.
#' It gets the number of levels, the data to be generated in each level and, within each level, it samples from the mixture defined by \code{K}, \code{mu}, \code{sd} with
#' weights defined in \code{prob}, which may vary in different levels
#' @inheritParams generate_data
#' @param prob [matrix] of dimension \code{d x K} where each row contains the weights for the mixture model in each level. Some may be zero.
#' @return [list] with a matrix of size \code{d x max(n_j)} named \code{data} containing the data to be fed to \code{\link{GDFMM_sampler}}
#' and a vector named \code{real_partition} with the cluster membership of each data point.
#' @export
generate_data_prob <- function(d, K = 3, p = prob, mu= c(-20,0,20), sd = c(1,1,1), n_j = rep(200, d), seed = 124123 )
{
  warning("VECCHIA, NON è DA USARE")
  if(length(mu) != K || length(sd) != K ) stop("The length of mu and sd must be equal to K")
  if(length(n_j) != d ) stop("The length of n_j must be equal to d")

  n = sum(n_j) #total number of data points
  #p = matrix(0, nrow = d, ncol = K) # matrix with components weights

  set.seed(seed)
  Kgruppo = c() # used to save the number of generated clusters in each level
  componenti_gruppo = NULL # used to state what components are used to generate data in each level

  data = matrix(NA, nrow = d, ncol = max(n_j))     # d x max(n_j) matrix
  cluster = matrix(NA, nrow = d, ncol = max(n_j))  # d x max(n_j) matrix
  real_partition = c()      # real_partition is a vector of length sum(n_j), it collects all the group membership.
  # values are collected level by level, so first all the values in level 1, the all values in level 2 and so on
  # cluster label must always start from 0!
  for(j in 1:d){
    Kgruppo[j] = 3 # number of clusters in each level
    componenti_gruppo[[j]] = 1:3 # choose the components
    p[j,1:Kgruppo[j]] = prob[j,]
    #p[j,1:Kgruppo[j]] = rep(1/Kgruppo[j], Kgruppo[j]) # set the weights all equals
    appoggio = genera_mix_gas(n = n_j[j], pro = p[j,1:Kgruppo[j]], means = mu,
                              sds = sd )

    data[j, 1:n_j[j]] = appoggio$y
    #cluster[j, 1:n_j[j]] = appoggio$clu, #errore, genera_mix_gas usa sempre indici che partono da 1!
    cluster[j, 1:n_j[j]] = unlist(lapply(1:n_j[j], function(h){componenti_gruppo[[j]][appoggio$clu[h]]}))
    real_partition = c(real_partition, cluster[j, 1:n_j[j]])
  }
  # In real partition devo avere valore da 0 a K senza buchi. Per esempio, se ho solo 0 e 2 non va bene!
  # quella che viene modificata dentro il sampler è
  # partion_within_sampler = arrange_partition(real_partition)
  return( list(data = data, real_partition = real_partition) )
}

#' simulate_data
#'
#' This function generate data to be used in the \code{\link{GDFMM_sampler}} or \code{\link{GDFMM_marginal_sampler}}.
#' It gets the number of levels, the data to be generated in each level and, within each level, it samples from the mixture defined by \code{K}, \code{mu}, \code{sd} with
#' weights defined in \code{prob}, which may vary in different levels
#' @inheritParams generate_data
#' @param prob [matrix] of dimension \code{d x K} where each row contains the weights for the mixture model in each level. Some may be zero.
#' @return [list] with a matrix of size \code{d x max(n_j)} named \code{data} containing the data to be fed to \code{\link{GDFMM_sampler}}
#' and a vector named \code{real_partition} with the cluster membership of each data point.
#' @export
simulate_data <- function(d, K = 3, prob=NULL, mu= c(-20,0,20), sd = c(1,1,1), n_j = rep(200, d), seed = 124123 )
{

  set.seed(seed)
  if(length(mu) != K || length(sd) != K ) stop("The length of mu and sd must be equal to K")
  if(length(n_j) != d ) stop("The length of n_j must be equal to d")
  if(is.null(prob)){ # set equal weights
    prob = matrix(1/K,nrow = d,ncol= K)
  }

  n = sum(n_j) #total number of data points

  # used to save the number of generated clusters in each level
  #Kgruppo = apply(prob, MARGIN = 1, FUN = function(level_probs){sum(level_probs>0)})
  # used to state what components are used to generate data in each level
  #componenti_gruppo = vector("list",length = d)

  data = matrix(NA, nrow = d, ncol = max(n_j))     # d x max(n_j) matrix
  #cluster = matrix(NA, nrow = d, ncol = max(n_j))  # d x max(n_j) matrix
  real_partition = c()      # real_partition is a vector of length sum(n_j), it collects all the group membership.
  # values are collected level by level, so first all the values in level 1, the all values in level 2 and so on

  for(j in 1:d){

    #componenti_gruppo[[j]] = which(prob[j,]>0)
    #p[j,1:Kgruppo[j]] = prob[j,]
    #p[j,1:Kgruppo[j]] = rep(1/Kgruppo[j], Kgruppo[j]) # set the weights all equals

    # generate mixture in level j
    temp = genera_mix_gas(n = n_j[j], pro = prob[j,], means = mu, sds = sd )

    # save data
    data[j, 1:n_j[j]] = temp$y

    # save clustering
    #cluster[j, 1:n_j[j]] = temp$clu
    #cluster[j, 1:n_j[j]] = unlist(lapply(1:n_j[j], function(h){componenti_gruppo[[j]][temp$clu[h]]}))
    #real_partition = c(real_partition, cluster[j, 1:n_j[j]])
    real_partition = c(real_partition, temp$clu )
  }
  # In real partition devo avere valore da 0 a K senza buchi. Per esempio, se ho solo 0 e 2 non va bene!
  # quella che viene modificata dentro il sampler è
  # partion_within_sampler = arrange_partition(real_partition)
  return( list(data = data, real_partition = real_partition) )
}

#' predictive_new_group
#'
#' This function computes the predictive distribution for group \code{d+1} that has not been observed.
#' @inheritParams predictive
#' @inheritParams set_options
#'
#' @return [matrix] of size \code{n x length(grid)} containing the quantiles of level \code{0.025,0.5,0.975}.
#' @export
predictive_new_group <- function(grid, fit, burnin = 1, alpha_gamma, beta_gamma){
    n_iter <- length(fit$mu) #number of iterations
    l_grid <- length(grid)
                            #MIX    <- matrix(0, nrow=n_iter, ncol=l_grid)
                            # MIX is a n_iter x l_grid matrix

    # This loop computes the predictive
    MIX = t(sapply(burnin:n_iter, simplify = "matrix",
                    function(it){
                      # Get sampled values
                      M_it  = length(fit$mu[[it]]) # compute the number of components (allocated or not)
                      mu_it   <- fit$mu[[it]]           # get the mean, (mu_{1}^{(it)}, ..., mu_{M}^{(it)})
                      sig2_it <- fit$sigma[[it]]        # get the variances, (sigma^2_{1}^{(it)}, ..., sigma^2_{M}^{(it)})
                      gamma_new_it <- rgamma(n=1,shape=alpha_gamma,rate=beta_gamma) # draw gamma_d+1 from the prior
                      S_new_it <- rgamma(n=M_it, shape = gamma_new_it, rate = 1) # draw unnormalized weights from the prior
                      T_new_it <- sum(S_new_it) # needed to normalize the weigths

                      # Important remark. If alpha_gamma and beta_gamma are too small, it is possible that the sampled gamma is
                      # so close to zero that T is barely equal to 0.
                      if(T_new_it < 1e-8)
                        w_it = rep(0,M_it)
                      else
                        w_it = S_new_it/T_new_it  # get normalized weights
                      # XX is a l_grid x M_it matrix, it contains the Normal kernels evauated over the grid
                      # XX[i,m] = Norm(grid[i] | mu_{m}^{(it)}, sigma^2_{m}^{(it)})
                      XX = t(sapply(1:M_it, simplify = "matrix",
                                    function(m){
                                      dnorm( x = grid, mean=mu_it[m], sd=sqrt(sig2_it[m]) ) # returns a vector of length equal to l_grid
                                    }
                                  ))
                      if(any(is.na(XX))){
                        stop('Trovato NAN in XX')
                      }
                      # Compute predicted density at iteration it
                      w_it %*% XX
                    }
                ))


    # Density estimation and credible bounds
    pred_est <- apply(MIX,2,quantile,prob=c(0.025,0.5,0.975))
    return(pred_est)
}



#' set_options_marginal
#'
#' Use this function to set up the options for the conditional Gibbs sampler.
#' @param partition [vector] of length equal to the number of data. If \code{NULL}, all data points are put in the same cluster.
#' @param Lambda0 [double] the initial value for Lambda.
#' @param gamma0 [double] the initial value for the gamma parameters.
#' @param mu0 [double] the mean parameter in the prior of mu.
#' @param k0 [double] the parameter in the prior of mu.
#' @param sigma0 [double] the rate parameter in the prior of sigma.
#' @param nu0 [double] the shape parameter in the prior of sigma.
#' @param Adapt_MH_hyp1 [double] default is 0.7.
#' @param Adapt_MH_hyp2 [double] default is 0.234.
#' @param alpha_gamma [double] the shape parameter in the prior of gamma.
#' @param beta_gamma [double] the rate parameter in the prior of gamma.
#' @param alpha_lambda [double] the shape parameter in the prior of lambda.
#' @param beta_lambda [double] the rate parameter in the prior of lambda.
#' @param UpdateU [bool] set \code{TRUE} if U must be updated. Set \code{FALSE} to fix it to a common value.
#' @param UpdateGamma [bool] set \code{TRUE} if gamma must be updated. Set \code{FALSE} to fix it to a common value.
#' @param UpdateTau [bool] set \code{TRUE} if tau must be updated. Set \code{FALSE} to fix it to a common value.
#' @param UpdateLambda [bool] set \code{TRUE} if Lambda must be updated. Set \code{FALSE} to fix it to a common value.
#'
#' @export
set_options_marginal = function( partition = NULL,
                        Lambda0 = 3, mu0 = 0, sigma0 = 1, gamma0 = 1,
                        Adapt_MH_hyp1 = 0.7,Adapt_MH_hyp2 = 0.234,
                        sp_mala_U = 0.01, sp_mala_gamma=0.01,
                        k0 = 1/10, nu0 = 10, alpha_gamma = 1, beta_gamma = 1, alpha_lambda = 1, beta_lambda = 1,
                        init_mean_cluster = NULL, init_var_cluster = NULL,
                        UpdateU = T, UpdateGamma = T, UpdateTau = T, UpdateLambda = T
                      )
{
  option<-list("Lambda0" = Lambda0, "mu0" = mu0,"sigma0"= sigma0, "gamma0" = gamma0,
               "Adapt_MH_hyp1"= Adapt_MH_hyp1,"Adapt_MH_hyp2"= Adapt_MH_hyp2,
               "sp_mala_U"=sp_mala_U,"sp_mala_gamma"=sp_mala_gamma,
               "k0"= k0, "nu0"=nu0, "alpha_gamma"=alpha_gamma,
               "beta_gamma"=beta_gamma, "alpha_lambda"=alpha_lambda, "beta_lambda"=beta_lambda,
               "init_mean_cluster" = init_mean_cluster, "init_var_cluster" = init_var_cluster,
               "UpdateU" = UpdateU, "UpdateGamma" = UpdateGamma,
               "UpdateTau" = UpdateTau, "UpdateLambda" = UpdateLambda, "partition" = partition
              )
  return (option)
}

#' GDFMM Marginal Gibbs Sampler:
#'
#' function to run the GDFMM marginal Gibbs Sampler. There is the possibility to fix
#' the partition, passing TRUE to FixPartition and specifying the
#' partion in the option. Default prior for P0 is an inverse gamma
#'
#' @param data input data
#' @param niter number of iterations
#' @param burnin burnin period
#' @param thin thinning value
#' @param seed seed for GSL random engine (0 is for random seed)
#' @param P0.prior string with the prior to be used as P0
#' @param FixPartition TRUE if we want to fix the partition
#' @param option the output of \code{\link{set_options_marginal}} function
#' @return results of Gibbs Sampler
#' @export
GDFMM_marginal_sampler <- function( data, niter, burnin, thin, seed,
                                    P0.prior = "Normal-InvGamma", FixPartition = F, option = NULL)
{

  n = ncol(data)*nrow(data) - sum(is.na(data)) #get number of data points

  #Check option to be in the correct form
  option_temp = set_options_marginal(partition = NULL)
  if(is.null(option)) # no option, set default
    option = option_temp

  if(length(option) != length(option_temp))
    stop("option parameter is malformed. Its length is not the expected one. Use set_options() function to set it correctely.")
  if(!all(names(option) == names(option_temp) ))
    stop("option parameter is malformed. The names are not the expected ones. Use set_options() function to set them correctely.")

  #Check partiton
  if(is.null(option$partition)){ # set empty partition
    if(FixPartition)
        stop("If FixPartition is selected, a partition must be provided in option$partition")
    option$partition = rep(0,n)
  }else{
    cat("\n Check that provided partition is well formed. It must start from 0 and all values must be contiguous \n")
    option$partition = arrange_partition(option$partition)

    # check that partiton and data are coherent
    if(n != length(option$partition))
      stop("The number of points in the data is not coherent with the length of the partition. Are there missing values in the data? Such implementation is not able to deal with them")
  }

  # Check initial values for tau
  K_init = length(table(option$partition)) # compute initial number of clusters
  if(is.null(option$init_var_cluster))
    option$init_var_cluster = rgamma(n=K_init,
                                     shape = option$nu0/2,
                                     rate  = option$nu0*option$sigma0/2 )
  if(length(option$init_var_cluster)!=K_init)
    stop("The length of option$init_var_cluster must be equal to the initial number of clusters deduced from the initial partition ")
  if(is.null(option$init_mean_cluster))
    option$init_mean_cluster = rnorm(n=K_init,
                                     option$mu0, sqrt(option$init_var_cluster/option$k0))
  if(length(option$init_mean_cluster)!=K_init)
    stop("The length of option$init_mean_cluster must be equal to the initial number of clusters deduced from the initial partition ")

  #if( any(is.na(data)) )
    #stop("There are nan in data") --> per come sto passando i dati non posso fare questo controllo. malissimo in ottica missing data




  return( GDFMM:::GDFMM_marginal_sampler_c(data, niter, burnin, thin, seed, P0.prior, FixPartition, option))
}


#' empirical_bayes_normalinvgamma
#'
#' function to set normal-inversegamma parameters. If \code{data} is not \code{NULL}, an empirical bayes procedure
#' is automatically applied by setting \code{barmu = mean(data)}, \code{barsig2 = var(as.vector(data))/3}
#' \code{varmu} and \code{varsig2} are not set using \code{data}
#' @param data [matrix] input data.
#' @param barmu [scalar] desired value of the marginal mean of the normally distribuited parameter. Used only if \code{data} is \code{NULL}.
#' @param varmu [scalar] desired value of the marginal variance of the normally distribuited parameter
#' @param barsig2 [scalar] desired value of the mean of the inverse-gamma distribuited parameter. Used only if \code{data} is \code{NULL}.
#' @param varsig2 [scalar] desired value of the variance of the inverse-gamma distribuited parameter
#' @param correction [scalar] if \code{data} are provided, the expected value of \code{sigma} is set equal to \code{var(data)/correction}.
#' @export
empirical_bayes_normalinvgamma <- function( data = NULL, barmu = 0, varmu = 1, barsig2 = 10, varsig2 = 5, correction = 3  )
{
  if(!is.null(data)){
    barmu   = mean(data, na.rm = T)
    barsig2 = var(as.vector(data), na.rm = T) / correction
  }
  # Initialize return value
  res = list("mu0" = 0, "sigma0"= 1.0, "k0" = 1.0, "nu0" = 10 )

  # Compute useful quantities
  sum_barsig2_varsig2 = barsig2*barsig2 + varsig2
  sum_barsig2_2varsig2 = barsig2*barsig2 + 2*varsig2

  # Compute nu0
  res$nu0 = 2*( barsig2*barsig2/varsig2 + 2 )

  # Compute sigma0^2
  res$sigma0 = (sum_barsig2_varsig2/sum_barsig2_2varsig2) * barsig2

  # Compute mu0
  res$mu0 = barmu

  # Compute k0
  res$k0 = barsig2/varmu

  return(res)

}

#' data_mat2list
#'
#' function to transform the data matrix in long form in a list of vectors form.
#' @param data [matrix] input data in long form
#' @export
data_mat2list <- function( data )
{
  library(tidyverse)
  data = as_tibble(data) %>% mutate(level = as.factor(V1), index = as.integer(V2), value = as.double(V3)) %>%
         select(level,index,value)

  d = length(unique(data$level)) #get number of levels
  data_list = vector("list",length = d) #initialize list for data
  n_j = vector("numeric",length = d)
  for(j in 1:d){ #for each level
    data_nj = data %>% filter(level == j) # filter data in level j
    n_j[j] = nrow(data_nj) # compute number of data in each level
    data_list[[j]] = data_nj$value # fill the list of vectors
  }
  res = list("data" = data_list,
             "d" = d,
             "n_j" = n_j)

  return(res)
}


#' Non Central Student t - Density
#'
#' @param n0 the degree of fredoom.
#' @param mu0 the location parameter.
#' @param gamma0 the scale parameter.
#' @return scalar representing the evaluation of the density at x
#' @export
dnct = function( x, n0, mu0, gamma0 )
{
  if(gamma0 <= 0)
    stop("The scale parameter has to be strictly positive.")
  return( 1/gamma0 * dt(x = (x-mu0)/gamma0, df = n0 ) )
}


#' log_stable_sum
#'
#' This functions computes log(sum_i(a_i)) using a stable formula for log values. Let us denote a* to the the maximum value of vector a which is attained when i = i*.
#' \eqn{log(sum_i(a_i)) = log(a*) + log[ 1 + sum_{i not i*}(exp{log(a_i) - log(a*)}) ]}
#' See that only the logarithm of the elements of a are needed. Hence, it is likely that one has already computed them in log scale. If so, set is_log = T
#' @param a [vector] vector whose values must be summed.
#' @param is_log [bool] states if elements of a are already in log scale or not. If not, elements of a must be strictly positive.
#' @export
log_stable_sum = function(a, is_log = T)
{

  # Check vector length
  if(length(a) == 0)
    return (0.0)

  # Compute log of each element. They must be positve
  if(!is_log) {
    if(any(a<0))
      stop("Elements of a must be positive.")
    a = log(a)
  }

  # Find the maximum
  a_star = max(a)

  # Handle degenerate case
  if(a_star == -Inf)
    return (-Inf)

  # Compute log stable sum formula
  return(  a_star + log(sum( exp(a - a_star) ))   )
}


#' predictive_marginal
#'
#' This function computes the predictive distribution for group \code{idx_group} generated from the \code{\link{GDFMM_marginal_sampler}}.
#' @inheritParams predictive
#' @param option [list] the output of \code{\link{set_options_marginal}} function used to fit \code{\link{GDFMM_marginal_sampler}}.
#'
#' @return [matrix] of size \code{n x length(grid)} containing the quantiles of level \code{0.025,0.5,0.975}.
#' @export
predictive_marginal <- function(idx_group, grid, fit, option, burnin = 0)
{

  n_iter <- length(fit$K) #number of iterations
  n_save <- n_iter - burnin #number of saved iterations to take into account. The first burnin must be discarded
  l_grid <- length(grid)  #length of the grid
  MIX    <- matrix(0, nrow=n_save, ncol=l_grid)

  # This loop computes the predictive distribution over a grid
  for(it in (burnin + 1):n_iter){

     # Get sampled values
     M_it <- fit$K[it]                        # get the number of clusters (allocated or not)
     mu_it   <- fit$mu[[it]]                  # get the mean, (mu_{1}^{(it)}, ..., mu_{M}^{(it)})
     sig2_it <- fit$sigma[[it]]               # get the variances, (sigma^2_{1}^{(it)}, ..., sigma^2_{M}^{(it)})
     log_q_it <- fit$log_q[[it]][idx_group,]  # get the weights of required group

    if(length(log_q_it)!=(M_it+1))
      stop("Error in predictive_marginal, length of log_q_it must be M_it + 1")
    # Kernel_grid is a  M_it x l_grid matrix, it contains the Normal kernels evauated over the grid
    # Kernel_grid[m,i] = Norm(grid[i] | mu_{m}^{(it)}, sigma^2_{m}^{(it)})
    Kernel_grid = t(sapply(1:M_it, simplify = "matrix",
                     function(m){
                         dnorm( x = grid, mean=mu_it[m], sd=sqrt(sig2_it[m]) )
                     }
                   ))

    # Prior_grid is a vector of length l_grid, it contains the marginal prior evauated over the grid
    # Prior_grid[i] = nct(grid[i] | dof = nu0, loc = mu0, scale = sqrt(k0/(k0+1)*sigma0) ), where nct is the non-central student-t distribution
    scale = sqrt( (option$k0)/(option$k0 + 1) * option$sigma0 )
    Prior_grid = GDFMM:::dnct(x = grid, n0 = option$nu0, mu0 = option$mu0, gamma0 = scale)

    # Compute normalized weigths
    log_stable_sum = log_stable_sum(log_q_it)
    q = exp( log_q_it - log_stable_sum )
    # Compute predicted density at iteration it
    MIX[it-burnin,] <- q[M_it+1] * Prior_grid +  q[1:M_it] %*% Kernel_grid
  }
    # Density estimation and credible bounds
  pred_est <- apply(MIX,2,quantile,prob=c(0.025,0.5,0.975))
  return(pred_est)

}




#' predictive_marginal_all_groups
#'
#' This function computes the predictive distribution for all groups generated from the \code{\link{GDFMM_marginal_sampler}}.
#' @inheritParams predictive
#' @return [list] of length \code{d} where each element is the return object of \code{\link{predictive_marginal}}.
#' @export
predictive_marginal_all_groups <- function(grid, fit, option, burnin = 0){
  d = nrow(fit$gamma)
  lapply(1:d, predictive_marginal, grid = grid, fit = fit, option = option, burnin = burnin)
}







#' Compute_coclust_error
#'
#' This function computes the coclustering error and the coclustering error star.
#' Being errors, low values are to be preferred and 1 is the maximum, worst, value.
#' See Bassetti,Casarin et al. 2020 for proper definition
#' @param real_partition [vector] of length \code{n}, the total number of points, containing the true cluster membership of each data point
#' @param psm [matrix] of size \code{n x n} containing the posterior similarity matrix, i.e each element represents the relative frequency of two data points being assigned to the same cluster
#'
#' @return [list] \code{coclust_err} and \code{coclust_err_star}
#' @export
Compute_coclust_error = function(real_partition, psm){

  Ndata = length(real_partition)

  # Compute the true pairwise co-clustering matrix
  expg = expand.grid(real_partition, real_partition)
  expg$match = as.numeric(expg$Var1==expg$Var2)

  coclust_true = matrix( expg$match,
                         nrow = Ndata,
                         ncol = Ndata,
                         byrow = T
  )

  # Co-clustering error
  coclust_error = (1/(Ndata^2)) * sum(abs(coclust_true - psm))

  # Co-clustering error star
  psm[psm<0.5]  = 0
  psm[psm>=0.5] = 1
  coclust_error_star = (1/(Ndata^2)) * sum(abs(coclust_true - psm))

  return(  list("coclust_err"=coclust_error,
                "coclust_err_star"=coclust_error_star)  )

}


#' Compute_L1_dist
#'
#' This function computes the L1 distance between the true density and the predictive density within each level
#' and the mean value across all levels.
#' Low values are to be preferred.
#' @param Pred [list] the predictive distributions of all levels evaluated  \code{grid}
#' @inheritParams simulate_data
#' @inheritParams predictive
#'
#' @return [list] the first element, ??? , is a vector with the L1 error within each level.
#' The second element, ???, is the mean value across all levels.
#'
#' @export
#'
#' @examples Compute_L1_dist(list(Pred_all[[1]][2,], Pred_all[[2]][2,]),mix_probs,mu,sd)
Compute_L1_dist = function(Pred, p_mix, mu, sigma, grid ){

  d = length(Pred)
  K = length(mu)
  L1_err = rep(0,d)

  if(d != nrow(p_mix))
    stop("Mismatch of dimensions. Pred must be a list of d elements,
          p_mix must be a matrix with d rows")
  if(K != length(sigma) || K != ncol(p_mix))
    stop("Mismatch of dimensions. mu and sigma must be vector of length K.
          p_mix must be a matrix with K columns")

  # Evaluate the true densities over a grid of points
  true_dens_eval = apply(p_mix, 1, FUN = function(x){
    GDFMM::dmix(x = grid, w_j = x, mu_vec = mu, sigma_vec = sd)
  })


  for(j in 1:d){

    # Compute the absolute value of the differences between the true and the predicted densities
    ydiff = abs( true_dens_eval[,j] - Pred[[j]] )
    # Compute L1 distance
    L1_err[j] = pracma::trapz(x = grid, y = ydiff)

  }

  return(  list( "L1err_per_level" = L1_err,
                 "L1err_average" = mean(L1_err)
               )
         )

}


#' Handle input
#'
#' This function gets a tibble with three columns containing the data in long form. Data are handled and a list is returned.
#' @param tb [tibble] a tibble with three columns. First column contains the ID of the data, 
#' second column contains the level membership of each data. Finally, the third column contains the numerical value of the variabile
#' @export
input_handle = function(tb){
  
  # set names
  ncol_tb = ncol(tb)
  r = ncol_tb - 3 # number of covariates
  names(tb)[1:3] = c("ID","level","value")
  cov_names = names(tb)[4:ncol_tb]

  tb = tb %>% ungroup()
  IDs  = tb %>% distinct(ID) %>% pull(ID)
  # compute number of individuals in each level
  n_j = tb %>% distinct(ID,level) %>% 
               group_by(level) %>% 
               summarise(count = n()) %>% 
               pull(count)
  
  # compute number of levels
  d = length(n_j)
  
  # compute number of individuals
  n = tb %>% group_by(ID) %>% 
             summarise(n()) %>% 
             nrow()
  
  # compute quantities for each individual i in level j
  N_ji = matrix(0,nrow = d, ncol = n) # number of observations for individual i in level j
  mean_ji = matrix(0,nrow = d, ncol = n) # mean of observations for individual i in level j
  var_ji  = matrix(0,nrow = d, ncol = n) # variance of observations for individual i in level j
  #s_i  = rep(0,n) # number of levels in which individual i apprears

  data = vector("list", length = d) # list containing all observed values, for each level and for each individual
  data = lapply(1:d, FUN = function(s){data[[s]] = vector("list", length = n) })

  cov_list = vector("list", length = d) # list containing all covariates, for each level and for each individual
  cov_list = lapply(1:d, FUN = function(s){cov_list[[s]] = vector("list", length = n) })

  if(r > 0)
    formula <- as.formula( paste("value ~ ", paste(cov_names, collapse = "+")) )
  
  for(i in 1:n) {
    filter_i = tb %>% filter(ID == IDs[i]) # filter for i-th individual
    
              ### compute mean-var-number of observations for i-th individual in all d levels
              ##temp_i = filter_i %>% group_by(level) %>% 
                ##summarise(count = n(), mean = mean(value), var = var(value)) %>% 
                ##select(count, mean, var,level) %>% mutate(level = as.integer(level)) %>% as.matrix()
          ##    
              ##temp_i[is.na(temp_i[,3]),3] = 0
              ##s_i[i] = nrow(temp_i)
          ##    
              ### save
              ##N_ji[as.numeric(temp_i[,4]),i]    = as.numeric(temp_i[,1])
              ##mean_ji[as.numeric(temp_i[,4]),i] = as.numeric(temp_i[,2])
              ##var_ji[as.numeric(temp_i[,4]),i]  = as.numeric(temp_i[,3])
    
    # fill data structure 
    for(j in 1:d){
      temp_ji = filter_i %>% filter(level == j)
      data_ji = temp_ji %>% pull(value) # get all observation for individual i in all d levels
      data[[j]][[i]] = data_ji
      N_ji[j,i] = length(data_ji) # compute number of observations

      if(N_ji[j,i] > 0)
        mean_ji[j,i] = mean(data_ji) # compute mean
      
      # compute variance, if not defined set 0
      if(N_ji[j,i] > 1)
        var_ji[j,i] = var(data_ji)

      # get design matrix for observation i in level j
      if(r > 0){
        if(N_ji[j,i] == 0){
          cov_list[[j]][[i]] = NA
        }else if(N_ji[j,i] == 1){
          covariates = model.matrix( formula, data = temp_ji )[,-1]
          cov_list[[j]][[i]] = matrix(covariates, nrow = 1, ncol = r)
        }else{
          cov_list[[j]][[i]] = model.matrix( formula, data = temp_ji )[,-1]
        }
      }
      
      
    }  
  
  }
  if(r == 0)
    cov_list = NULL
  return( list("n"=n,
               "d"=d,
               "r"=r,
               "n_j"=n_j,
               "ID_i" = as.character(IDs),
               "observations"=data,#"s_i"=s_i,
               "covariates" = cov_list,
               "N_ji"=N_ji,
               "mean_ji"=mean_ji,
               "var_ji"=var_ji)
        )
  
}


#' Conditional Sampler: function to run the GDFMM model. There is the possibility to fix
#'                      the partition, passing TRUE to FixPartition and specifying the
#'                      partion in the option. Default prior for P0 is an inverse gamma
#'
#' @param data [tibble], the input data. This must be the return object of \code{\link{handle_input}}
#' @param niter [integer], the number of iterations
#' @param burnin [integer], the burnin period
#' @param thin [integer], the thinning value
#' @param seed [integer], the seed for GSL random engine (0 ==> random seed)
#' @param P0.prior [string] with the prior to be used as P0
#' @param FixPartition TRUE if we want to fix the partition
#' @param option [list] the output of \code{\link{set_options}} function
#' @export
ConditionalSampler <- function(data, niter, burnin, thin, seed,
                               P0.prior = "Normal-InvGamma", FixPartition = F, option = NULL) {

  # check input data
  names_data_input = c("n","d","r","n_j","ID_i","observations","covariates","N_ji","mean_ji","var_ji")
  if(length(data) != length(names_data_input))
    stop("data input is malformed. Its length is not the expected one. Use set_options() function to set it correctely.")
  if(!all(names(data) == names_data_input ))
    stop("data input parameter is malformed. The names are not the expected ones. Use set_options() function to set them correctely.")

  # get number of observations that will be clustered
  nobs = sum(data$n_j)

  #get number of individuals
  n = data$n

  #Check option to be in the correct form
  option_temp = set_options(partition = NULL)
  if(is.null(option)) # no option, set default
    option = option_temp

  if(length(option) != length(option_temp))
    stop("option parameter is malformed. Its length is not the expected one. Use set_options() function to set it correctely.")
  if(!all(names(option) == names(option_temp) ))
    stop("option parameter is malformed. The names are not the expected ones. Use set_options() function to set them correctely.")

  #Check partiton
  if(is.null(option$partition)){ # set empty partition
    if(FixPartition)
        stop("If FixPartition is selected, a partition must be provided in option$partition")
    option$partition = rep(0,n)
  }else{
    cat("\n Check that provided partition is well formed. It must start from 0 and all values must be contiguous \n")
    option$partition = arrange_partition(option$partition)

    # check that partiton and data are coherent
    if(nobs != length(option$partition))
      stop("The number of points in the data is not coherent with the length of the partition. Are there missing values in the data? Such implementation is not able to deal with them")
  }

  # Check initial values for tau
  K_init = length(table(option$partition)) # compute initial number of clusters
  if(is.null(option$init_var_cluster))
    option$init_var_cluster = 1/rgamma(n=K_init+option$Mstar0,
                                       shape = option$nu0/2,
                                       rate  = option$nu0*option$sigma0/2 )
  if(length(option$init_var_cluster)!=K_init+option$Mstar0)
    stop("The length of option$init_var_cluster must be equal to the initial number of clusters deduced from the initial partition plus Mstar0 ")
  if(is.null(option$init_mean_cluster))
    option$init_mean_cluster = rnorm(n=K_init+option$Mstar0,
                                     option$mu0, sqrt(option$init_var_cluster/option$k0))
  if(length(option$init_mean_cluster)!=K_init+option$Mstar0)
    stop("The length of option$init_mean_cluster must be equal to the initial number of clusters deduced from the initial partition plus Mstar0")

  # Check proposal for Mstar
  option$proposal_Mstar = floor(option$proposal_Mstar)  
  if(option$proposal_Mstar <= 0)
    stop("proposal_Mstar must be a strictly positive integer")

  # check covariates
  if(option$IncludeCovariates)
  {
    r = data$r
    if(length(option$beta0)!= r || nrow(option$Sigma0) != r || ncol(option$Sigma0) != r  )
      stop("The number of covariates r that is defined in data is not coherent with the size of beta0 or Sigma0 provided in options.")  
    if(  min(eigen(option$Sigma0)$values) <= 1e-14  )
      stop("Sigma0 is not positive definite or it is very ill conditioned.")
  }else{
    if(data$r != 0 )
      warning(" option$IncludeCovariates is set to FALSE but r > 0 ")
  }
  
  return( GDFMM:::MCMC_conditional_c(data, niter, burnin, thin, seed, P0.prior, FixPartition, option) )
}
