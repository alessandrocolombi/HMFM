#' point_density : compute the point estimate with credible intervals for
#'
#' @param x point where density is computed
#' @param mu_vec vector of mean for the gaussian distribution
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
  #if(sum(w_j) != 1)
  # stop("weigths don't sum to one")

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

#' GDFMM Gibbs Sampler: options for fixing the Partition (FixPartition) and for changing prior (P0.prior)
#'
#' @param data input data
#' @param niter number of iterations
#' @param burnin burnin period
#' @param thin thinning value
#' @param seed seed for GSL random engine
#' @param P0.prior string with the prior to be used as P0
#' @param FixPartition TRUE if we want to fix the partition
#' @param option list with initial values, hyperparameters and other options
#' @return results of Gibbs Sampler
#' @export
GDFMM_sampler <- function(data, niter, burnin, thin, seed,
                            P0.prior = "Normal-InvGamma", FixPartition = F, option) {

  cat('\n This is the R function: ')
  #Data check and pre-processing
  #--> handle here different types of input types. User you be able to pass the data the simplest possible ways. For example, this
  #    function should be able to handle both matrixes and data.frames (or others if needed).

  #Check number of iterations

  #Check P0.prior

  #Check options
      # if(FixPartition){ check che una partition sia passata}
  cat('Call the c++ function passing the preprocessed data, you can only pass types that can be traslated from R. \n')
  # This is just an example, of course you can save the c++ output and perform further operations in R
  return( GDFMM:::GDFMM_sampler_c(data, niter, burnin, thin, seed, P0.prior, FixPartition, option))
}



