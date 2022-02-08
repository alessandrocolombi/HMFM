#' point_density : compute the point estimate with credible intervals for 
#'
#' @param x point where density is computed
#' @param mu_vec vector of mean for the gaussian distribution
#' @param sigma_vec vector of standard deviation for the gaussian distribution
#' @return c(inf, mean, sup) of point estimate of the density in x 
#' @export
point_density <- function(x, mu_vec, sigma_vec) {
  estimated_density = mean(dnorm(x, mu_vec, sigma_vec))
  err = 1/sd(dnorm(x, mu_vec, sigma_vec))*qnorm(0.975)
  ret = c(estimated_density - err, estimated_density, estimated_density + err)
  d_estimate = w_j[1]*dnorm(grid, mu[1], sigma[1])
  return(ret)
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



