#' Title R function
#'
#' @param x input parameter
#' @return returns x plus 10
#' @export
hello2 <- function(x) {
  print("Hello, world!")
  return(x+10)
}

#' The R function
#'
#' @param data input data.
#' @param niter number of iterations
#' @param burnin burnin period
#' @param thin thinning value
#' @param P0.prior string with the prior to be used as P0
#' @param option list with initial values, hyperparameters and other options
#' @return what it returns
#' @export
example_GDFMM_sampler <- function(data, niter, burnin, thin, P0.prior = "Normal-InvGamma", option) {

  cat('\n This is the R function: ')
  #Data check and pre-processing
  #--> handle here different types of input types. User you be able to pass the data the simplest possible ways. For example, this
  #    function should be able to handle both matrixes and data.frames (or others if needed).

  #Check number of iterations

  #Check P0.prior

  #Check options

  cat('Call the c++ function passing the preprocessed data, you can only pass types that can be traslated from R. \n')
  # This is just an example, of course you can save the c++ output and perform further operations in R
  return( GDFMM:::example_GDFMM_sampler_c(data, niter, burnin, thin, P0.prior))
}




