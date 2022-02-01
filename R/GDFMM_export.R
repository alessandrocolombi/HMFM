#' Title R function
#'
#' @param x input parameter
#' @return returns x plus 10
#' @export
hello2 <- function(data, niter, burnin, thin,seed, P0.prior = "Normal-InvGamma", option) {
  print("Hello, world!")
  return(x+10)
}


#' GDFMM Gibbs Sampler, all updates
#'
#' @param data input data.
#' @param niter number of iterations
#' @param burnin burnin period
#' @param thin thinning value
#' @param P0.prior string with the prior to be used as P0
#' @param option list with initial values, hyperparameters and other options
#' @return what it returns
#' @export
GDFMM_sampler <- function(data, niter, burnin, thin,seed, P0.prior = "Normal-InvGamma", option) {


  maxind <- function(vout) {
    maxind<-as.character(as.data.frame(table(vout))[which(as.data.frame(table(vout))$Freq==max(table(vout))),1])
    return(maxind)
  }
  cat('\n This is the R function: ')
  #Data check and pre-processing
  #--> handle here different types of input types. User you be able to pass the data the simplest possible ways. For example, this
  #    function should be able to handle both matrixes and data.frames (or others if needed).

  #Check number of iterations

  #Check P0.prior
  output<-GDFMM:::GDFMM_sampler_c(data, niter, burnin, thin,seed, P0.prior, option)
  #Check options
  K<-maxind(output$K)
  Mstar<-maxind(output$Mstar)
  Lambda<-mean(output$lambda) #shall we use median or mean

  # vec<-'vec'
  # for ( i in 1:length(output$C)){
  #   for ( j in j )
  #     for ( k in k )
  #   vec<-c(vec,output$C[[i]][j,k])
  #   maxind(vec)
  # }

  cat('Estimated K:',K,'\n')
  cat('Estimated Mstar:', Mstar, '\n')
  cat("Estimated Lambda:", Lambda, '\n')
  # This is just an example, of course you can save the c++ output and perform further operations in R
  return(output)
}




