# Libraries

suppressWarnings(suppressPackageStartupMessages(library(GDFMM)))
suppressWarnings(suppressPackageStartupMessages(library(hdp)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(salso)))
suppressWarnings(suppressPackageStartupMessages(library(mcclust.ext)))

col_type = c("chartreuse3","orange","darkred","royalblue","cyan3")
# Functions ----------------------------------------------------------------

AM_density_estimation <- function(grid, fit, burnin = 1){
  n_iter <- length(fit$K) #number of iterations
  l_grid <- length(grid)
  #MIX    <- matrix(0, nrow=n_iter, ncol=l_grid)
  # MIX is a n_iter x l_grid matrix

  # This loop computes the predictive
  MIX = t(sapply(burnin:n_iter, simplify = "matrix",
                 function(it){
                   #cat("\n it = ",it,"\n")
                   # Get sampled values
                   M_it  = fit$M[it] # compute the number of components (allocated or not)
                   w_it = fit$W[[it]]
                   mu_it   <- unlist(fit$mu[[it]])          # get the mean, (mu_{1}^{(it)}, ..., mu_{M}^{(it)})
                   sig2_it <- unlist(fit$sig2[[it]])        # get the variances, (sigma^2_{1}^{(it)}, ..., sigma^2_{M}^{(it)})

                   # XX is a l_grid x M_it matrix, it contains the Normal kernels evauated over the grid
                   # XX[i,m] = Norm(grid[i] | mu_{m}^{(it)}, sigma^2_{m}^{(it)})
                   XX = t(sapply(1:M_it, simplify = "matrix",
                                 function(m){
                                   dnorm( x = grid, mean=mu_it[m], sd=sqrt(sig2_it[m]) ) # returns a vector of length equal to l_grid
                                 }
                   ))
                   # Compute predicted density at iteration it
                   (w_it) %*% XX
                 }
  ))


  # Density estimation and credible bounds
  pred_est <- apply(MIX,2,quantile,prob=c(0.025,0.5,0.975))
  return(pred_est)
}
SimStudy_Exp1 = function(seed){
  suppressWarnings(suppressPackageStartupMessages(library(GDFMM)))
  suppressWarnings(suppressPackageStartupMessages(library(hdp)))
  suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
  suppressWarnings(suppressPackageStartupMessages(library(salso)))
  suppressWarnings(suppressPackageStartupMessages(library(mcclust.ext)))

  # Initialize return object
  names_results = c("HMFM_marg","HDP","HMFM_cond","MFM","pooled")
  names_save = c("ARI_est_part","ARI_est_part_group1","ARI_est_part_group2",
                 "K_mode", "K1_mode", "K2_mode",
                 "K_ARI", "K1_ARI", "K2_ARI",
                 "err_coclust","err_coclust_group1","err_coclust_group2",
                 "err_L1_group1", "err_L1_group2")

  results = lapply(1:length(names_results),
                   FUN = function(x){temp = vector("list",length = length(names_save));names(temp) = names_save;temp})
  names(results) = names_results


  # Data Generation -------------------------------------------------------
  d = 2                   # number of groups
  K = 3                   # number of global clusters
  mu=c(-3,0,1.75)         # vectors of means
  sd = c(sqrt(0.1),sqrt(0.5),sqrt(1.5))       # vector of sd
  n_j=c(300,300)           # set cardinality of the groups
  mix_probs = matrix(c(0.5,0.5,0,
                       0,0.20,0.80), nrow = d, ncol = K, byrow = T)

  mix_obs = mix_probs * n_j

  genD = simulate_data(d = d, K = K, p = mix_probs, mu = mu, sd = sd, n_j = n_j, seed = seed)
  data = genD$data
  real_partition = genD$real_partition

  set.seed(seed)
  data[1, 1:mix_obs[1,1] ] = rnorm(n=mix_obs[1,1],mean = mu[1],sd = sd[1])
  real_partition[1:mix_obs[1,1]] = 1
  data[1, (mix_obs[1,1] + 1):( mix_obs[1,1] + mix_obs[1,2] )] = rnorm(n=mix_obs[1,2],mean = mu[2],sd = sd[2])
  real_partition[(mix_obs[1,1] + 1):( mix_obs[1,1] + mix_obs[1,2] ) ] = 2


  data[2, 1:mix_obs[2,2] ] = rnorm(n=mix_obs[2,2],mean = mu[2],sd = sd[2])
  real_partition[ (n_j[1] + 1): (n_j[1] + mix_obs[2,2] ) ] = 2
  data[2, (mix_obs[2,2] + 1 ):( n_j[2] )] = rnorm(n=mix_obs[2,3],mean = mu[3],sd = sd[3])
  real_partition[(n_j[1] + mix_obs[2,2] + 1 ):( n_j[1] + n_j[2] ) ] = 3


  data_list = vector("list",length = d)
  for(j in 1:d)
    data_list[[j]] = data[j,1:n_j[j]]

  # Set number of iterations
  niter  <-  30000
  burnin <-  10000
  thin   <-      1

  # a) HMFM - marginal ------------------------------------------------------
  ## Hyperparam

  ### $P_0$
  Range = range(unlist(data_list))
  mu0 = mean(unlist(data_list))
  R = Range[2] - Range[1]
  k0  = 1/R^2
  nu0 = 4
  sigma0 = 0.5
  scale = sqrt( (k0 + 1)/(k0) * sigma0 )
  mean_marginal = mu0
  var_marginal  = nu0/(nu0-2) * scale^2

  ### Process
  Exp_Lambda   = 5
  Var_Lambda   = 5
  gamma_guess  = 0.5
  Lambda_guess = Exp_Lambda
  b_lambda = Exp_Lambda/Var_Lambda
  a_lambda = Exp_Lambda * b_lambda
  a_gamma = a_lambda/d
  b_gamma = a_gamma / (gamma_guess * Lambda_guess)


  ## Run
  gamma0 = rep(0.001,d)
  Lambda0 = 5
  Ncenters = 20
  Kmeans0 = kmeans(x = unlist(data_list), centers = Ncenters, iter.max = 50, nstart = 10 )
  part0 = Kmeans0$cluster#centers = Kmeans0$centers
  gamma0 = rep(0.001,d)
  Lambda0 = 5
  Ncenters = 20
  Kmeans0 = kmeans(x = unlist(data_list), centers = Ncenters, iter.max = 50, nstart = 10 )
  part0 = Kmeans0$cluster#centers = Kmeans0$centers
  option = set_options_marginal(
    "mu0" = mu0,"sigma0"= sigma0, "k0"= k0, "nu0"=nu0,
    "Adapt_MH_hyp1"= 0.7,"Adapt_MH_hyp2"= 0.234,
    "sp_mala_U" = 0.01, "sp_mala_gamma"=0.01,
    "gamma0" = gamma0,"Lambda0" = Lambda0,
    "alpha_gamma"=a_gamma, "beta_gamma"=b_gamma,
    "alpha_lambda"=a_lambda, "beta_lambda"=b_lambda,
    "init_mean_cluster" = NULL,#unlist(mean_data_per_cluster),
    "init_var_cluster" = NULL,#unlist(var_data_per_cluster),
    "UpdateU" = T, "UpdateGamma" = T, "UpdateTau" = T, "UpdateLambda" = T,
    "partition" = part0
  )

  GDFMM = GDFMM_marginal_sampler(data, niter, burnin, thin, seed = seed, FixPartition = F, option = option)

  ## Clustering
  part_matrix <- GDFMM$Partition
  sim_matrix <- psm(part_matrix)

  VI_sara = minVI(sim_matrix)
  ARI_VI = arandi(VI_sara$cl,real_partition, adjust = TRUE)
  if(is.na(ARI_VI)){
    ARI_VI = arandi(VI_sara$cl,real_partition, adjust = FALSE)
    if(ARI_VI != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }
  ARI_VI_1  = arandi(VI_sara$cl[1:n_j[1]], real_partition[1:n_j[1]], adjust = TRUE)
  if(is.na(ARI_VI_1)){
    ARI_VI_1 = arandi(VI_sara$cl[1:n_j[1]], real_partition[1:n_j[1]], adjust = FALSE)
    if(ARI_VI_1 != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }
  ARI_VI_2  = arandi(VI_sara$cl[(n_j[1]+1):sum(n_j)],real_partition[(n_j[1]+1):sum(n_j)], adjust = TRUE)
  if(is.na(ARI_VI_2)){
    ARI_VI_2  = arandi(VI_sara$cl[(n_j[1]+1):sum(n_j)],real_partition[(n_j[1]+1):sum(n_j)], adjust = FALSE)
    if(ARI_VI_2 != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }

  K_ARI  = length(table(VI_sara$cl))
  K1_ARI = length(table(VI_sara$cl[1:n_j[1]]))
  K2_ARI = length(table(VI_sara$cl[(n_j[1]+1):sum(n_j)]))

  Local_Clustering = list("list", length = d)
  idx_start = c(1,cumsum(n_j)[1:(d-1)]+1)
  idx_end = cumsum(n_j)
  Kj = matrix(0,nrow = niter, ncol = d)
  for(j in 1:d){
    Local_Clustering[[j]] = GDFMM$Partition[ , idx_start[j]:idx_end[j] ]
    Kj[,j] = apply(Local_Clustering[[j]], 1, FUN = function(Part_it){length(table(Part_it))})
  }
  K_mod = as.numeric(which.max(table(GDFMM$K)))
  K1_mod = as.numeric(which.max(table(Kj[,1])))
  K2_mod = as.numeric(which.max(table(Kj[,2])))

  ## Density estimation
  xrange = c(-5,5)
  l_grid = 200
  grid = seq(xrange[1],xrange[2],length.out = l_grid)
  Pred_all = predictive_marginal_all_groups(grid = grid, fit = GDFMM, burnin = 0, option = option)

  ## Indici
  coclus_list  = Compute_coclust_error(real_partition, sim_matrix)
  coclus_list1 = Compute_coclust_error(real_partition[1:n_j[1]], sim_matrix[1:n_j[1],1:n_j[1]])
  coclus_list2 = Compute_coclust_error(real_partition[(1+n_j[1]):sum(n_j)], sim_matrix[(1+n_j[1]):sum(n_j),(1+n_j[1]):sum(n_j)])

  Pred_median = vector("list", length = d)
  for(j in 1:d){
    Pred_median[[j]] = Pred_all[[j]][2,]
  }

  L1_list = Compute_L1_dist(Pred = Pred_median,
                            p_mix = mix_probs, mu = mu, sigma = sd,
                            grid = grid)

  # save
  results$HMFM_marg$ARI_est_part = ARI_VI
  results$HMFM_marg$ARI_est_part_group1 = ARI_VI_1
  results$HMFM_marg$ARI_est_part_group2 = ARI_VI_2
  results$HMFM_marg$K_mode = K_mod
  results$HMFM_marg$K1_mode = K1_mod
  results$HMFM_marg$K2_mode = K2_mod
  results$HMFM_marg$K_ARI = K_ARI
  results$HMFM_marg$K1_ARI = K1_ARI
  results$HMFM_marg$K2_ARI = K2_ARI
  results$HMFM_marg$err_coclust = coclus_list$coclust_err
  results$HMFM_marg$err_coclust_group1 = coclus_list1$coclust_err
  results$HMFM_marg$err_coclust_group2 = coclus_list2$coclust_err
  results$HMFM_marg$err_L1_group1 = L1_list$L1err_per_level[1]
  results$HMFM_marg$err_L1_group2 = L1_list$L1err_per_level[2]

  # remove and finish
  rm(GDFMM, sim_matrix, Kj, Pred_all,Local_Clustering)


  # b) HDP ------------------------------------------------------------------
  ## Hyperparam

  Range = range(unlist(data_list))
  priorMean = mean(unlist(data_list))
  priorLambda = 1/(Range[2]-Range[1])^2
  priorA = 2
  priorB = 1

  sigma0 = priorB/priorA
  nu0 = 2*priorA


  scale = sqrt( (priorLambda + 1)/(priorLambda) * sigma0 )
  mean_marginal = priorMean
  var_marginal  = nu0/(nu0-2) * scale^2

  a_alpha = 1
  b_alpha = 1
  a_gamma = 1
  b_gamma = 0.1


  alpha0 = 1
  gamma0 = 1
  Niter   = niter
  Nburnin =  burnin

  fit = HDPMarginalSampler(Niter = Niter, Nburnin = Nburnin, d = d, n_j = n_j,
                           data_list = data_list,
                           priorMean = priorMean, priorA = priorA, priorB = priorB,
                           priorLambda = priorLambda,
                           a_gamma = a_gamma, b_gamma = b_gamma,
                           a_alpha = a_alpha, b_alpha = b_alpha,
                           alpha_init = alpha0, gamma_init = gamma0,
                           UpdateConc = TRUE, precompute_Stirling = TRUE)

  ## Clustering
  part_matrix <- fit$Partition
  sim_matrix <- psm(part_matrix)
  VI_sara = minVI(sim_matrix)
  ARI_VI = arandi(VI_sara$cl,real_partition, adjust = TRUE)
  if(is.na(ARI_VI)){
    ARI_VI = arandi(VI_sara$cl,real_partition, adjust = FALSE)
    if(ARI_VI != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }
  ARI_VI_1  = arandi(VI_sara$cl[1:n_j[1]], real_partition[1:n_j[1]], adjust = TRUE)
  if(is.na(ARI_VI_1)){
    ARI_VI_1 = arandi(VI_sara$cl[1:n_j[1]], real_partition[1:n_j[1]], adjust = FALSE)
    if(ARI_VI_1 != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }
  ARI_VI_2  = arandi(VI_sara$cl[(n_j[1]+1):sum(n_j)],real_partition[(n_j[1]+1):sum(n_j)], adjust = TRUE)
  if(is.na(ARI_VI_2)){
    ARI_VI_2  = arandi(VI_sara$cl[(n_j[1]+1):sum(n_j)],real_partition[(n_j[1]+1):sum(n_j)], adjust = FALSE)
    if(ARI_VI_2 != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }

  K_ARI = length(table(VI_sara$cl))
  K1_ARI = length(table(VI_sara$cl[1:n_j[1]]))
  K2_ARI = length(table(VI_sara$cl[(n_j[1]+1):sum(n_j)]))

  Local_Clustering = list("list", length = d)
  idx_start = c(1,cumsum(n_j)[1:(d-1)]+1)
  idx_end = cumsum(n_j)
  Kj = matrix(0,nrow = niter, ncol = d)
  for(j in 1:d){
    Local_Clustering[[j]] = fit$Partition[ , idx_start[j]:idx_end[j] ]
    Kj[,j] = apply(Local_Clustering[[j]], 1, FUN = function(Part_it){length(table(Part_it))})
  }
  K_mod = as.numeric(which.max(table(fit$K)))
  K1_mod = as.numeric(which.max(table(Kj[,1])))
  K2_mod = as.numeric(which.max(table(Kj[,2])))


  ## Density estimation
  xrange = c(-5,5)
  l_grid = 200
  grid = seq(xrange[1],xrange[2],length.out = l_grid)
  # Predictive in all groups
  Pred_all = hdp::predictive_all_groups(grid = grid, fit = fit,
                                        priorMean = priorMean, priorA = priorA, priorB = priorB,
                                        priorLambda = priorLambda, burnin = Niter/2)
  ## Indici
  coclus_list = Compute_coclust_error(real_partition, sim_matrix)
  coclus_list1 = Compute_coclust_error(real_partition[1:n_j[1]], sim_matrix[1:n_j[1],1:n_j[1]])
  coclus_list2 = Compute_coclust_error(real_partition[(1+n_j[1]):sum(n_j)], sim_matrix[(1+n_j[1]):sum(n_j),(1+n_j[1]):sum(n_j)])

  Pred_median = vector("list", length = d)
  for(j in 1:d){
    Pred_median[[j]] = Pred_all[[j]][2,]
  }

  L1_list = Compute_L1_dist(Pred = Pred_median,
                            p_mix = mix_probs, mu = mu, sigma = sd,
                            grid = grid)

  # save
  results$HDP$ARI_est_part = ARI_VI
  results$HDP$ARI_est_part_group1 = ARI_VI_1
  results$HDP$ARI_est_part_group2 = ARI_VI_2
  results$HDP$K_mode = K_mod
  results$HDP$K1_mode = K1_mod
  results$HDP$K2_mode = K2_mod
  results$HDP$K_ARI = K_ARI
  results$HDP$K1_ARI = K1_ARI
  results$HDP$K2_ARI = K2_ARI
  results$HDP$err_coclust = coclus_list$coclust_err
  results$HDP$err_coclust_group1 = coclus_list1$coclust_err
  results$HDP$err_coclust_group2 = coclus_list2$coclust_err
  results$HDP$err_L1_group1 = L1_list$L1err_per_level[1]
  results$HDP$err_L1_group2 = L1_list$L1err_per_level[2]

  # remove and finish
  rm(fit, sim_matrix, Kj, Pred_all, Local_Clustering)


  # c) HMFM - conditional ---------------------------------------------------
  names = c()
  SeasonNumber = c()
  for(j in 1:d){
    names = c(names, seq(1,n_j[j]))
    SeasonNumber = c(SeasonNumber, rep(j,n_j[j]))
  }

  data_small = tibble( ID = names,
                       SeasonNumber = SeasonNumber,
                       TrueClustering = real_partition,
                       Result = unlist(data_list))


  ## Hyperparameters

  ### $P_0$
  Res_range = range( data_small$Result )
  R = Res_range[2] - Res_range[1]
  mu0 = mean(data_small$Result) # should be 0
  k0  = 1/R^2
  nu0 = 4
  sigma0 = 0.5#*(100*R/1)

  scale = sqrt( (k0 + 1)/(k0) * sigma0 )
  mean_marginal = mu0
  var_marginal  = nu0/(nu0-2) * scale^2

  ### Process
  Exp_Lambda   =  5
  Var_Lambda   =  5
  gamma_guess  =  0.5
  Lambda_guess = Exp_Lambda

  b_lambda = Exp_Lambda/Var_Lambda
  a_lambda = Exp_Lambda * b_lambda

  a_gamma = a_lambda/d
  b_gamma = a_gamma / (gamma_guess * Lambda_guess)

  ## Initial values
  data_med4season = data_small %>% group_by(ID,SeasonNumber) %>%
    mutate(MedResult = Result) %>%
    select(ID,SeasonNumber,Result, MedResult) %>%
    distinct(ID,SeasonNumber, .keep_all = TRUE) %>%
    ungroup() %>% arrange(SeasonNumber)

  Ncenters = 20
  ymedian = data_med4season %>% pull(MedResult)
  Kmeans0 = kmeans(x = ymedian, centers = Ncenters, iter.max = 50, nstart = 10 )
  KmeansCl = Kmeans0$cluster
  centers = Kmeans0$centers
  data_med4season = data_med4season %>% cbind("Kmeans" = KmeansCl)

  data_with_init = data_small %>% left_join(data_med4season %>%
                                              select("ID","SeasonNumber","MedResult","Kmeans"),
                                            by = c("ID","SeasonNumber"))

  ## Setup
  dt = input_handle(data_with_init[,c(1,2,4,6)], intercept = FALSE)
  n = dt$n
  d = dt$d
  r = dt$r
  n_j = dt$n_j

  # initial values
  beta0 = rep(0,dt$r)
  Sigma0 = 100*diag(dt$r)

  Lambda0 = 10
  gamma0 = rep(0.25,d)
  Mstar0 = 0

  cluster_mean = data_med4season %>% group_by(Kmeans) %>%
    summarise(ClusterMean = mean(Result)) %>% ungroup() %>% pull(ClusterMean)
  cluster_var  = data_med4season %>% group_by(Kmeans) %>%
    summarise(ClusterVar = var(Result)) %>% ungroup() %>% pull(ClusterVar)
  initial_partition = unlist(unlist(dt$initialPartition))


  option = set_options( "mu0" = mu0,"sigma0"= sigma0, "k0"= k0, "nu0"=nu0,
                        "Adapt_MH_hyp2" = 0.234, "Adapt_MH_var0"=0.1,
                        "proposal_Mstar" = 1,
                        "Lambda0" = Lambda0, "gamma0" = gamma0, "Mstar0" = Mstar0,
                        "beta0" = beta0, "Sigma0" = Sigma0,
                        "alpha_gamma" = a_gamma, "beta_gamma" = b_gamma,
                        "alpha_lambda" = a_lambda, "beta_lambda" = b_lambda,
                        "init_mean_cluster" = c(centers, rep(0,Mstar0)),
                        "init_var_cluster" = c(cluster_var, rep(1,Mstar0)),
                        "partition" = initial_partition,
                        "IncludeCovariates" = FALSE,
                        "UpdateU" = T, "UpdateM" = T, "UpdateGamma" = T,
                        "UpdateS" = T, "UpdateTau" = T, "UpdateLambda" = T,
                        "UpdateBeta" = F )

  prior = "Normal-InvGamma"


  GDFMM = ConditionalSampler(dt[1:11], niter, burnin, thin, seed = 123, option = option, FixPartition = F,
                             P0.prior = prior)

  ## Clustering
  part_matrix <- GDFMM$Partition
  sim_matrix <- psm(part_matrix)
  VI_sara = minVI(sim_matrix)
  ARI_VI = arandi(VI_sara$cl,real_partition, adjust = TRUE)
  if(is.na(ARI_VI)){
    ARI_VI = arandi(VI_sara$cl,real_partition, adjust = FALSE)
    if(ARI_VI != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }
  ARI_VI_1  = arandi(VI_sara$cl[1:n_j[1]], real_partition[1:n_j[1]], adjust = TRUE)
  if(is.na(ARI_VI_1)){
    ARI_VI_1 = arandi(VI_sara$cl[1:n_j[1]], real_partition[1:n_j[1]], adjust = FALSE)
    if(ARI_VI_1 != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }
  ARI_VI_2  = arandi(VI_sara$cl[(n_j[1]+1):sum(n_j)],real_partition[(n_j[1]+1):sum(n_j)], adjust = TRUE)
  if(is.na(ARI_VI_2)){
    ARI_VI_2  = arandi(VI_sara$cl[(n_j[1]+1):sum(n_j)],real_partition[(n_j[1]+1):sum(n_j)], adjust = FALSE)
    if(ARI_VI_2 != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }

  K_ARI = length(table(VI_sara$cl))
  K1_ARI = length(table(VI_sara$cl[1:n_j[1]]))
  K2_ARI = length(table(VI_sara$cl[(n_j[1]+1):sum(n_j)]))

  ## Level dependent clustering
  Local_Clustering = list("list", length = d)
  idx_start = c(1,cumsum(n_j)[1:(d-1)]+1)
  idx_end = cumsum(n_j)
  Kj = matrix(0,nrow = niter, ncol = d)
  for(j in 1:d){
    Local_Clustering[[j]] = GDFMM$Partition[ , idx_start[j]:idx_end[j] ]
    Kj[,j] = apply(Local_Clustering[[j]], 1, FUN = function(Part_it){length(table(Part_it))})
  }
  K_mod = as.numeric(names(which.max(table(GDFMM$K))))
  K1_mod = as.numeric(names(which.max(table(Kj[,1]))))
  K2_mod = as.numeric(names(which.max(table(Kj[,2]))))

  ## Density estimation
  xrange = c(-5,5)
  l_grid = 200
  grid = seq(xrange[1],xrange[2],length.out = l_grid)
  # Predictive in all groups
  Pred_all = GDFMM::predictive_all_groups(grid = grid, fit = GDFMM, burnin = 1)
  ## Indici
  coclus_list = Compute_coclust_error(real_partition, sim_matrix)
  coclus_list1 = Compute_coclust_error(real_partition[1:n_j[1]], sim_matrix[1:n_j[1],1:n_j[1]])
  coclus_list2 = Compute_coclust_error(real_partition[(1+n_j[1]):sum(n_j)], sim_matrix[(1+n_j[1]):sum(n_j),(1+n_j[1]):sum(n_j)])

  Pred_median = vector("list", length = d)
  for(j in 1:d){
    Pred_median[[j]] = Pred_all[[j]][2,]
  }

  L1_list = Compute_L1_dist(Pred = Pred_median,
                            p_mix = mix_probs, mu = mu, sigma = sd,
                            grid = grid)

  # save
  results$HMFM_cond$ARI_est_part = ARI_VI
  results$HMFM_cond$ARI_est_part_group1 = ARI_VI_1
  results$HMFM_cond$ARI_est_part_group2 = ARI_VI_2
  results$HMFM_cond$K_mode = K_mod
  results$HMFM_cond$K1_mode = K1_mod
  results$HMFM_cond$K2_mode = K2_mod
  results$HMFM_cond$K_ARI = K_ARI
  results$HMFM_cond$K1_ARI = K1_ARI
  results$HMFM_cond$K2_ARI = K2_ARI
  results$HMFM_cond$err_coclust = coclus_list$coclust_err
  results$HMFM_cond$err_coclust_group1 = coclus_list1$coclust_err
  results$HMFM_cond$err_coclust_group2 = coclus_list2$coclust_err
  results$HMFM_cond$err_L1_group1 = L1_list$L1err_per_level[1]
  results$HMFM_cond$err_L1_group2 = L1_list$L1err_per_level[2]
  # remove and finish
  rm(GDFMM, sim_matrix, Kj, Pred_all,Local_Clustering)

  # d) MFM ------------------------------------------------------------------
  # First gruppo
  ## Hyperparam

  ### $P_0$
  Range = range(data[1,1:n_j[1]])
  mu0 = mean(data[1,1:n_j[1]])
  R = Range[2] - Range[1]
  k0  = 1/R^2
  nu0 = 4#10
  sigma0 = 0.5#10#*(100*R/1)
  scale = sqrt( (k0 + 1)/(k0) * sigma0 )
  mean_marginal = mu0
  var_marginal  = nu0/(nu0-2) * scale^2

  mixture_uvn_params = AntMAN::AM_mix_hyperparams_uninorm  (m0=mu0, k0=k0, nu0=nu0, sig02=sigma0)
  mcmc_params        = AntMAN::AM_mcmc_parameters(niter=niter, burnin=burnin, thin=10, verbose=1)
  components_prior   = AntMAN::AM_mix_components_prior_pois (init=3,  a=1, b=1)
  weights_prior      = AntMAN::AM_mix_weights_prior_gamma(init=2, a=1, b=1)

  fit <- AntMAN::AM_mcmc_fit(
    y = data[1,1:n_j[1]],
    mix_kernel_hyperparams = mixture_uvn_params,
    mix_components_prior =components_prior,
    mix_weight_prior = weights_prior,
    mcmc_parameters = mcmc_params)

  ## Clustering
  sim_matrix <- AntMAN::AM_coclustering(fit)

  VI_sara = minVI(sim_matrix)
  ARI_VI_1  = arandi(VI_sara$cl,real_partition[1:n_j[1]], adjust = TRUE)
  if(is.na(ARI_VI_1)){
    ARI_VI_1 = arandi(VI_sara$cl[1:n_j[1]], real_partition[1:n_j[1]], adjust = FALSE)
    if(ARI_VI_1 != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }

  K1_ARI = length(table(VI_sara$cl))
  K1_mod = as.numeric(names(which.max(table(fit$K))))

  ## Density estimation
  xrange = c(-5,5)
  l_grid = 200
  grid = seq(xrange[1],xrange[2],length.out = l_grid)
  Pred_all = AM_density_estimation(grid = grid, fit = fit, burnin = 1)

  ## Indici
  coclus_list1 = Compute_coclust_error(real_partition[1:n_j[1]], sim_matrix)
  Pred_median = vector("list", length = 1)
  Pred_median[[1]] = Pred_all[2,]

  L1_list1 = Compute_L1_dist(Pred = Pred_median,
                             p_mix = matrix(mix_probs[1,],nrow=1,ncol=K), mu = mu, sigma = sd,
                             grid = grid)

  # Second group
  ## Hyperparam
  ### $P_0$
  Range = range(data[2,1:n_j[2]])
  mu0 = mean(data[2,1:n_j[2]])
  R = Range[2] - Range[1]
  k0  = 1/R^2
  nu0 = 4#10
  sigma0 = 0.5#10#*(100*R/1)
  scale = sqrt( (k0 + 1)/(k0) * sigma0 )
  mean_marginal = mu0
  var_marginal  = nu0/(nu0-2) * scale^2

  mixture_uvn_params = AntMAN::AM_mix_hyperparams_uninorm  (m0=mu0, k0=k0, nu0=nu0, sig02=sigma0)
  mcmc_params        = AntMAN::AM_mcmc_parameters(niter=niter, burnin=burnin, thin=10, verbose=1)
  components_prior   = AntMAN::AM_mix_components_prior_pois (init=3,  a=1, b=1)
  weights_prior      = AntMAN::AM_mix_weights_prior_gamma(init=2, a=1, b=1)

  fit <- AntMAN::AM_mcmc_fit(
    y = data[2,1:n_j[2]],
    mix_kernel_hyperparams = mixture_uvn_params,
    mix_components_prior =components_prior,
    mix_weight_prior = weights_prior,
    mcmc_parameters = mcmc_params)

  ## Clustering
  sim_matrix <- AntMAN::AM_coclustering(fit)

  VI_sara = minVI(sim_matrix)
  ARI_VI_2  =  arandi(VI_sara$cl,real_partition[(n_j[1]+1):sum(n_j)], adjust = TRUE)
  if(is.na(ARI_VI_2)){
    ARI_VI_2  =  arandi(VI_sara$cl,real_partition[(n_j[1]+1):sum(n_j)], adjust = FALSE)
    if(ARI_VI_2 != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }

  K2_ARI = length(table(VI_sara$cl))
  K2_mod = as.numeric(names(which.max(table(fit$K))))

  ## Density estimation
  xrange = c(-5,5)
  l_grid = 200
  grid = seq(xrange[1],xrange[2],length.out = l_grid)
  Pred_all = AM_density_estimation(grid = grid, fit = fit, burnin = 1)

  ## Indici
  coclus_list2 = Compute_coclust_error(real_partition[(n_j[1]+1):sum(n_j)], sim_matrix)
  Pred_median = vector("list", length = 1)
  Pred_median[[1]] = Pred_all[2,]

  L1_list2 = Compute_L1_dist(Pred = Pred_median,
                             p_mix = matrix(mix_probs[2,],nrow=1,ncol=K), mu = mu, sigma = sd,
                             grid = grid)

  # save
  results$MFM$ARI_est_part = NaN
  results$MFM$ARI_est_part_group1 = ARI_VI_1
  results$MFM$ARI_est_part_group2 = ARI_VI_2
  results$MFM$K_mode = NaN
  results$MFM$K1_mode = K1_mod
  results$MFM$K2_mode = K2_mod
  results$MFM$K_ARI = NaN
  results$MFM$K1_ARI = K1_ARI
  results$MFM$K2_ARI = K2_ARI
  results$MFM$err_coclust = NaN
  results$MFM$err_coclust_group1 = coclus_list1$coclust_err
  results$MFM$err_coclust_group2 = coclus_list2$coclust_err
  results$MFM$err_L1_group1 = L1_list1$L1err_per_level
  results$MFM$err_L1_group2 = L1_list2$L1err_per_level

  # remove and finish
  rm(fit, sim_matrix, Pred_all)

  # e) MFM - pooled ---------------------------------------------------------
  ## Hyperparam
  ### $P_0$
  data_pooled = unlist(data_list)
  Range = range(data_pooled)
  mu0 = mean(data_pooled)
  R = Range[2] - Range[1]
  k0  = 1/R^2
  nu0 = 4#10
  sigma0 = 0.5#10#*(100*R/1)
  scale = sqrt( (k0 + 1)/(k0) * sigma0 )
  mean_marginal = mu0
  var_marginal  = nu0/(nu0-2) * scale^2

  mixture_uvn_params = AntMAN::AM_mix_hyperparams_uninorm  (m0=mu0, k0=k0, nu0=nu0, sig02=sigma0)
  mcmc_params        = AntMAN::AM_mcmc_parameters(niter=niter, burnin=burnin, thin=10, verbose=1)
  components_prior   = AntMAN::AM_mix_components_prior_pois (init=3,  a=1, b=1)
  weights_prior      = AntMAN::AM_mix_weights_prior_gamma(init=2, a=1, b=1)

  fit <- AntMAN::AM_mcmc_fit(
    y = data_pooled,
    mix_kernel_hyperparams = mixture_uvn_params,
    mix_components_prior =components_prior,
    mix_weight_prior = weights_prior,
    mcmc_parameters = mcmc_params)

  ## Clustering
  sim_matrix <- AntMAN::AM_coclustering(fit)

  VI_sara = minVI(sim_matrix)
  ARI_VI  = arandi(VI_sara$cl,real_partition, adjust = TRUE)
  if(is.na(ARI_VI)){
    ARI_VI  = arandi(VI_sara$cl,real_partition, adjust = TRUE)
    if(ARI_VI != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }
  ARI_VI_1  = arandi(VI_sara$cl[1:n_j[1]], real_partition[1:n_j[1]], adjust = TRUE)
  if(is.na(ARI_VI_1)){
    ARI_VI_1 = arandi(VI_sara$cl[1:n_j[1]], real_partition[1:n_j[1]], adjust = FALSE)
    if(ARI_VI_1 != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }
  ARI_VI_2  = arandi(VI_sara$cl[(n_j[1]+1):sum(n_j)],real_partition[(n_j[1]+1):sum(n_j)], adjust = TRUE)
  if(is.na(ARI_VI_2)){
    ARI_VI_2  = arandi(VI_sara$cl[(n_j[1]+1):sum(n_j)],real_partition[(n_j[1]+1):sum(n_j)], adjust = FALSE)
    if(ARI_VI_2 != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }

  K_ARI = length(table(VI_sara$cl))
  K_mod = as.numeric(names(which.max(table(fit$K))))

  ## Indici
  coclus_list = Compute_coclust_error(real_partition, sim_matrix)
  coclus_list1 = Compute_coclust_error(real_partition[1:n_j[1]], sim_matrix[1:n_j[1],1:n_j[1]])
  coclus_list2 = Compute_coclust_error(real_partition[(1+n_j[1]):sum(n_j)], sim_matrix[(1+n_j[1]):sum(n_j),(1+n_j[1]):sum(n_j)])

  # save
  results$pooled$ARI_est_part = ARI_VI
  results$pooled$ARI_est_part_group1 = ARI_VI_1
  results$pooled$ARI_est_part_group2 = ARI_VI_2
  results$pooled$K_mode = K_mod
  results$pooled$K1_mode = NaN
  results$pooled$K2_mode = NaN
  results$pooled$K_ARI = K_ARI
  results$pooled$K1_ARI = NaN
  results$pooled$K2_ARI = NaN
  results$pooled$err_coclust = coclus_list$coclust_err
  results$pooled$err_coclust_group1 = coclus_list1$coclust_err
  results$pooled$err_coclust_group2 = coclus_list2$coclust_err
  results$pooled$err_L1_group1 = NaN
  results$pooled$err_L1_group2 = NaN


  # remove and finish
  rm(fit, sim_matrix)


  # return ------------------------------------------------------------------
  return(results)

}


# Example plot -----------------------------------------------------------------
seed = 1234
d = 2                   # number of groups
K = 3                   # number of global clusters
mu=c(-3,0,1.75)         # vectors of means
sd = c(sqrt(0.1),sqrt(0.5),sqrt(1.5))       # vector of sd
n_j=c(300,300)           # set cardinality of the groups
mix_probs = matrix(c(0.5,0.5,0,
                     0,0.20,0.80), nrow = d, ncol = K, byrow = T)

mix_obs = mix_probs * n_j

genD = simulate_data(d = d, K = K, p = mix_probs, mu = mu, sd = sd, n_j = n_j, seed = seed)
data = genD$data
real_partition = genD$real_partition

set.seed(seed)
data[1, 1:mix_obs[1,1] ] = rnorm(n=mix_obs[1,1],mean = mu[1],sd = sd[1])
real_partition[1:mix_obs[1,1]] = 1
data[1, (mix_obs[1,1] + 1):( mix_obs[1,1] + mix_obs[1,2] )] = rnorm(n=mix_obs[1,2],mean = mu[2],sd = sd[2])
real_partition[(mix_obs[1,1] + 1):( mix_obs[1,1] + mix_obs[1,2] ) ] = 2


data[2, 1:mix_obs[2,2] ] = rnorm(n=mix_obs[2,2],mean = mu[2],sd = sd[2])
real_partition[ (n_j[1] + 1): (n_j[1] + mix_obs[2,2] ) ] = 2
data[2, (mix_obs[2,2] + 1 ):( n_j[2] )] = rnorm(n=mix_obs[2,3],mean = mu[3],sd = sd[3])
real_partition[(n_j[1] + mix_obs[2,2] + 1 ):( n_j[1] + n_j[2] ) ] = 3


data_list = vector("list",length = d)
for(j in 1:d)
  data_list[[j]] = data[j,1:n_j[j]]


names = c()
GroupNumber = c()
for(j in 1:d){
  names = c(names, seq(1,n_j[j]))
  GroupNumber = c(GroupNumber, rep(j,n_j[j]))
}

data_small = tibble( ID = names,
                     Group = GroupNumber,
                     TrueClustering = real_partition,
                     Value = unlist(data_list))




mycol = hcl.colors(n=3,palette = "Zissou1")
mycol_cluster = hcl.colors(n=K, palette = "Temps")

idx_start = 1
idx_end = n_j[1]
xrange = c(-5,5)
l_grid = 200
grid = seq(xrange[1],xrange[2],length.out = l_grid)


ylim_list = vector("list",2)
ylim_list[[1]] = c(0,1.5)
ylim_list[[2]] = c(0,0.75)

par(mfrow = c(1,2), mar = c(2,2,2,1), bty = "l")
for(j in 1:d){
  plot(0,0,main = paste0("Group ",j),xlab = " ", type = "n", xlim = xrange, ylim = ylim_list[[j]])
  grid(lty = 1,lwd = 1, col = "gray90" )
  for(k in 1:K){
    res = data_small %>% filter(Group == j) %>% filter(TrueClustering == k) %>% pull(Value)
    if(length(res)>0){
      hist(res, freq = FALSE, nclass = "fd",
           col= ACutils::t_col(mycol_cluster[k], percent = 60), add = T)
      points( x = res, y = rep(0,length(res) ), pch = 16, col = mycol_cluster[k] )
    }
  }
  points(grid, GDFMM::dmix(x = grid, w_j = mix_probs[j,], mu_vec = mu, sigma_vec = sd),
         col = "red", lwd = 2, type = "l")

}

# Pooled data
par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
plot(0,0,main = paste0("Pooled data"),xlab = " ", type = "n", xlim = xrange, ylim = c(0,1.5))
grid(lty = 1,lwd = 1, col = "gray90" )
for(k in 1:K){
  res = data_small %>% filter(TrueClustering == k) %>% pull(Value)
  if(length(res)>0){
    hist(res, freq = FALSE, nclass = "fd",
         col= ACutils::t_col(mycol_cluster[k], percent = 60), add = T)
    points( x = res, y = rep(0,length(res) ), pch = 16, col = mycol_cluster[k] )
  }
}
dens_pooled = density(data_small %>% pull(Value))
points(dens_pooled$x, dens_pooled$y, col = "red", lwd = 2, type = "l")


# Run ----------------------------------------------------------

Nrep  = 50

seed0 = 1605
set.seed(seed0)
seeds = sample(1:999999, size = Nrep)
num_cores = 7

tictoc::tic()
  cluster <- parallel::makeCluster(num_cores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)
  parallel::clusterExport(cluster, list("AM_density_estimation"))
  res = parallel::parLapply( cl = cluster, seeds,
                             fun = SimStudy_Exp1)
  parallel::stopCluster(cluster)
tictoc::toc()


## PS1
name = "err_L1_group1"
ylabel = "PS group 1"

HDP_res = sapply(res, function(x){x$HDP$err_L1_group1})
HMFMmarg_res = sapply(res, function(x){x$HMFM_marg$err_L1_group1})
HMFMcond_res = sapply(res, function(x){x$HMFM_cond$err_L1_group1})
MFM_res = sapply(res, function(x){x$MFM$err_L1_group1})

exp_temp = tibble("err_L1_group1" = HDP_res, "type" = as_factor("HDP"))
exp_temp = exp_temp %>%
  rbind(tibble("err_L1_group1" = HMFMmarg_res, "type" = as_factor("HMFM-marg"))) %>%
  rbind(tibble("err_L1_group1" = HMFMcond_res, "type" = as_factor("HMFM-cond"))) %>%
  rbind(tibble("err_L1_group1" = MFM_res, "type" = as_factor("MFM")))


PS1_plot1 = exp_temp %>% select(type,!!name) %>%
  ggplot(aes(y=!!sym(name), x=type, fill=type)) + geom_boxplot(fill = col_type[1:4]) +
  labs(y=paste0(ylabel,", n = ",sum(n_j)), x = " ") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        text = element_text(size = 10))  + ylim(c(0,0.2))


## PS2
name = "err_L1_group2"
ylabel = "PS group 2"

HDP_res = sapply(res, function(x){x$HDP$err_L1_group2})
HMFMmarg_res = sapply(res, function(x){x$HMFM_marg$err_L1_group2})
HMFMcond_res = sapply(res, function(x){x$HMFM_cond$err_L1_group2})
MFM_res = sapply(res, function(x){x$MFM$err_L1_group2})

exp_temp = tibble("err_L1_group2" = HDP_res, "type" = as_factor("HDP"))
exp_temp = exp_temp %>%
  rbind(tibble("err_L1_group2" = HMFMmarg_res, "type" = as_factor("HMFM-marg"))) %>%
  rbind(tibble("err_L1_group2" = HMFMcond_res, "type" = as_factor("HMFM-cond"))) %>%
  rbind(tibble("err_L1_group2" = MFM_res, "type" = as_factor("MFM")))


PS2_plot1 = exp_temp %>% select(type,!!name) %>%
  ggplot(aes(y=!!sym(name), x=type, fill=type)) + geom_boxplot(fill = col_type[1:4]) +
  labs(y=paste0(ylabel,", n = ",sum(n_j)), x = " ") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        text = element_text(size = 10))  + ylim(c(0,0.2))

ggpubr::ggarrange(PS1_plot1, PS2_plot1, ncol = 2)

## CCE group 1
name = "err_coclust_group1"
ylabel = "CCE group 1"

HDP_res = sapply(res, function(x){x$HDP$err_coclust_group1})
HMFMmarg_res = sapply(res, function(x){x$HMFM_marg$err_coclust_group1})
HMFMcond_res = sapply(res, function(x){x$HMFM_cond$err_coclust_group1})
MFM_res = sapply(res, function(x){x$MFM$err_coclust_group1})
pooled_res = sapply(res, function(x){x$pooled$err_coclust_group1})


exp_temp = tibble("err_coclust_group1" = HDP_res, "type" = as_factor("HDP"))
exp_temp = exp_temp %>%
  rbind(tibble("err_coclust_group1" = HMFMmarg_res, "type" = as_factor("HMFM-marg"))) %>%
  rbind(tibble("err_coclust_group1" = HMFMcond_res, "type" = as_factor("HMFM-cond"))) %>%
  rbind(tibble("err_coclust_group1" = MFM_res, "type" = as_factor("MFM"))) %>%
  rbind(tibble("err_coclust_group1" = pooled_res,   "type" = as_factor("pooled")))

CCE1_plot1 = exp_temp %>% select(type,!!name) %>%
  ggplot(aes(y=!!sym(name), x=type, fill=type)) + geom_boxplot(fill = col_type) +
  labs(y=paste0(ylabel,", n = ",sum(n_j)), x = " ") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        text = element_text(size = 10))  + ylim(c(0,0.5))



## CCE group 2
name = "err_coclust_group2"
ylabel = "CCE group 2"

HDP_res = sapply(res, function(x){x$HDP$err_coclust_group2})
HMFMmarg_res = sapply(res, function(x){x$HMFM_marg$err_coclust_group2})
HMFMcond_res = sapply(res, function(x){x$HMFM_cond$err_coclust_group2})
MFM_res = sapply(res, function(x){x$MFM$err_coclust_group2})
pooled_res = sapply(res, function(x){x$pooled$err_coclust_group2})


exp_temp = tibble("err_coclust_group2" = HDP_res, "type" = as_factor("HDP"))
exp_temp = exp_temp %>%
  rbind(tibble("err_coclust_group2" = HMFMmarg_res, "type" = as_factor("HMFM-marg"))) %>%
  rbind(tibble("err_coclust_group2" = HMFMcond_res, "type" = as_factor("HMFM-cond"))) %>%
  rbind(tibble("err_coclust_group2" = MFM_res, "type" = as_factor("MFM"))) %>%
  rbind(tibble("err_coclust_group2" = pooled_res,   "type" = as_factor("pooled")))

CCE2_plot1 = exp_temp %>% select(type,!!name) %>%
  ggplot(aes(y=!!sym(name), x=type, fill=type)) + geom_boxplot(fill = col_type) +
  labs(y=paste0(ylabel,", n = ",sum(n_j)), x = " ") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        text = element_text(size = 10))  + ylim(c(0,0.5))


ggpubr::ggarrange(CCE1_plot1, CCE2_plot1, ncol = 2)

## K_ARI
name = "K_ARI"
ylabel = "Est. K"

HDP_res = sapply(res, function(x){x$HDP$K_ARI})
HMFMmarg_res = sapply(res, function(x){x$HMFM_marg$K_ARI})
HMFMcond_res = sapply(res, function(x){x$HMFM_cond$K_ARI})
pooled_res = sapply(res, function(x){x$pooled$K_ARI})

exp_temp = tibble("K_ARI" = HDP_res, "type" = as_factor("HDP"))
exp_temp = exp_temp %>%
  rbind(tibble("K_ARI" = HMFMmarg_res, "type" = as_factor("HMFM-marg"))) %>%
  rbind(tibble("K_ARI" = HMFMcond_res, "type" = as_factor("HMFM-cond"))) %>%
  rbind(tibble("K_ARI" = pooled_res,   "type" = as_factor("pooled")))


K_ARI_plot1 = ggplot(exp_temp, aes(x = K_ARI, fill = type)) +
  geom_bar(position = "dodge") + theme_bw() +
  scale_fill_manual(values = c("HDP" = col_type[1],
                               "HMFM-marg" = col_type[2],
                               "HMFM-cond" = col_type[3],
                               "pooled" = col_type[5])) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        text = element_text(size = 10)) +
  labs(y=paste0(ylabel,", n = ",sum(n_j)), x = " ")
# K_ARI_plot1


## K_mode
name = "K_mode"
ylabel = "Est. K"

HDP_res = sapply(res, function(x){x$HDP$K_mode})
HMFMmarg_res = sapply(res, function(x){x$HMFM_marg$K_mode})
HMFMcond_res = sapply(res, function(x){x$HMFM_cond$K_mode})
pooled_res   = sapply(res, function(x){x$pooled$K_mode})


exp_temp = tibble("K_mode" = HDP_res, "type" = as_factor("HDP"))
exp_temp = exp_temp %>%
  rbind(tibble("K_mode" = HMFMmarg_res, "type" = as_factor("HMFM-marg"))) %>%
  rbind(tibble("K_mode" = HMFMcond_res, "type" = as_factor("HMFM-cond"))) %>%
  rbind(tibble("K_mode" = HMFMcond_res, "type" = as_factor("pooled")))


K_mode_plot1 = ggplot(exp_temp, aes(x = K_mode, fill = type)) +
  geom_bar(position = "dodge") + theme_bw() +
  scale_fill_manual(values = c("HDP" = col_type[1],
                               "HMFM-marg" = col_type[2],
                               "HMFM-cond" = col_type[3],
                               "pooled" = col_type[5])) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        text = element_text(size = 10)) +
  labs(y=paste0(ylabel,", n = ",sum(n_j)), x = " ")
# K_mode_plot1


# Table group 1
HDP_res = sapply(res, function(x){x$HDP$ARI_est_part_group1})
HMFMmarg_res = sapply(res, function(x){x$HMFM_marg$ARI_est_part_group1})
HMFMcond_res = sapply(res, function(x){x$HMFM_cond$ARI_est_part_group1})
MFM_res = sapply(res, function(x){x$MFM$ARI_est_part_group1})
pooled_res = sapply(res, function(x){x$pooled$ARI_est_part_group1})

cat("\n GROUP 1 \n")
cat("\n","HDP ARI: mean = ",mean(HDP_res),"; sd = ",sd(HDP_res),"\n")
cat("\n","HMFM-marg ARI: mean = ",mean(HMFMmarg_res),"; sd = ",sd(HMFMmarg_res),"\n")
cat("\n","HMFM-cond ARI: mean = ",mean(HMFMcond_res),"; sd = ",sd(HMFMcond_res),"\n")
cat("\n","MFM ARI: mean = ",mean(MFM_res),"; sd = ",sd(MFM_res),"\n")
cat("\n","pooled ARI: mean = ",mean(pooled_res),"; sd = ",sd(pooled_res),"\n")

# Table group 2
HDP_res = sapply(res, function(x){x$HDP$ARI_est_part_group2})
HMFMmarg_res = sapply(res, function(x){x$HMFM_marg$ARI_est_part_group2})
HMFMcond_res = sapply(res, function(x){x$HMFM_cond$ARI_est_part_group2})
MFM_res = sapply(res, function(x){x$MFM$ARI_est_part_group2})
pooled_res = sapply(res, function(x){x$pooled$ARI_est_part_group2})

cat("\n GROUP 2 \n")
cat("\n","HDP ARI: mean = ",mean(HDP_res),"; sd = ",sd(HDP_res),"\n")
cat("\n","HMFM-marg ARI: mean = ",mean(HMFMmarg_res),"; sd = ",sd(HMFMmarg_res),"\n")
cat("\n","HMFM-cond ARI: mean = ",mean(HMFMcond_res),"; sd = ",sd(HMFMcond_res),"\n")
cat("\n","MFM ARI: mean = ",mean(MFM_res),"; sd = ",sd(MFM_res),"\n")
cat("\n","pooled ARI: mean = ",mean(pooled_res),"; sd = ",sd(pooled_res),"\n")


# 1 Cluster %
HDP_res = sum(sapply(res, function(x){x$HDP$K2_ARI == 1}))/Nrep
HMFMmarg_res = sum(sapply(res, function(x){x$HMFM_marg$K2_ARI == 1}))/Nrep
HMFMcond_res = sum(sapply(res, function(x){x$HMFM_cond$K2_ARI == 1}))/Nrep
MFM_res = sum(sapply(res, function(x){x$MFM$K2_ARI == 1}))/Nrep

cat("\n % 1 cluster \n")
cat("\n HDP : ",HDP_res * 100, "% \n")
cat("\n HMFM-marg : ",HMFMmarg_res * 100, "% \n")
cat("\n HMFM-cond : ",HMFMcond_res * 100, "% \n")
cat("\n MFM : ",MFM_res * 100, "% \n")
cat("\n pooled : NOT DEFINED \n")


beepr::beep()
