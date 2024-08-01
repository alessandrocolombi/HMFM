# Libraries

suppressWarnings(suppressPackageStartupMessages(library(GDFMM)))
suppressWarnings(suppressPackageStartupMessages(library(hdp)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(salso)))
suppressWarnings(suppressPackageStartupMessages(library(mcclust.ext)))

# col_type = c("chartreuse3","orange","darkred","cyan3")
col_type = c("darkred","orange","chartreuse3","cyan3")

# Functions ----------------------------------------------------------------


SimStudy_Exp3 = function(seed){
  suppressWarnings(suppressPackageStartupMessages(library(GDFMM)))
  suppressWarnings(suppressPackageStartupMessages(library(hdp)))
  suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
  suppressWarnings(suppressPackageStartupMessages(library(salso)))
  suppressWarnings(suppressPackageStartupMessages(library(mcclust.ext)))

  # Initialize return object
  names_results = c("HMFM_marg","HDP","HMFM_cond","pooled")
  names_save = c("ARI_est_part","ARI_est_part_local",
                 "K_mode", "K_mode_local",
                 "K_ARI", "K_ARI_local",
                 "err_coclust","err_coclust_local",
                 "err_L1_local", "err_L1_mean")

  results = lapply(1:length(names_results),
                   FUN = function(x){temp = vector("list",length = length(names_save));names(temp) = names_save;temp})
  names(results) = names_results




  # Data Generation -------------------------------------------------------
  #1) Data Generation
  d = 15
  K = 5
  n_j = rep(30,d)
  d1 = 12
  d2 = d-d1
  set.seed(seed)
  # first mixture
  mu1 = c(-3,0,3)         # vectors of means
  sd1 = c(sqrt(0.5),sqrt(0.5),sqrt(0.5))       # vector of sd
  mix_probs1_vec = c(0.25,0.5,0.25)
  mix_probs1 = matrix(0,nrow = d1, ncol = length(mu1))
  K_j = rep(0,d)
  set.seed(123545) # This ensures that K_j does not change! Either components_j will change
  for(j in 1:d1){
    K_j[j] = sample(2:3, size = 1)
    components_j = sort(sample(1:3, size = K_j[j]))
    mix_probs_j = rep(0,3)
    mix_probs_j[components_j] = mix_probs1_vec[components_j]
    mix_probs_j = mix_probs_j/sum(mix_probs_j)
    mix_probs1[j,] = mix_probs_j
  }
  genD1 = simulate_data(d = d1, K = 3, p = mix_probs1,
                        mu = mu1, sd = sd1, n_j = n_j[1:d1],
                        seed = seed)
  data1 = genD1$data
  real_partition1 = genD1$real_partition


  # second mixture
  mu2 = c(-1.5,1.5)         # vectors of means
  sd2 = c(sqrt(0.5),sqrt(0.5))       # vector of sd
  mix_probs2_vec = c(0.5,0.5)
  mix_probs2 = matrix(0,nrow = d2, ncol = length(mu2))
  K_j[(d1+1):d] = 2
  for(j in 1:d2){
    mix_probs2[j,] = mix_probs2_vec
  }

  genD2 = simulate_data(d = d2, K = 2, p = mix_probs2,
                        mu = mu2, sd = sd2, n_j = n_j[(d1+1):d],
                        seed = seed)
  data2 = genD2$data
  real_partition2 = genD2$real_partition
  real_partition2 = real_partition2 + length(mu1)

  data = rbind(data1,data2)
  real_partition = c(real_partition1,real_partition2)

  data_list = vector("list",length = d)
  for(j in 1:d)
    data_list[[j]] = data[j,1:n_j[j]]


  # Set number of iterations
  niter  <-  35000
  burnin <-  25000
  thin   <-      1

  set.seed(seed) # just in case
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
  Exp_Lambda   = 15
  Var_Lambda   =  3
  gamma_guess  =  0.05
  Lambda_guess = Exp_Lambda
  b_lambda = Exp_Lambda/Var_Lambda
  a_lambda = Exp_Lambda * b_lambda
  a_gamma = a_lambda/d
  b_gamma = a_gamma / (gamma_guess * Lambda_guess)


  ## Run
  gamma0 = rep(0.05,d)
  Lambda0 = 15
  Ncenters = 20
  Kmeans0 = kmeans(x = unlist(data_list), centers = Ncenters, iter.max = 50, nstart = 10 )
  part0 = Kmeans0$cluster#centers = Kmeans0$centers
  option = set_options_marginal(
    "mu0" = mu0,"sigma0"= sigma0, "k0"= k0, "nu0"=nu0,
    "Adapt_MH_hyp1"= 0.7,"Adapt_MH_hyp2"= 0.234,
    "sp_mala_U" = 0.05, "sp_mala_gamma"=0.01,
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
  ARI_VI  = arandi(VI_sara$cl,real_partition, adjust = TRUE)
  if(is.na(ARI_VI)){
    ARI_VI  = arandi(VI_sara$cl,real_partition, adjust = FALSE)
    if(ARI_VI != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }

  Local_Clustering = list("list", length = d)
  idx_start = c(1,cumsum(n_j)[1:(d-1)]+1)
  idx_end = cumsum(n_j)
  Kj = matrix(0,nrow = niter, ncol = d)
  ARI_j = vector(length = d)
  Kj_ARI = vector(length = d)
  for(j in 1:d){
    Local_Clustering[[j]] = GDFMM$Partition[ , idx_start[j]:idx_end[j] ]
    Kj[,j] = apply(Local_Clustering[[j]], 1, FUN = function(Part_it){length(table(Part_it))})
    ARI_j[j]  = arandi(  VI_sara$cl[idx_start[j]:idx_end[j]],
                         real_partition[idx_start[j]:idx_end[j]],
                         adjust = TRUE  )
    if(is.na(ARI_j[j])){
      ARI_j[j]  = arandi(  VI_sara$cl[idx_start[j]:idx_end[j]],
                           real_partition[idx_start[j]:idx_end[j]],
                           adjust = FALSE  )
      if(ARI_j[j] != 1)
        stop("Adjusted RI is NaN but Unadjusted RI is not 1")
    }

    Kj_ARI[j] = length(table(VI_sara$cl[idx_start[j]:idx_end[j]]))

  }
  K_mod  = as.numeric(names(which.max(table(GDFMM$K))))
  Kj_mod = unlist(lapply(1:d, function(jj){as.numeric(names(which.max(table(Kj[,jj]))))}))
  K_ARI  = length(table(VI_sara$cl))

  ## Density estimation
  xrange = c(-7,7)
  l_grid = 200
  grid = seq(xrange[1],xrange[2],length.out = l_grid)
  Pred_all = predictive_marginal_all_groups(grid = grid, fit = GDFMM, burnin = 0, option = option)

  ## Indici
  coclus_list = Compute_coclust_error(real_partition, sim_matrix)
  local_coclus_list = vector("list", length = d)
  Pred_median = vector("list", length = d)
  for(j in 1:d){
    Pred_median[[j]] = Pred_all[[j]][2,]
    take_indices = idx_start[j]:idx_end[j]
    local_sim_matrix = psm(GDFMM$Partition[,take_indices])
    local_coclus_list[[j]] = Compute_coclust_error( real_partition[take_indices],local_sim_matrix )
  }
  local_coclus_j = sapply(local_coclus_list, FUN = function(xx){xx[[1]]})

  # Compute L1 distance wrt to true density, which is the same for all first d1 groups
  L1_list1 = Compute_L1_dist(Pred = Pred_median[1:d1],
                             p_mix = mix_probs1, mu = mu1, sigma = sd1,
                             grid = grid)
  # Compute L1 distance wrt to true density, which is the same for all final d-d1 groups
  L1_list2 = Compute_L1_dist(Pred = Pred_median[(1+d1):d],
                             p_mix = mix_probs2, mu = mu2, sigma = sd2,
                             grid = grid)

  # save
  results$HMFM_marg$ARI_est_part = ARI_VI
  results$HMFM_marg$ARI_est_part_local = ARI_j
  results$HMFM_marg$K_mode = K_mod
  results$HMFM_marg$K_mode_local = Kj_mod
  results$HMFM_marg$K_ARI = K_ARI
  results$HMFM_marg$K_ARI_local = Kj_ARI
  results$HMFM_marg$err_coclust = coclus_list$coclust_err
  results$HMFM_marg$err_coclust_local = local_coclus_j
  results$HMFM_marg$err_L1_local = c(L1_list1[[1]],L1_list2[[1]])
  results$HMFM_marg$err_L1_mean = (d1 * L1_list1[[2]] + d2 * L1_list2[[2]])/d

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
    ARI_VI  = arandi(VI_sara$cl,real_partition, adjust = FALSE)
    if(ARI_VI != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }

  Local_Clustering = list("list", length = d)
  idx_start = c(1,cumsum(n_j)[1:(d-1)]+1)
  idx_end = cumsum(n_j)
  Kj = matrix(0,nrow = niter, ncol = d)
  ARI_j = vector(length = d)
  Kj_ARI = vector(length = d)
  for(j in 1:d){
    Local_Clustering[[j]] = fit$Partition[ , idx_start[j]:idx_end[j] ]
    Kj[,j] = apply(Local_Clustering[[j]], 1, FUN = function(Part_it){length(table(Part_it))})
    ARI_j[j]  = arandi(  VI_sara$cl[idx_start[j]:idx_end[j]],
                         real_partition[idx_start[j]:idx_end[j]],
                         adjust = TRUE  )
    if(is.na(ARI_j[j])){
      ARI_j[j]  = arandi(  VI_sara$cl[idx_start[j]:idx_end[j]],
                           real_partition[idx_start[j]:idx_end[j]],
                           adjust = FALSE  )
      if(ARI_j[j] != 1)
        stop("Adjusted RI is NaN but Unadjusted RI is not 1")
    }
    Kj_ARI[j] = length(table(VI_sara$cl[idx_start[j]:idx_end[j]]))

  }
  K_mod  = as.numeric(names(which.max(table(fit$K))))
  Kj_mod = unlist(lapply(1:d, function(jj){as.numeric(names(which.max(table(Kj[,jj]))))}))
  K_ARI  = length(table(VI_sara$cl))
  ## Density estimation
  xrange = c(-7,7)
  l_grid = 200
  grid = seq(xrange[1],xrange[2],length.out = l_grid)
  # Predictive in all groups
  Pred_all = hdp::predictive_all_groups(grid = grid, fit = fit,
                                        priorMean = priorMean, priorA = priorA, priorB = priorB,
                                        priorLambda = priorLambda, burnin = 1)


  ## Indici
  coclus_list = Compute_coclust_error(real_partition, sim_matrix)
  local_coclus_list = vector("list", length = d)
  Pred_median = vector("list", length = d)
  for(j in 1:d){
    Pred_median[[j]] = Pred_all[[j]][2,]
    take_indices = idx_start[j]:idx_end[j]
    local_sim_matrix = psm(fit$Partition[,take_indices])
    local_coclus_list[[j]] = Compute_coclust_error( real_partition[take_indices],local_sim_matrix )
  }

  # Compute L1 distance wrt to true density, which is the same for all first d1 groups
  L1_list1 = Compute_L1_dist(Pred = Pred_median[1:d1],
                             p_mix = mix_probs1, mu = mu1, sigma = sd1,
                             grid = grid)
  # Compute L1 distance wrt to true density, which is the same for all final d-d1 groups
  L1_list2 = Compute_L1_dist(Pred = Pred_median[(1+d1):d],
                             p_mix = mix_probs2, mu = mu2, sigma = sd2,
                             grid = grid)


  # save
  results$HDP$ARI_est_part = ARI_VI
  results$HDP$ARI_est_part_local = ARI_j
  results$HDP$K_mode = K_mod
  results$HDP$K_mode_local = Kj_mod
  results$HDP$K_ARI = K_ARI
  results$HDP$K_ARI_local = Kj_ARI
  results$HDP$err_coclust = coclus_list$coclust_err
  results$HDP$err_coclust_local = local_coclus_j
  results$HDP$err_L1_local = c(L1_list1[[1]],L1_list2[[1]])
  results$HDP$err_L1_mean = (d1 * L1_list1[[2]] + d2 * L1_list2[[2]])/d

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
  sigma0 = 0.5

  scale = sqrt( (k0 + 1)/(k0) * sigma0 )
  mean_marginal = mu0
  var_marginal  = nu0/(nu0-2) * scale^2

  ### Process
  Exp_Lambda   =  15
  Var_Lambda   =  3
  gamma_guess  =  0.05
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
  ## Run
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
    ARI_VI  = arandi(VI_sara$cl,real_partition, adjust = FALSE)
    if(ARI_VI != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }

  ## Level dependent clustering
  Local_Clustering = list("list", length = d)
  idx_start = c(1,cumsum(n_j)[1:(d-1)]+1)
  idx_end = cumsum(n_j)
  Kj = matrix(0,nrow = niter, ncol = d)
  ARI_j = vector(length = d)
  Kj_ARI = vector(length = d)
  for(j in 1:d){
    Local_Clustering[[j]] = GDFMM$Partition[ , idx_start[j]:idx_end[j] ]
    Kj[,j] = apply(Local_Clustering[[j]], 1, FUN = function(Part_it){length(table(Part_it))})
    ARI_j[j]  = arandi(  VI_sara$cl[idx_start[j]:idx_end[j]],
                         real_partition[idx_start[j]:idx_end[j]],
                         adjust = TRUE  )
    if(is.na(ARI_j[j])){
      ARI_j[j]  = arandi(  VI_sara$cl[idx_start[j]:idx_end[j]],
                           real_partition[idx_start[j]:idx_end[j]],
                           adjust = FALSE  )
      if(ARI_j[j] != 1)
        stop("Adjusted RI is NaN but Unadjusted RI is not 1")
    }
    Kj_ARI[j] = length(table(VI_sara$cl[idx_start[j]:idx_end[j]]))

  }
  K_mod  = as.numeric(names(which.max(table(GDFMM$K))))
  Kj_mod = unlist(lapply(1:d, function(jj){as.numeric(names(which.max(table(Kj[,jj]))))}))
  K_ARI  = length(table(VI_sara$cl))
  ## Density estimation
  xrange = c(-7,7)
  l_grid = 200
  grid = seq(xrange[1],xrange[2],length.out = l_grid)
  # Predictive in all groups
  Pred_all = GDFMM::predictive_all_groups(grid = grid, fit = GDFMM, burnin = 1)
  ## Indici
  coclus_list = Compute_coclust_error(real_partition, sim_matrix)
  local_coclus_list = vector("list", length = d)
  Pred_median = vector("list", length = d)
  for(j in 1:d){
    Pred_median[[j]] = Pred_all[[j]][2,]
    take_indices = idx_start[j]:idx_end[j]
    local_sim_matrix = psm(GDFMM$Partition[,take_indices])
    local_coclus_list[[j]] = Compute_coclust_error( real_partition[take_indices],local_sim_matrix )
  }

  # Compute L1 distance wrt to true density, which is the same for all first d1 groups
  L1_list1 = Compute_L1_dist(Pred = Pred_median[1:d1],
                             p_mix = mix_probs1, mu = mu1, sigma = sd1,
                             grid = grid)
  # Compute L1 distance wrt to true density, which is the same for all final d-d1 groups
  L1_list2 = Compute_L1_dist(Pred = Pred_median[(1+d1):d],
                             p_mix = mix_probs2, mu = mu2, sigma = sd2,
                             grid = grid)


  # save
  results$HMFM_cond$ARI_est_part = ARI_VI
  results$HMFM_cond$ARI_est_part_local = ARI_j
  results$HMFM_cond$K_mode = K_mod
  results$HMFM_cond$K_mode_local = Kj_mod
  results$HMFM_cond$K_ARI = K_ARI
  results$HMFM_cond$K_ARI_local = Kj_ARI
  results$HMFM_cond$err_coclust = coclus_list$coclust_err
  results$HMFM_cond$err_coclust_local = local_coclus_j
  results$HMFM_cond$err_L1_local = c(L1_list1[[1]],L1_list2[[1]])
  results$HMFM_cond$err_L1_mean = (d1 * L1_list1[[2]] + d2 * L1_list2[[2]])/d

  # remove and finish
  rm(GDFMM, sim_matrix, Kj, Pred_all,Local_Clustering)

  # d) MFM - pooled ---------------------------------------------------------
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
    ARI_VI  = arandi(VI_sara$cl,real_partition, adjust = FALSE)
    if(ARI_VI != 1)
      stop("Adjusted RI is NaN but Unadjusted RI is not 1")
  }
  idx_start = c(1,cumsum(n_j)[1:(d-1)]+1)
  idx_end = cumsum(n_j)
  ARI_j = vector(length = d)
  for(j in 1:d){
    ARI_j[j]  = arandi(  VI_sara$cl[idx_start[j]:idx_end[j]],
                         real_partition[idx_start[j]:idx_end[j]],
                         adjust = TRUE  )
    if(is.na(ARI_j[j])){
      ARI_j[j]  = arandi(  VI_sara$cl[idx_start[j]:idx_end[j]],
                           real_partition[idx_start[j]:idx_end[j]],
                           adjust = FALSE  )
      if(ARI_j[j] != 1)
        stop("Adjusted RI is NaN but Unadjusted RI is not 1")
    }

  }
  K_mod  = as.numeric(names(which.max(table(fit$K))))
  K_ARI  = length(table(VI_sara$cl))


  ## Indici
  coclus_list = Compute_coclust_error(real_partition, sim_matrix)
  local_coclus_list = vector("list", length = d)
  for(j in 1:d){
    take_indices = idx_start[j]:idx_end[j]
    local_coclus_list[[j]] = NaN#Compute_coclust_error( real_partition[take_indices],sim_matrix[take_indices,take_indices] )
  }
  local_coclus_j = NaN#sapply(local_coclus_list, FUN = function(xx){xx[[1]]})

  # save
  results$pooled$ARI_est_part = ARI_VI
  results$pooled$ARI_est_part_local = ARI_j
  results$pooled$K_mode = K_mod
  results$pooled$K_mode_local = NaN
  results$pooled$K_ARI = K_ARI
  results$pooled$K_ARI_local = NaN
  results$pooled$err_coclust = coclus_list$coclust_err
  results$pooled$err_coclust_local = NaN#local_coclus_j
  results$pooled$err_L1_local = NaN
  results$pooled$err_L1_mean = NaN


  # remove and finish
  rm(fit, sim_matrix)


  # return ------------------------------------------------------------------
  return(results)

}


# Example plot -----------------------------------------------------------------
d = 15
K = 5
n_j = rep(30,d)
d1 = 12
d2 = d-d1
seed = 123545
# prima mistura
mu1 = c(-3,0,3)         # vectors of means
sd1 = c(sqrt(0.5),sqrt(0.5),sqrt(0.5))       # vector of sd
mix_probs1_vec = c(0.25,0.5,0.25)
mix_probs1 = matrix(0,nrow = d1, ncol = length(mu1))
K_j = rep(0,d)
set.seed(seed)
for(j in 1:d1){
  K_j[j] = sample(2:3, size = 1)
  components_j = sort(sample(1:3, size = K_j[j]))
  mix_probs_j = rep(0,3)
  mix_probs_j[components_j] = mix_probs1_vec[components_j]
  mix_probs_j = mix_probs_j/sum(mix_probs_j)
  mix_probs1[j,] = mix_probs_j
}
genD1 = simulate_data(d = d1, K = 3, p = mix_probs1,
                      mu = mu1, sd = sd1, n_j = n_j[1:d1],
                      seed = seed)
data1 = genD1$data
real_partition1 = genD1$real_partition


# seconda mistura
mu2 = c(-1.5,1.5)         # vectors of means
sd2 = c(sqrt(0.5),sqrt(0.5))       # vector of sd
mix_probs2_vec = c(0.5,0.5)
mix_probs2 = matrix(0,nrow = d2, ncol = length(mu2))
K_j[(d1+1):d] = 2
for(j in 1:d2){
  mix_probs2[j,] = mix_probs2_vec
}

genD2 = simulate_data(d = d2, K = 2, p = mix_probs2,
                      mu = mu2, sd = sd2, n_j = n_j[(d1+1):d],
                      seed = seed)
data2 = genD2$data
real_partition2 = genD2$real_partition

data = rbind(data1,data2)
real_partition = c(real_partition1,real_partition2)

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
yrange = c(0,1)
l_grid = 200
grid = seq(xrange[1],xrange[2],length.out = l_grid)



par( mfrow = c(3,5), mar = c(2,2,2,1), bty = "l")
for(j in 1:d1){
  plot(0,0,main = paste0("Group ",j), xlab = " ", type = "n", xlim = xrange, ylim = yrange)
  grid(lty = 1,lwd = 1, col = "gray90" )
  for(k in 1:K){
    res = data_small %>% filter(Group == j) %>% filter(TrueClustering == k) %>% pull(Value)
    if(length(res)>0){
      hist(res, freq = FALSE, breaks = 4, #nclass = "fd",
           col= ACutils::t_col(mycol_cluster[k], percent = 60), add = T)
      points( x = res, y = rep(0,length(res) ), pch = 16, col = mycol_cluster[k] )
    }
  }
  points(grid, GDFMM::dmix(x = grid, w_j = mix_probs1[j,], mu_vec = mu1, sigma_vec = sd1),
         col = "red", lwd = 2, type = "l")

}
for(j in (d1+1):d){
  #par( mar = c(2,2,2,1), bty = "l")
  plot(0,0,main = paste0("Group ",j), xlab = " ", type = "n", xlim = xrange, ylim = yrange)
  grid(lty = 1,lwd = 1, col = "gray90" )
  for(k in 1:K){
    res = data_small %>% filter(Group == j) %>% filter(TrueClustering == k) %>% pull(Value)
    if(length(res)>0){
      hist(res, freq = FALSE, breaks = 4, #nclass = "fd",
           col= ACutils::t_col(mycol_cluster[3+k], percent = 60), add = T)
      points( x = res, y = rep(0,length(res) ), pch = 16, col = mycol_cluster[3+k] )
    }
  }
  points(grid, GDFMM::dmix(x = grid, w_j = mix_probs2[j-d1,], mu_vec = mu2, sigma_vec = sd2),
         col = "red", lwd = 2, type = "l")

}


# Pooled data
par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
plot(0,0,main = paste0("Pooled data"),xlab = " ", type = "n", xlim = xrange, ylim = c(0,0.75))
grid(lty = 1,lwd = 1, col = "gray90" )
for(k in 1:3){
  res = data_small %>% filter(Group %in% 1:d1) %>% filter(TrueClustering == k) %>% pull(Value)
  if(length(res)>0){
    hist(res, freq = FALSE, breaks = 4, #nclass = "fd",
         col= ACutils::t_col(mycol_cluster[k], percent = 60), add = T)
    points( x = res, y = rep(0,length(res) ), pch = 16, col = mycol_cluster[k] )
  }
}
for(k in 4:5){
  res = data_small %>% filter(Group %in% (d1+1):d) %>% filter(TrueClustering == k-3) %>% pull(Value)
  if(length(res)>0){
    hist(res, freq = FALSE, breaks = 4, #nclass = "fd",
         col= ACutils::t_col(mycol_cluster[k], percent = 60), add = T)
    points( x = res, y = rep(0,length(res) ), pch = 16, col = mycol_cluster[k] )
  }
}

dens_pooled = density(data_small %>% pull(Value))
points(dens_pooled$x, dens_pooled$y, col = "red", lwd = 2, type = "l")




# Run ----------------------------------------------------------

Nrep  = 50

seed0 = 290696
set.seed(seed0)
seeds = sample(1:999999, size = Nrep)
num_cores = 3

tictoc::tic()
  cluster <- parallel::makeCluster(num_cores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)
  parallel::clusterExport(cluster, list())
  res = parallel::parLapply( cl = cluster, seeds,
                             fun = SimStudy_Exp3)
  parallel::stopCluster(cluster)
tictoc::toc()

beepr::beep()

save(res,file = "Exp3.Rdat")

## Predictive score
name = "err_L1_local"
ylabel = "PS group "
PS_plot = vector("list", length = d)
for(j in 1:d){
  HDP_res = sapply(res, function(x){x$HDP$err_L1_local[j]})
  HMFMmarg_res = sapply(res, function(x){x$HMFM_marg$err_L1_local[j]})
  HMFMcond_res = sapply(res, function(x){x$HMFM_cond$err_L1_local[j]})

  exp_temp = tibble("var" = HMFMcond_res, "type" = as_factor("HMFM-cond"))
  exp_temp = exp_temp %>%
    rbind(tibble("var" = HMFMmarg_res, "type" = as_factor("HMFM-marg"))) %>%
    rbind(tibble("var" = HDP_res, "type" = as_factor("HDP")))


  PS_plot[[j]] = exp_temp %>% select(type,var) %>%
    ggplot(aes(y=var, x=type, fill=type)) + geom_boxplot(fill = col_type[1:3]) +
    labs(y=paste0(ylabel,j,", n = ",sum(n_j)), x = " ") + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), legend.position="none",
          text = element_text(size = 10))  + ylim(c(0,0.6))
}

# Combine all plots together
combined_plot <- do.call(ggpubr::ggarrange, c(PS_plot, nrow = 3, ncol = 5))

# Display the combined plot
print(combined_plot)


## ARI
name = "ARI_est_part"
ylabel = "ARI"

HDP_res = sapply(res, function(x){x$HDP$ARI_est_part})
HMFMmarg_res = sapply(res, function(x){x$HMFM_marg$ARI_est_part})
HMFMcond_res = sapply(res, function(x){x$HMFM_cond$ARI_est_part})
pooled_res = sapply(res, function(x){x$pooled$ARI_est_part})


exp_temp = tibble("ARI_est_part" = HMFMcond_res, "type" = as_factor("HMFM-cond"))
exp_temp = exp_temp %>%
  rbind(tibble("ARI_est_part" = HMFMmarg_res, "type" = as_factor("HMFM-marg"))) %>%
  rbind(tibble("ARI_est_part" = HDP_res, "type" = as_factor("HDP"))) %>%
  rbind(tibble("ARI_est_part" = pooled_res,   "type" = as_factor("MFM-pooled")))

ARI_plot = exp_temp %>% select(type,!!name) %>%
  ggplot(aes(y=!!sym(name), x=type, fill=type)) + geom_boxplot(fill = col_type) +
  labs(y=paste0(ylabel,", n = ",sum(n_j)), x = " ") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        text = element_text(size = 10))


## CCE
name = "err_coclust"
ylabel = "CCE"

HDP_res = sapply(res, function(x){x$HDP$err_coclust})
HMFMmarg_res = sapply(res, function(x){x$HMFM_marg$err_coclust})
HMFMcond_res = sapply(res, function(x){x$HMFM_cond$err_coclust})
pooled_res = sapply(res, function(x){x$pooled$err_coclust})


exp_temp = tibble("err_coclust" = HMFMcond_res, "type" = as_factor("HMFM-cond"))
exp_temp = exp_temp %>%
  rbind(tibble("err_coclust" = HMFMmarg_res, "type" = as_factor("HMFM-marg"))) %>%
  rbind(tibble("err_coclust" = HDP_res, "type" = as_factor("HDP"))) %>%
  rbind(tibble("err_coclust" = pooled_res,   "type" = as_factor("MFM-pooled")))

CCE_plot = exp_temp %>% select(type,!!name) %>%
  ggplot(aes(y=!!sym(name), x=type, fill=type)) + geom_boxplot(fill = col_type) +
  labs(y=paste0(ylabel,", n = ",sum(n_j)), x = " ") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        text = element_text(size = 10))


## K_ARI
name = "K_ARI"
ylabel = "Est. K"

HDP_res = sapply(res, function(x){x$HDP$K_ARI})
HMFMmarg_res = sapply(res, function(x){x$HMFM_marg$K_ARI})
HMFMcond_res = sapply(res, function(x){x$HMFM_cond$K_ARI})
pooled_res = sapply(res, function(x){x$pooled$K_ARI})

exp_temp = tibble("K_ARI" = HMFMcond_res, "type" = as_factor("HMFM-cond"))
exp_temp = exp_temp %>%
  rbind(tibble("K_ARI" = HMFMmarg_res, "type" = as_factor("HMFM-marg"))) %>%
  rbind(tibble("K_ARI" = HDP_res, "type" = as_factor("HDP"))) %>%
  rbind(tibble("K_ARI" = pooled_res,   "type" = as_factor("MFM-pooled")))


K_ARI_plot = ggplot(exp_temp, aes(x = K_ARI, fill = type)) +
  geom_bar(position = "dodge") + theme_bw() +
  scale_fill_manual(values = c("HDP" = col_type[3],
                               "HMFM-marg" = col_type[2],
                               "HMFM-cond" = col_type[1],
                               "MFM-pooled" = col_type[4])) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        text = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(min(exp_temp$K_ARI), max(exp_temp$K_ARI), by = 1)) +
  labs(y=paste0(ylabel,", n = ",sum(n_j)), x = " ")
# K_ARI_plot

beepr::beep()
