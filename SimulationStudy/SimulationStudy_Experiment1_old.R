# Libraries
suppressWarnings(suppressPackageStartupMessages(library(GDFMM)))
suppressWarnings(suppressPackageStartupMessages(library(hdp)))
suppressWarnings(suppressPackageStartupMessages(library(ACutils)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(library(salso)))
suppressWarnings(suppressPackageStartupMessages(library(AntMAN)))
suppressWarnings(suppressPackageStartupMessages(library(transport)))
suppressWarnings(suppressPackageStartupMessages(library(wesanderson)))
suppressWarnings(suppressPackageStartupMessages(library(mcclust.ext)))
suppressWarnings(suppressPackageStartupMessages(library(abind)))

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


Nrep =  50
seed0 = 1605

save_data = vector("list",length = Nrep)
save_res  = tibble("type" = as.character(), "time" = as.numeric(),
                   "ARI_est_part" = as.numeric(), "ARI_est_part_group1" = as.numeric(), "ARI_est_part_group2" = as.numeric(),
                   "K_mode" = as.integer(), "K1_mode" = as.integer(), "K2_mode" = as.integer(),
                   "K_ARI" = as.integer(), "K1_ARI" = as.integer(), "K2_ARI" = as.integer(),
                   "err_coclust" = as.numeric(), "err_coclust_star" = as.numeric(),
                   "err_coclust_group1" = as.numeric(), "err_coclust_star_group1" = as.numeric(),
                   "err_coclust_group2" = as.numeric(), "err_coclust_star_group2" = as.numeric(),
                   "err_L1_group1" = as.numeric(), "err_L1_group2" = as.numeric(), "err_L1_mean" = as.numeric(),
                   "err_wass2_group1" = as.numeric(), "err_wass2_group2" = as.numeric(), "err_wass2_mean" = as.numeric())


for(rep in 1:Nrep){
  cat("\n ----------------------------- \n")
  cat(" rep =  ",rep)
  cat("\n ----------------------------- \n")
  #1) Data Generation
  d = 2                   # number of groups
  K = 3                   # number of global clusters
  mu=c(-3,0,1.75)         # vectors of means
  sd = c(sqrt(0.1),sqrt(0.5),sqrt(1.5))       # vector of sd
  n_j=c(300,300)           # set cardinality of the groups
  mix_probs = matrix(c(0.5,0.5,0,
                       0,0.20,0.80), nrow = d, ncol = K, byrow = T)

  mix_obs = mix_probs * n_j

  seed = seed0 * rep

  # qua inutile, serve solo per definire le strutture con le dimensioni giuste
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

  save_data[[rep]] = list("data" = data, "real_partition" = real_partition, "seed" = seed)

  #2) Vec-FDP
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
  niter  <-  30000
  burnin <-  10000
  thin   <- 1
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

  tictoc::tic()
    GDFMM = GDFMM_marginal_sampler(data, niter, burnin, thin, seed = seed, FixPartition = F, option = option)
  time = tictoc::toc()
  time = time$toc - time$tic

  ## Clustering
  part_matrix <- GDFMM$Partition
  sim_matrix <- psm(part_matrix)

  VI_sara = minVI(sim_matrix)
  ARI_VI  = arandi(VI_sara$cl,real_partition, adjust = TRUE)
  ARI_VI_1  = arandi(VI_sara$cl[1:n_j[1]],real_partition[1:n_j[1]], adjust = TRUE)
  ARI_VI_2  = arandi(VI_sara$cl[(n_j[1]+1):sum(n_j)],real_partition[(n_j[1]+1):sum(n_j)], adjust = TRUE)

  K_ARI = length(table(VI_sara$cl))
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
  K_mod = as.numeric(names(which.max(table(GDFMM$K)))) #as.numeric(which.max(table(GDFMM$K)))
  K1_mod = as.numeric(names(which.max(table(Kj[,1]))))
  K2_mod = as.numeric(names(which.max(table(Kj[,2]))))

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

  # wasserstein 2
  wass2 = vector("list",2)
  wass2[[1]] = rep(0,d)
  for(j in 1:d){
    f1_true = GDFMM::dmix(x = grid, w_j = mix_probs[j,], mu_vec = mu, sigma_vec = sd)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f1_true, p = 2)
  }
  wass2[[2]] = mean(wass2[[1]])
  names(wass2) = c("wass2_per_level","wass2_average")


  save_res_temp  = tibble("type" = "VecFDP", "ARI_est_part" = ARI_VI, "time" = time,
                          "ARI_est_part_group1" = ARI_VI_1, "ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = K_mod, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = K_ARI, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = coclus_list$coclust_err, "err_coclust_star" = coclus_list$coclust_err_star,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list$L1err_per_level[1], "err_L1_group2" = L1_list$L1err_per_level[2],
                          "err_L1_mean" = L1_list$L1err_average,
                          "err_wass2_group1" = wass2$wass2_per_level[1], "err_wass2_group2" = wass2$wass2_per_level[2],
                          "err_wass2_mean" = wass2$wass2_average, )

  save_res = save_res %>% rbind(save_res_temp)
  rm(GDFMM, sim_matrix, Kj, Pred_all,Local_Clustering)

  #3) HDP

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

  Niter   = 30000
  Nburnin = 10000

  tictoc::tic()
  fit = HDPMarginalSampler(Niter = Niter, Nburnin = Nburnin, d = d, n_j = n_j,
                           data_list = data_list,
                           priorMean = priorMean, priorA = priorA, priorB = priorB,
                           priorLambda = priorLambda,
                           a_gamma = a_gamma, b_gamma = b_gamma,
                           a_alpha = a_alpha, b_alpha = b_alpha,
                           alpha_init = alpha0, gamma_init = gamma0, UpdateConc = TRUE)
  time = tictoc::toc()
  time = time$toc - time$tic

  ## Clustering
  part_matrix <- fit$Partition
  sim_matrix <- psm(part_matrix)
  VI_sara = minVI(sim_matrix)
  ARI_VI = arandi(VI_sara$cl,real_partition)
  ARI_VI_1  = arandi(VI_sara$cl[1:n_j[1]],real_partition[1:n_j[1]], adjust = TRUE)
  ARI_VI_2  = arandi(VI_sara$cl[(n_j[1]+1):sum(n_j)],real_partition[(n_j[1]+1):sum(n_j)], adjust = TRUE)

  K_ARI = length(table(VI_sara$cl))
  K1_ARI = length(table(VI_sara$cl[1:n_j[1]]))
  K2_ARI = length(table(VI_sara$cl[(n_j[1]+1):sum(n_j)]))

  Local_Clustering = list("list", length = d)
  idx_start = c(1,cumsum(n_j)[1:(d-1)]+1)
  idx_end = cumsum(n_j)
  Kj = matrix(0,nrow = Niter, ncol = d)
  for(j in 1:d){
    Local_Clustering[[j]] = fit$Partition[ , idx_start[j]:idx_end[j] ]
    Kj[,j] = apply(Local_Clustering[[j]], 1, FUN = function(Part_it){length(table(Part_it))})
  }
  K_mod = as.numeric(names(which.max(table(fit$K))))
  K1_mod = as.numeric(names(which.max(table(Kj[,1]))))
  K2_mod = as.numeric(names(which.max(table(Kj[,2]))))


  ## Density estimation
  xrange = c(-5,5)
  l_grid = 200
  grid = seq(xrange[1],xrange[2],length.out = l_grid)
  # Predictive in all groups
  Pred_all = hdp::predictive_all_groups(grid = grid, fit = fit,
                                        priorMean = priorMean, priorA = priorA, priorB = priorB,
                                        priorLambda = priorLambda, burnin = 1)
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


  # wasserstein 2
  wass2 = vector("list",2)
  wass2[[1]] = rep(0,d)
  for(j in 1:d){
    f1_true = GDFMM::dmix(x = grid, w_j = mix_probs[j,], mu_vec = mu, sigma_vec = sd)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f1_true, p = 2)
  }
  wass2[[2]] = mean(wass2[[1]])
  names(wass2) = c("wass2_per_level","wass2_average")


  save_res_temp  = tibble("type" = "HDP", "ARI_est_part" = ARI_VI, "time" = time,
                          "ARI_est_part_group1" = ARI_VI_1, "ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = K_mod, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = K_ARI, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = coclus_list$coclust_err, "err_coclust_star" = coclus_list$coclust_err_star,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list$L1err_per_level[1], "err_L1_group2" = L1_list$L1err_per_level[2],
                          "err_L1_mean" = L1_list$L1err_average,
                          "err_wass2_group1" = wass2$wass2_per_level[1], "err_wass2_group2" = wass2$wass2_per_level[2],
                          "err_wass2_mean" = wass2$wass2_average, )

  save_res = save_res %>% rbind(save_res_temp)
  rm(fit, sim_matrix, Kj, Pred_all,Local_Clustering)

  #3) Conditional Vec-FDP
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
  ## Run
  niter  <-  30000
  burnin <-  10000
  thin   <-      1

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

  tictoc::tic()
  GDFMM = ConditionalSampler(dt[1:11], niter, burnin, thin, seed = 123, option = option, FixPartition = F,
                             P0.prior = prior)
  time = tictoc::toc()
  time = time$toc - time$tic

  ## Clustering
  part_matrix <- GDFMM$Partition
  sim_matrix <- psm(part_matrix)
  VI_sara = minVI(sim_matrix)
  ARI_VI = arandi(VI_sara$cl,real_partition)
  ARI_VI_1  = arandi(VI_sara$cl[1:n_j[1]],real_partition[1:n_j[1]], adjust = TRUE)
  ARI_VI_2  = arandi(VI_sara$cl[(n_j[1]+1):sum(n_j)],real_partition[(n_j[1]+1):sum(n_j)], adjust = TRUE)

  K_ARI  = length(table(VI_sara$cl))
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

  # wasserstein 2
  wass2 = vector("list",2)
  wass2[[1]] = rep(0,d)
  for(j in 1:d){
    f1_true = GDFMM::dmix(x = grid, w_j = mix_probs[j,], mu_vec = mu, sigma_vec = sd)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f1_true, p = 2)
  }
  wass2[[2]] = mean(wass2[[1]])
  names(wass2) = c("wass2_per_level","wass2_average")


  save_res_temp  = tibble("type" = "VecFDP-Cond", "ARI_est_part" = ARI_VI, "time" = time,
                          "ARI_est_part_group1" = ARI_VI_1, "ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = K_mod, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = K_ARI, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = coclus_list$coclust_err, "err_coclust_star" = coclus_list$coclust_err_star,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list$L1err_per_level[1], "err_L1_group2" = L1_list$L1err_per_level[2],
                          "err_L1_mean" = L1_list$L1err_average,
                          "err_wass2_group1" = wass2$wass2_per_level[1], "err_wass2_group2" = wass2$wass2_per_level[2],
                          "err_wass2_mean" = wass2$wass2_average, )

  save_res = save_res %>% rbind(save_res_temp)
  rm(GDFMM, sim_matrix, Kj, Pred_all,Local_Clustering)

  #5) AntMAN - Independent analysis

  # Primo gruppo
  ## Hyperparam

  ### $P_0$
  Range = range(data[1,1:n_j[1]])
  mu0 = mean(data[1,1:n_j[1]])
  R = Range[2] - Range[1]
  k0  = 1/R^2
  nu0 = 4
  sigma0 = 0.5
  scale = sqrt( (k0 + 1)/(k0) * sigma0 )
  mean_marginal = mu0
  var_marginal  = nu0/(nu0-2) * scale^2

  mixture_uvn_params = AntMAN::AM_mix_hyperparams_uninorm  (m0=mu0, k0=k0, nu0=nu0, sig02=sigma0)
  mcmc_params        = AntMAN::AM_mcmc_parameters(niter=niter, burnin=burnin, thin=10, verbose=1)
  components_prior   = AntMAN::AM_mix_components_prior_pois (init=3,  a=1, b=1)
  weights_prior      = AntMAN::AM_mix_weights_prior_gamma(init=2, a=1, b=1)

  tictoc::tic()
  fit <- AntMAN::AM_mcmc_fit(
    y = data[1,1:n_j[1]],
    mix_kernel_hyperparams = mixture_uvn_params,
    mix_components_prior =components_prior,
    mix_weight_prior = weights_prior,
    mcmc_parameters = mcmc_params)
  time1 = tictoc::toc()
  time1 = time1$toc - time1$tic

  ## Clustering
  sim_matrix <- AntMAN::AM_coclustering(fit)

  VI_sara = minVI(sim_matrix)
  ARI_VI_1  = arandi(VI_sara$cl,real_partition[1:n_j[1]], adjust = TRUE)

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
  for(j in 1:1){
    Pred_median[[j]] = Pred_all[2,]
  }

  L1_list1 = Compute_L1_dist(Pred = Pred_median,
                             p_mix = matrix(mix_probs[1,],nrow=1,ncol=3), mu = mu, sigma = sd,
                             grid = grid)

  wass2_1 = 0
  for(j in 1:1){
    f1_true = GDFMM::dmix(x = grid, w_j = mix_probs[j,], mu_vec = mu, sigma_vec = sd)
    wass2_1 = transport::wasserstein1d(a = Pred_median[[j]], b = f1_true, p = 2)
  }

  rm(fit, sim_matrix, Pred_all)

  # Secondo gruppo
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

  tictoc::tic()
  fit <- AntMAN::AM_mcmc_fit(
    y = data[2,1:n_j[2]],
    mix_kernel_hyperparams = mixture_uvn_params,
    mix_components_prior =components_prior,
    mix_weight_prior = weights_prior,
    mcmc_parameters = mcmc_params)
  time2 = tictoc::toc()
  time2 = time2$toc - time2$tic

  ## Clustering
  sim_matrix <- AntMAN::AM_coclustering(fit)

  VI_sara = minVI(sim_matrix)
  ARI_VI_2  =  arandi(VI_sara$cl,real_partition[(n_j[1]+1):sum(n_j)], adjust = TRUE)

  K2_ARI = length(table(VI_sara$cl))
  K2_mod = as.numeric(names(which.max(table(fit$K))))

  ## Density estimation
  xrange = c(-7,7)
  l_grid = 200
  grid = seq(xrange[1],xrange[2],length.out = l_grid)
  Pred_all = AM_density_estimation(grid = grid, fit = fit, burnin = 1)

  ## Indici
  coclus_list2 = Compute_coclust_error(real_partition[(n_j[1]+1):sum(n_j)], sim_matrix)
  Pred_median = vector("list", length = 1)
  for(j in 1:1){
    Pred_median[[j]] = Pred_all[2,]
  }

  L1_list2 = Compute_L1_dist(Pred = Pred_median,
                             p_mix = matrix(mix_probs[2,],nrow=1,ncol=3), mu = mu, sigma = sd,
                             grid = grid)

  wass2_2 = 0
  for(j in 1:1){
    f1_true = GDFMM::dmix(x = grid, w_j = mix_probs[2,], mu_vec = mu, sigma_vec = sd)
    wass2_2 = transport::wasserstein1d(a = Pred_median[[j]], b = f1_true, p = 2)
  }


  save_res_temp  = tibble("type" = "AntMAN", "ARI_est_part" = NA, "time" = time1+time2,
                          "ARI_est_part_group1" = ARI_VI_1, "ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = NA, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = NA, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = NA,
                          "err_coclust_star" = NA,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list1$L1err_per_level[1], "err_L1_group2" = L1_list2$L1err_per_level[1],
                          "err_L1_mean" = mean(L1_list1$L1err_average, L1_list2$L1err_average),
                          "err_wass2_group1" = wass2_1,"err_wass2_group2" = wass2_2, "err_wass2_mean" = mean(wass2_1,wass2_2) )

  save_res = save_res %>% rbind(save_res_temp)
  rm(fit, sim_matrix, Pred_all)
}


save_res$type = as.factor(save_res$type)
# save(save_data, file = "data_exp2_cond_indep_n300.Rdat")
# save(save_res,  file = "res_exp2_cond_indep_n300.Rdat")
#save(save_data, file = "data_exp2_cond_indep_n30.Rdat")
#save(save_res,  file = "res_exp2_cond_indep_n30.Rdat")

save(save_data, file = "data_exp2_n300_final_newmarginal.Rdat")
save(save_res,  file = "res_exp2_n300_final_newmarginal.Rdat")
