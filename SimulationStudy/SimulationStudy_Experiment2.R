# Libraries

suppressWarnings(suppressPackageStartupMessages(library(GDFMM)))
suppressWarnings(suppressPackageStartupMessages(library(hdp)))
suppressWarnings(suppressPackageStartupMessages(library(ACutils)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(library(salso)))
suppressWarnings(suppressPackageStartupMessages(library(wesanderson)))
suppressWarnings(suppressPackageStartupMessages(library(mcclust.ext)))
suppressWarnings(suppressPackageStartupMessages(library(abind)))



# n_j = c(25,25) ----------------------------------------------------------

Nrep =  50
seed0 = 1605

save_data = vector("list",length = Nrep)
save_res  =  tibble("type" = as.character(),  "time" = as.numeric(),
                    "ARI_est_part" = as.numeric(),"ARI_est_part_group1" = as.numeric(),"ARI_est_part_group2" = as.numeric(),
                    "K_mode" = as.integer(), "K1_mode" = as.integer(), "K2_mode" = as.integer(),
                    "K_ARI" = as.integer(), "K1_ARI" = as.integer(), "K2_ARI" = as.integer(),
                    "err_coclust" = as.numeric(), "err_coclust_star" = as.numeric(),
                    "err_coclust_group1" = as.numeric(), "err_coclust_star_group1" = as.numeric(),
                    "err_coclust_group2" = as.numeric(), "err_coclust_star_group2" = as.numeric(),
                    "err_L1_group1" = as.numeric(), "err_L1_group2" = as.numeric(),
                    "err_L1_mean" = as.numeric(),
                    "err_wass2_group1" = as.numeric(), "err_wass2_group2" = as.numeric(),
                    "err_wass2_mean" = as.numeric() )

for(rep in 1:Nrep){
  cat("\n ----------------------------- \n")
  cat(" rep =  ",rep)
  cat("\n ----------------------------- \n")
  #1) Data Generation
  d = 2                   # number of groups
  K = 2                   # number of global clusters
  mu=c(0,1)              # vectors of means
  sd = c(1,1)            # vector of sd
  n_j=c(25,25)           # set cardinality of the groups

  mix_probs = matrix(c(1,0,0,1), nrow = d, ncol = K, byrow = T)

  seed = seed0 * rep
  genD = simulate_data(d = d, K = K, p = mix_probs, mu = mu, sd = sd, n_j = n_j, seed = seed)
  data = genD$data
  real_partition = genD$real_partition
  data_list = vector("list",length = d)
  for(j in 1:d)
    data_list[[j]] = genD$data[j,1:n_j[j]]

  save_data[[rep]] = list("data" = data, "real_partition" = real_partition, "seed" = seed)

  #2) Vec-FDP (marginal)
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
  Exp_Lambda   = 10
  Var_Lambda   = 2
  gamma_guess  = 0.01
  Lambda_guess = Exp_Lambda
  b_lambda = Exp_Lambda/Var_Lambda
  a_lambda = Exp_Lambda * b_lambda
  a_gamma = a_lambda/d
  b_gamma = a_gamma / (gamma_guess * Lambda_guess)

  ## Run
  niter  <-  10000
  burnin <-   5000
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
  ARI_VI_1  = arandi(VI_sara$cl[1:n_j[1]], real_partition[1:n_j[1]], adjust = TRUE)
  ARI_VI_2  = arandi(VI_sara$cl[(n_j[1]+1):sum(n_j)],real_partition[(n_j[1]+1):sum(n_j)], adjust = TRUE)

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


  save_res_temp  = tibble("type" = "VecFDP", "time" = time,
                          "ARI_est_part" = ARI_VI,"ARI_est_part_group1" =ARI_VI_1,"ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = K_mod, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = K_ARI, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = coclus_list$coclust_err, "err_coclust_star" = coclus_list$coclust_err_star,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list$L1err_per_level[1], "err_L1_group2" = L1_list$L1err_per_level[2],
                          "err_L1_mean" = L1_list$L1err_average,
                          "err_wass2_group1" = wass2$wass2_per_level[1], "err_wass2_group2" = wass2$wass2_per_level[2],
                          "err_wass2_mean" = wass2$wass2_average )

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

  Niter   = niter
  Nburnin =  burnin

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

  # wasserstein 2
  wass2 = vector("list",2)
  wass2[[1]] = rep(0,d)
  for(j in 1:d){
    f1_true = GDFMM::dmix(x = grid, w_j = mix_probs[j,], mu_vec = mu, sigma_vec = sd)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f1_true, p = 2)
  }
  wass2[[2]] = mean(wass2[[1]])
  names(wass2) = c("wass2_per_level","wass2_average")


  save_res_temp  = tibble("type" = "HDP", "time" = time,
                          "ARI_est_part" = ARI_VI,"ARI_est_part_group1" =ARI_VI_1,"ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = K_mod, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = K_ARI, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = coclus_list$coclust_err, "err_coclust_star" = coclus_list$coclust_err_star,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list$L1err_per_level[1], "err_L1_group2" = L1_list$L1err_per_level[2],
                          "err_L1_mean" = L1_list$L1err_average,
                          "err_wass2_group1" = wass2$wass2_per_level[1], "err_wass2_group2" = wass2$wass2_per_level[2],
                          "err_wass2_mean" = wass2$wass2_average )

  save_res = save_res %>% rbind(save_res_temp)
  rm(fit, sim_matrix, Kj, Pred_all, Local_Clustering)

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
  sigma0 = 0.5

  scale = sqrt( (k0 + 1)/(k0) * sigma0 )
  mean_marginal = mu0
  var_marginal  = nu0/(nu0-2) * scale^2

  ### Process
  Exp_Lambda   = 10
  Var_Lambda   = 2
  gamma_guess  = 0.01

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

  # wasserstein 2
  wass2 = vector("list",2)
  wass2[[1]] = rep(0,d)
  for(j in 1:d){
    f1_true = GDFMM::dmix(x = grid, w_j = mix_probs[j,], mu_vec = mu, sigma_vec = sd)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f1_true, p = 2)
  }
  wass2[[2]] = mean(wass2[[1]])
  names(wass2) = c("wass2_per_level","wass2_average")


  save_res_temp  = tibble("type" = "VecFDP-Cond", "time" = time,
                          "ARI_est_part" = ARI_VI,"ARI_est_part_group1" =ARI_VI_1,"ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = K_mod, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = K_ARI, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = coclus_list$coclust_err, "err_coclust_star" = coclus_list$coclust_err_star,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list$L1err_per_level[1], "err_L1_group2" = L1_list$L1err_per_level[2],
                          "err_L1_mean" = L1_list$L1err_average,
                          "err_wass2_group1" = wass2$wass2_per_level[1], "err_wass2_group2" = wass2$wass2_per_level[2],
                          "err_wass2_mean" = wass2$wass2_average )

  save_res = save_res %>% rbind(save_res_temp)
  rm(GDFMM, sim_matrix, Kj, Pred_all,Local_Clustering)


}


save_res$type = as.factor(save_res$type)
save(save_data, file = "data_exp1_n25_newmarginal.Rdat")
save(save_res,  file = "res_exp1_n25_newmarginal.Rdat")






# n_j = c(50,50) ----------------------------------------------------------

Nrep =  50
seed0 = 1605

save_data = vector("list",length = Nrep)
save_res  =  tibble("type" = as.character(), "ARI_est_part" = as.numeric(), "time" = as.numeric(),
                    "K_mode" = as.integer(), "K1_mode" = as.integer(), "K2_mode" = as.integer(),
                    "err_coclust" = as.numeric(), "err_coclust_star" = as.numeric(),
                    "err_coclust_group1" = as.numeric(), "err_coclust_star_group1" = as.numeric(),
                    "err_coclust_group2" = as.numeric(), "err_coclust_star_group2" = as.numeric(),
                    "err_L1_group1" = as.numeric(), "err_L1_group2" = as.numeric(),
                    "err_L1_mean" = as.numeric(),
                    "err_wass2_group1" = as.numeric(), "err_wass2_group2" = as.numeric(),
                    "err_wass2_mean" = as.numeric() )

for(rep in 1:Nrep){
  cat("\n ----------------------------- \n")
  cat(" rep =  ",rep)
  cat("\n ----------------------------- \n")
  #1) Data Generation
  d = 2                   # number of groups
  K = 2                   # number of global clusters
  mu=c(0,1)              # vectors of means
  sd = c(1,1)            # vector of sd
  n_j=c(50,50)           # set cardinality of the groups

  mix_probs = matrix(c(1,0,0,1), nrow = d, ncol = K, byrow = T)

  seed = seed0 * rep
  genD = simulate_data(d = d, K = K, p = mix_probs, mu = mu, sd = sd, n_j = n_j, seed = seed)
  data = genD$data
  real_partition = genD$real_partition
  data_list = vector("list",length = d)
  for(j in 1:d)
    data_list[[j]] = genD$data[j,1:n_j[j]]

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
  Exp_Lambda   = 10
  Var_Lambda   = 2
  gamma_guess  = 0.01
  Lambda_guess = Exp_Lambda
  b_lambda = Exp_Lambda/Var_Lambda
  a_lambda = Exp_Lambda * b_lambda
  a_gamma = a_lambda/d
  b_gamma = a_gamma / (gamma_guess * Lambda_guess)

  ## Run
  niter  <-  10000
  burnin <-   5000
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
  K_mod = as.numeric(which.max(table(GDFMM$K)))
  K1_mod = as.numeric(which.max(table(Kj[,1])))
  K2_mod = as.numeric(which.max(table(Kj[,2])))

  ## Density estimation
  xrange = c(-5,5)
  l_grid = 200
  grid = seq(xrange[1],xrange[2],length.out = l_grid)
  Pred_all = predictive_marginal_all_groups(grid = grid, fit = GDFMM, burnin = 0, option = option)


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


  save_res_temp  = tibble("type" = "VecFDP", "time" = time,
                          "ARI_est_part" = ARI_VI,"ARI_est_part_group1" =ARI_VI_1,"ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = K_mod, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = K_ARI, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = coclus_list$coclust_err, "err_coclust_star" = coclus_list$coclust_err_star,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list$L1err_per_level[1], "err_L1_group2" = L1_list$L1err_per_level[2],
                          "err_L1_mean" = L1_list$L1err_average,
                          "err_wass2_group1" = wass2$wass2_per_level[1], "err_wass2_group2" = wass2$wass2_per_level[2],
                          "err_wass2_mean" = wass2$wass2_average )

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

  Niter   = niter
  Nburnin =  burnin

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

  # wasserstein 2
  wass2 = vector("list",2)
  wass2[[1]] = rep(0,d)
  for(j in 1:d){
    f1_true = GDFMM::dmix(x = grid, w_j = mix_probs[j,], mu_vec = mu, sigma_vec = sd)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f1_true, p = 2)
  }
  wass2[[2]] = mean(wass2[[1]])
  names(wass2) = c("wass2_per_level","wass2_average")


  save_res_temp  = tibble("type" = "HDP", "time" = time,
                          "ARI_est_part" = ARI_VI,"ARI_est_part_group1" =ARI_VI_1,"ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = K_mod, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = K_ARI, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = coclus_list$coclust_err, "err_coclust_star" = coclus_list$coclust_err_star,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list$L1err_per_level[1], "err_L1_group2" = L1_list$L1err_per_level[2],
                          "err_L1_mean" = L1_list$L1err_average,
                          "err_wass2_group1" = wass2$wass2_per_level[1], "err_wass2_group2" = wass2$wass2_per_level[2],
                          "err_wass2_mean" = wass2$wass2_average )

  save_res = save_res %>% rbind(save_res_temp)
  rm(fit, sim_matrix, Kj, Pred_all, Local_Clustering)

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
  sigma0 = 0.5

  scale = sqrt( (k0 + 1)/(k0) * sigma0 )
  mean_marginal = mu0
  var_marginal  = nu0/(nu0-2) * scale^2

  ### Process
  Exp_Lambda   = 10
  Var_Lambda   =  3
  gamma_guess  = 0.01

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

  # wasserstein 2
  wass2 = vector("list",2)
  wass2[[1]] = rep(0,d)
  for(j in 1:d){
    f1_true = GDFMM::dmix(x = grid, w_j = mix_probs[j,], mu_vec = mu, sigma_vec = sd)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f1_true, p = 2)
  }
  wass2[[2]] = mean(wass2[[1]])
  names(wass2) = c("wass2_per_level","wass2_average")


  save_res_temp  = tibble("type" = "VecFDP-Cond", "time" = time,
                          "ARI_est_part" = ARI_VI,"ARI_est_part_group1" =ARI_VI_1,"ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = K_mod, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = K_ARI, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = coclus_list$coclust_err, "err_coclust_star" = coclus_list$coclust_err_star,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list$L1err_per_level[1], "err_L1_group2" = L1_list$L1err_per_level[2],
                          "err_L1_mean" = L1_list$L1err_average,
                          "err_wass2_group1" = wass2$wass2_per_level[1], "err_wass2_group2" = wass2$wass2_per_level[2],
                          "err_wass2_mean" = wass2$wass2_average )

  save_res = save_res %>% rbind(save_res_temp)
  rm(GDFMM, sim_matrix, Kj, Pred_all,Local_Clustering)

}


save_res$type = as.factor(save_res$type)
# save(save_data, file = "data_exp1_n50_prova2.Rdat")
# save(save_res,  file = "res_exp1_n50_prova2.Rdat")
save(save_data, file = "data_exp1_n50_newmarginal.Rdat")
save(save_res,  file = "res_exp1_n50_newmarginal.Rdat")





# n_j = c(100,100) ----------------------------------------------------------

Nrep =  50
seed0 = 1605

save_data = vector("list",length = Nrep)
save_res  =  tibble("type" = as.character(), "ARI_est_part" = as.numeric(), "time" = as.numeric(),
                    "K_mode" = as.integer(), "K1_mode" = as.integer(), "K2_mode" = as.integer(),
                    "err_coclust" = as.numeric(), "err_coclust_star" = as.numeric(),
                    "err_coclust_group1" = as.numeric(), "err_coclust_star_group1" = as.numeric(),
                    "err_coclust_group2" = as.numeric(), "err_coclust_star_group2" = as.numeric(),
                    "err_L1_group1" = as.numeric(), "err_L1_group2" = as.numeric(),
                    "err_L1_mean" = as.numeric(),
                    "err_wass2_group1" = as.numeric(), "err_wass2_group2" = as.numeric(),
                    "err_wass2_mean" = as.numeric() )

for(rep in 1:Nrep){
  cat("\n ----------------------------- \n")
  cat(" rep =  ",rep)
  cat("\n ----------------------------- \n")
  #1) Data Generation
  d = 2                   # number of groups
  K = 2                   # number of global clusters
  mu=c(0,1)              # vectors of means
  sd = c(1,1)            # vector of sd
  n_j=c(100,100)           # set cardinality of the groups

  mix_probs = matrix(c(1,0,0,1), nrow = d, ncol = K, byrow = T)

  seed = seed0 * rep
  genD = simulate_data(d = d, K = K, p = mix_probs, mu = mu, sd = sd, n_j = n_j, seed = seed)
  data = genD$data
  real_partition = genD$real_partition
  data_list = vector("list",length = d)
  for(j in 1:d)
    data_list[[j]] = genD$data[j,1:n_j[j]]

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
  Exp_Lambda   = 10
  Var_Lambda   =  3
  gamma_guess  = 0.01
  Lambda_guess = Exp_Lambda
  b_lambda = Exp_Lambda/Var_Lambda
  a_lambda = Exp_Lambda * b_lambda
  a_gamma = a_lambda/d
  b_gamma = a_gamma / (gamma_guess * Lambda_guess)

  ## Run
  niter  <-  10000
  burnin <-   5000
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
  K_mod = as.numeric(which.max(table(GDFMM$K)))
  K1_mod = as.numeric(which.max(table(Kj[,1])))
  K2_mod = as.numeric(which.max(table(Kj[,2])))

  ## Density estimation
  xrange = c(-5,5)
  l_grid = 200
  grid = seq(xrange[1],xrange[2],length.out = l_grid)
  Pred_all = predictive_marginal_all_groups(grid = grid, fit = GDFMM, burnin = 0, option = option)


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


  save_res_temp  = tibble("type" = "VecFDP", "time" = time,
                          "ARI_est_part" = ARI_VI,"ARI_est_part_group1" =ARI_VI_1,"ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = K_mod, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = K_ARI, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = coclus_list$coclust_err, "err_coclust_star" = coclus_list$coclust_err_star,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list$L1err_per_level[1], "err_L1_group2" = L1_list$L1err_per_level[2],
                          "err_L1_mean" = L1_list$L1err_average,
                          "err_wass2_group1" = wass2$wass2_per_level[1], "err_wass2_group2" = wass2$wass2_per_level[2],
                          "err_wass2_mean" = wass2$wass2_average )

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

  Niter   = niter
  Nburnin =  burnin

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

  # wasserstein 2
  wass2 = vector("list",2)
  wass2[[1]] = rep(0,d)
  for(j in 1:d){
    f1_true = GDFMM::dmix(x = grid, w_j = mix_probs[j,], mu_vec = mu, sigma_vec = sd)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f1_true, p = 2)
  }
  wass2[[2]] = mean(wass2[[1]])
  names(wass2) = c("wass2_per_level","wass2_average")


  save_res_temp  = tibble("type" = "HDP", "time" = time,
                          "ARI_est_part" = ARI_VI,"ARI_est_part_group1" =ARI_VI_1,"ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = K_mod, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = K_ARI, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = coclus_list$coclust_err, "err_coclust_star" = coclus_list$coclust_err_star,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list$L1err_per_level[1], "err_L1_group2" = L1_list$L1err_per_level[2],
                          "err_L1_mean" = L1_list$L1err_average,
                          "err_wass2_group1" = wass2$wass2_per_level[1], "err_wass2_group2" = wass2$wass2_per_level[2],
                          "err_wass2_mean" = wass2$wass2_average )

  save_res = save_res %>% rbind(save_res_temp)
  rm(fit, sim_matrix, Kj, Pred_all, Local_Clustering)

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
  sigma0 = 0.5

  scale = sqrt( (k0 + 1)/(k0) * sigma0 )
  mean_marginal = mu0
  var_marginal  = nu0/(nu0-2) * scale^2

  ### Process
  Exp_Lambda   = 10
  Var_Lambda   =  3
  gamma_guess  = 0.01

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

  # wasserstein 2
  wass2 = vector("list",2)
  wass2[[1]] = rep(0,d)
  for(j in 1:d){
    f1_true = GDFMM::dmix(x = grid, w_j = mix_probs[j,], mu_vec = mu, sigma_vec = sd)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f1_true, p = 2)
  }
  wass2[[2]] = mean(wass2[[1]])
  names(wass2) = c("wass2_per_level","wass2_average")


  save_res_temp  = tibble("type" = "VecFDP-Cond", "time" = time,
                          "ARI_est_part" = ARI_VI,"ARI_est_part_group1" =ARI_VI_1,"ARI_est_part_group2" = ARI_VI_2,
                          "K_mode" = K_mod, "K1_mode" = K1_mod, "K2_mode" = K2_mod,
                          "K_ARI" = K_ARI, "K1_ARI" = K1_ARI, "K2_ARI" = K2_ARI,
                          "err_coclust" = coclus_list$coclust_err, "err_coclust_star" = coclus_list$coclust_err_star,
                          "err_coclust_group1" = coclus_list1$coclust_err, "err_coclust_star_group1" = coclus_list1$coclust_err_star,
                          "err_coclust_group2" = coclus_list2$coclust_err, "err_coclust_star_group2" = coclus_list2$coclust_err_star,
                          "err_L1_group1" = L1_list$L1err_per_level[1], "err_L1_group2" = L1_list$L1err_per_level[2],
                          "err_L1_mean" = L1_list$L1err_average,
                          "err_wass2_group1" = wass2$wass2_per_level[1], "err_wass2_group2" = wass2$wass2_per_level[2],
                          "err_wass2_mean" = wass2$wass2_average )

  save_res = save_res %>% rbind(save_res_temp)
  rm(GDFMM, sim_matrix, Kj, Pred_all,Local_Clustering)


}


save_res$type = as.factor(save_res$type)
#save(save_data, file = "data_exp1_n100_prova2.Rdat")
#save(save_res,  file = "res_exp1_n100_prova2.Rdat")
save(save_data, file = "data_exp1_n100_newmarginal.Rdat")
save(save_res,  file = "res_exp1_n100_newmarginal.Rdat")




