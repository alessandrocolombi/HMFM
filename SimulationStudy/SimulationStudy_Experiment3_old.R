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
#suppressWarnings(suppressPackageStartupMessages(library(AntMAN)))



Nrep =  50
seed0 = 1605

save_data = vector("list",length = Nrep)
save_res  = tibble("type" = as.character(), "time" = as.numeric(),
                   "ARI_est_part" = as.numeric(),
                   "ARI_group1" = as.numeric(), "ARI_group2" = as.numeric(),
                   "ARI_group3" = as.numeric(), "ARI_group4" = as.numeric(),
                   "ARI_group5" = as.numeric(), "ARI_group6" = as.numeric(),
                   "ARI_group7" = as.numeric(), "ARI_group8" = as.numeric(),
                   "ARI_group9" = as.numeric(), "ARI_group10" = as.numeric(),
                   "ARI_group11" = as.numeric(), "ARI_group12" = as.numeric(),
                   "ARI_group13" = as.numeric(), "ARI_group14" = as.numeric(),
                   "ARI_group15" = as.numeric(),
                   "K_mode" = as.integer(),
                   "K_mode_group1" = as.numeric(), "K_mode_group2" = as.numeric(),
                   "K_mode_group3" = as.numeric(), "K_mode_group4" = as.numeric(),
                   "K_mode_group5" = as.numeric(), "K_mode_group6" = as.numeric(),
                   "K_mode_group7" = as.numeric(), "K_mode_group8" = as.numeric(),
                   "K_mode_group9" = as.numeric(), "K_mode_group10" = as.numeric(),
                   "K_mode_group11" = as.numeric(), "K_mode_group12" = as.numeric(),
                   "K_mode_group13" = as.numeric(), "K_mode_group14" = as.numeric(),
                   "K_mode_group15" = as.numeric(),
                   "K_ARI" = as.integer(),
                   "K_ARI_group1" = as.numeric(), "K_ARI_group2" = as.numeric(),
                   "K_ARI_group3" = as.numeric(), "K_ARI_group4" = as.numeric(),
                   "K_ARI_group5" = as.numeric(), "K_ARI_group6" = as.numeric(),
                   "K_ARI_group7" = as.numeric(), "K_ARI_group8" = as.numeric(),
                   "K_ARI_group9" = as.numeric(), "K_ARI_group10" = as.numeric(),
                   "K_ARI_group11" = as.numeric(), "K_ARI_group12" = as.numeric(),
                   "K_ARI_group13" = as.numeric(), "K_ARI_group14" = as.numeric(),
                   "K_ARI_group15" = as.numeric(),
                   "err_coclust" = as.numeric(),
                   "err_coclust_group1" = as.numeric(), "err_coclust_group2" = as.numeric(),
                   "err_coclust_group3" = as.numeric(), "err_coclust_group4" = as.numeric(),
                   "err_coclust_group5" = as.numeric(), "err_coclust_group6" = as.numeric(),
                   "err_coclust_group7" = as.numeric(), "err_coclust_group8" = as.numeric(),
                   "err_coclust_group9" = as.numeric(), "err_coclust_group10" = as.numeric(),
                   "err_coclust_group11" = as.numeric(), "err_coclust_group12" = as.numeric(),
                   "err_coclust_group13" = as.numeric(), "err_coclust_group14" = as.numeric(),
                   "err_coclust_group15" = as.numeric(),
                   "err_L1_group1" = as.numeric(), "err_L1_group2" = as.numeric(),
                   "err_L1_group3" = as.numeric(), "err_L1_group4" = as.numeric(),
                   "err_L1_group5" = as.numeric(), "err_L1_group6" = as.numeric(),
                   "err_L1_group7" = as.numeric(), "err_L1_group8" = as.numeric(),
                   "err_L1_group9" = as.numeric(), "err_L1_group10" = as.numeric(),
                   "err_L1_group11" = as.numeric(), "err_L1_group12" = as.numeric(),
                   "err_L1_group13" = as.numeric(), "err_L1_group14" = as.numeric(),
                   "err_L1_group15" = as.numeric(),
                   "err_L1_mean" = as.numeric(),
                   "err_wass2_group1" = as.numeric(), "err_wass2_group2" = as.numeric(),
                   "err_wass2_group3" = as.numeric(), "err_wass2_group4" = as.numeric(),
                   "err_wass2_group5" = as.numeric(), "err_wass2_group6" = as.numeric(),
                   "err_wass2_group7" = as.numeric(), "err_wass2_group8" = as.numeric(),
                   "err_wass2_group9" = as.numeric(), "err_wass2_group10" = as.numeric(),
                   "err_wass2_group11" = as.numeric(), "err_wass2_group12" = as.numeric(),
                   "err_wass2_group13" = as.numeric(), "err_wass2_group14" = as.numeric(),
                   "err_wass2_group15" = as.numeric(),
                   "err_wass2_mean" = as.numeric()
)



for(rep in 1:Nrep){
  cat("\n ----------------------------- \n")
  cat(" rep =  ",rep)
  cat("\n ----------------------------- \n")
  seed = seed0 * rep

  #1) Data Generation
  d = 15
  K = 5
  n_j = rep(30,d)
  d1 = 12
  d2 = d-d1
  seed = seed0 * rep
  set.seed(seed)
  # prima mistura
  mu1 = c(-3,0,3)         # vectors of means
  sd1 = c(sqrt(0.5),sqrt(0.5),sqrt(0.5))       # vector of sd
  mix_probs1_vec = c(0.25,0.5,0.25)
  mix_probs1 = matrix(0,nrow = d1, ncol = length(mu1))
  K_j = rep(0,d)
  set.seed(123545)
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

  save_data[[rep]] = list("data" = data, "real_partition" = real_partition, "seed" = seed)
  set.seed(seed)
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
  Exp_Lambda   = 15
  Var_Lambda   =  3
  gamma_guess  =  0.05
  Lambda_guess = Exp_Lambda
  b_lambda = Exp_Lambda/Var_Lambda
  a_lambda = Exp_Lambda * b_lambda
  a_gamma = a_lambda/d
  b_gamma = a_gamma / (gamma_guess * Lambda_guess)

  ## Run
  niter  <-  35000
  burnin <-  25000
  thin   <- 1
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

  tictoc::tic()
    GDFMM = GDFMM_marginal_sampler(data, niter, burnin, thin, seed = seed, FixPartition = F, option = option)
  time = tictoc::toc()
  time = time$toc - time$tic

  ## Clustering
  part_matrix <- GDFMM$Partition
  sim_matrix <- psm(part_matrix)

  VI_sara = minVI(sim_matrix)
  ARI_VI  = arandi(VI_sara$cl,real_partition, adjust = TRUE)


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
    local_coclus_list[[j]] = Compute_coclust_error( real_partition[take_indices],sim_matrix[take_indices,take_indices] )
  }

  sd = sd1 # c'? un bug nella funzione, ho scritto sd invece che sigma. da sistemare!!!!
  L1_list1 = Compute_L1_dist(Pred = Pred_median[1:d1],
                             p_mix = mix_probs1, mu = mu1, sigma = sd1,
                             grid = grid)
  sd = sd2
  L1_list2 = Compute_L1_dist(Pred = Pred_median[(1+d1):d],
                             p_mix = mix_probs2, mu = mu2, sigma = sd2,
                             grid = grid)
  # wasserstein 2
  wass2 = vector("list",2)
  wass2[[1]] = rep(0,d)
  for(j in 1:d1){
    f_true = GDFMM::dmix(x = grid, w_j = mix_probs1[j,], mu_vec = mu1, sigma_vec = sd1)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f_true, p = 2)
  }
  for(j in (d1+1):d){
    f_true = GDFMM::dmix(x = grid, w_j = mix_probs2[j-d1,], mu_vec = mu2, sigma_vec = sd2)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f_true, p = 2)
  }
  wass2[[2]] = mean(wass2[[1]])
  names(wass2) = c("wass2_per_level","wass2_average")

  save_res_temp  = tibble("type" = "VecFDP", "time" = time,
                          "ARI_est_part" = ARI_VI,
                          "ARI_group1" = ARI_j[1], "ARI_group2" = ARI_j[2],
                          "ARI_group3" = ARI_j[3], "ARI_group4" = ARI_j[4],
                          "ARI_group5" = ARI_j[5], "ARI_group6" = ARI_j[6],
                          "ARI_group7" = ARI_j[7], "ARI_group8" = ARI_j[8],
                          "ARI_group9" = ARI_j[9], "ARI_group10" = ARI_j[10],
                          "ARI_group11" = ARI_j[11], "ARI_group12" = ARI_j[12],
                          "ARI_group13" = ARI_j[13], "ARI_group14" = ARI_j[14],
                          "ARI_group15" = ARI_j[15],
                          "K_mode" = K_mod,
                          "K_mode_group1" = Kj_mod[1], "K_mode_group2" = Kj_mod[2],
                          "K_mode_group3" = Kj_mod[3], "K_mode_group4" = Kj_mod[4],
                          "K_mode_group5" = Kj_mod[5], "K_mode_group6" = Kj_mod[6],
                          "K_mode_group7" = Kj_mod[7], "K_mode_group8" = Kj_mod[8],
                          "K_mode_group9" = Kj_mod[9], "K_mode_group10" = Kj_mod[10],
                          "K_mode_group11" = Kj_mod[11], "K_mode_group12" = Kj_mod[12],
                          "K_mode_group13" = Kj_mod[13], "K_mode_group14" = Kj_mod[14],
                          "K_mode_group15" = Kj_mod[15],
                          "K_ARI" = K_ARI,
                          "K_ARI_group1" = Kj_ARI[1], "K_ARI_group2" = Kj_ARI[2],
                          "K_ARI_group3" = Kj_ARI[3], "K_ARI_group4" = Kj_ARI[4],
                          "K_ARI_group5" = Kj_ARI[5], "K_ARI_group6" = Kj_ARI[6],
                          "K_ARI_group7" = Kj_ARI[7], "K_ARI_group8" = Kj_ARI[8],
                          "K_ARI_group9" = Kj_ARI[9], "K_ARI_group10" = Kj_ARI[10],
                          "K_ARI_group11" = Kj_ARI[11], "K_ARI_group12" = Kj_ARI[12],
                          "K_ARI_group13" = Kj_ARI[13], "K_ARI_group14" = Kj_ARI[14],
                          "K_ARI_group15" = Kj_ARI[15],
                          "err_coclust" = coclus_list$coclust_err,
                          "err_coclust_group1" = local_coclus_list[[1]][[1]], "err_coclust_group2" = local_coclus_list[[2]][[1]],
                          "err_coclust_group3" = local_coclus_list[[3]][[1]], "err_coclust_group4" = local_coclus_list[[4]][[1]],
                          "err_coclust_group5" = local_coclus_list[[5]][[1]], "err_coclust_group6" = local_coclus_list[[6]][[1]],
                          "err_coclust_group7" = local_coclus_list[[7]][[1]], "err_coclust_group8" = local_coclus_list[[8]][[1]],
                          "err_coclust_group9" = local_coclus_list[[9]][[1]], "err_coclust_group10" = local_coclus_list[[10]][[1]],
                          "err_coclust_group11" = local_coclus_list[[11]][[1]], "err_coclust_group12" = local_coclus_list[[12]][[1]],
                          "err_coclust_group13" = local_coclus_list[[13]][[1]], "err_coclust_group14" = local_coclus_list[[14]][[1]],
                          "err_coclust_group15" = local_coclus_list[[15]][[1]],
                          "err_L1_group1" = L1_list1[[1]][1], "err_L1_group2" = L1_list1[[1]][2],
                          "err_L1_group3" = L1_list1[[1]][3], "err_L1_group4" = L1_list1[[1]][4],
                          "err_L1_group5" = L1_list1[[1]][5], "err_L1_group6" = L1_list1[[1]][6],
                          "err_L1_group7" = L1_list1[[1]][7], "err_L1_group8" = L1_list1[[1]][8],
                          "err_L1_group9" = L1_list1[[1]][9], "err_L1_group10" = L1_list1[[1]][10],
                          "err_L1_group11" = L1_list1[[1]][11], "err_L1_group12" = L1_list1[[1]][12],
                          "err_L1_group13" = L1_list2[[1]][1], "err_L1_group14" = L1_list2[[1]][2],
                          "err_L1_group15" = L1_list2[[1]][3],
                          "err_L1_mean" = (d1 * L1_list1[[2]] + d2 * L1_list2[[2]])/d,
                          "err_wass2_group1" = wass2[[1]][1], "err_wass2_group2" = wass2[[1]][2],
                          "err_wass2_group3" = wass2[[1]][3], "err_wass2_group4" = wass2[[1]][4],
                          "err_wass2_group5" = wass2[[1]][5], "err_wass2_group6" = wass2[[1]][6],
                          "err_wass2_group7" = wass2[[1]][7], "err_wass2_group8" = wass2[[1]][8],
                          "err_wass2_group9" = wass2[[1]][9], "err_wass2_group10" = wass2[[1]][10],
                          "err_wass2_group11" = wass2[[1]][11], "err_wass2_group12" = wass2[[1]][12],
                          "err_wass2_group13" = wass2[[1]][13], "err_wass2_group14" = wass2[[1]][14],
                          "err_wass2_group15" = wass2[[1]][15],
                          "err_wass2_mean" = wass2[[2]]
  )

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
  Nburnin = burnin

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
    local_coclus_list[[j]] = Compute_coclust_error( real_partition[take_indices],sim_matrix[take_indices,take_indices] )
  }

  sd = sd1 # c'? un bug nella funzione, ho scritto sd invece che sigma. da sistemare!!!!
  L1_list1 = Compute_L1_dist(Pred = Pred_median[1:d1],
                             p_mix = mix_probs1, mu = mu1, sigma = sd1,
                             grid = grid)
  sd = sd2
  L1_list2 = Compute_L1_dist(Pred = Pred_median[(1+d1):d],
                             p_mix = mix_probs2, mu = mu2, sigma = sd2,
                             grid = grid)
  # wasserstein 2
  wass2 = vector("list",2)
  wass2[[1]] = rep(0,d)
  for(j in 1:d1){
    f_true = GDFMM::dmix(x = grid, w_j = mix_probs1[j,], mu_vec = mu1, sigma_vec = sd1)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f_true, p = 2)
  }
  for(j in (d1+1):d){
    f_true = GDFMM::dmix(x = grid, w_j = mix_probs2[j-d1,], mu_vec = mu2, sigma_vec = sd2)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f_true, p = 2)
  }
  wass2[[2]] = mean(wass2[[1]])
  names(wass2) = c("wass2_per_level","wass2_average")

  save_res_temp  = tibble("type" = "HDP", "time" = time,
                          "ARI_est_part" = ARI_VI,
                          "ARI_group1" = ARI_j[1], "ARI_group2" = ARI_j[2],
                          "ARI_group3" = ARI_j[3], "ARI_group4" = ARI_j[4],
                          "ARI_group5" = ARI_j[5], "ARI_group6" = ARI_j[6],
                          "ARI_group7" = ARI_j[7], "ARI_group8" = ARI_j[8],
                          "ARI_group9" = ARI_j[9], "ARI_group10" = ARI_j[10],
                          "ARI_group11" = ARI_j[11], "ARI_group12" = ARI_j[12],
                          "ARI_group13" = ARI_j[13], "ARI_group14" = ARI_j[14],
                          "ARI_group15" = ARI_j[15],
                          "K_mode" = K_mod,
                          "K_mode_group1" = Kj_mod[1], "K_mode_group2" = Kj_mod[2],
                          "K_mode_group3" = Kj_mod[3], "K_mode_group4" = Kj_mod[4],
                          "K_mode_group5" = Kj_mod[5], "K_mode_group6" = Kj_mod[6],
                          "K_mode_group7" = Kj_mod[7], "K_mode_group8" = Kj_mod[8],
                          "K_mode_group9" = Kj_mod[9], "K_mode_group10" = Kj_mod[10],
                          "K_mode_group11" = Kj_mod[11], "K_mode_group12" = Kj_mod[12],
                          "K_mode_group13" = Kj_mod[13], "K_mode_group14" = Kj_mod[14],
                          "K_mode_group15" = Kj_mod[15],
                          "K_ARI" = K_ARI,
                          "K_ARI_group1" = Kj_ARI[1], "K_ARI_group2" = Kj_ARI[2],
                          "K_ARI_group3" = Kj_ARI[3], "K_ARI_group4" = Kj_ARI[4],
                          "K_ARI_group5" = Kj_ARI[5], "K_ARI_group6" = Kj_ARI[6],
                          "K_ARI_group7" = Kj_ARI[7], "K_ARI_group8" = Kj_ARI[8],
                          "K_ARI_group9" = Kj_ARI[9], "K_ARI_group10" = Kj_ARI[10],
                          "K_ARI_group11" = Kj_ARI[11], "K_ARI_group12" = Kj_ARI[12],
                          "K_ARI_group13" = Kj_ARI[13], "K_ARI_group14" = Kj_ARI[14],
                          "K_ARI_group15" = Kj_ARI[15],
                          "err_coclust" = coclus_list$coclust_err,
                          "err_coclust_group1" = local_coclus_list[[1]][[1]], "err_coclust_group2" = local_coclus_list[[2]][[1]],
                          "err_coclust_group3" = local_coclus_list[[3]][[1]], "err_coclust_group4" = local_coclus_list[[4]][[1]],
                          "err_coclust_group5" = local_coclus_list[[5]][[1]], "err_coclust_group6" = local_coclus_list[[6]][[1]],
                          "err_coclust_group7" = local_coclus_list[[7]][[1]], "err_coclust_group8" = local_coclus_list[[8]][[1]],
                          "err_coclust_group9" = local_coclus_list[[9]][[1]], "err_coclust_group10" = local_coclus_list[[10]][[1]],
                          "err_coclust_group11" = local_coclus_list[[11]][[1]], "err_coclust_group12" = local_coclus_list[[12]][[1]],
                          "err_coclust_group13" = local_coclus_list[[13]][[1]], "err_coclust_group14" = local_coclus_list[[14]][[1]],
                          "err_coclust_group15" = local_coclus_list[[15]][[1]],
                          "err_L1_group1" = L1_list1[[1]][1], "err_L1_group2" = L1_list1[[1]][2],
                          "err_L1_group3" = L1_list1[[1]][3], "err_L1_group4" = L1_list1[[1]][4],
                          "err_L1_group5" = L1_list1[[1]][5], "err_L1_group6" = L1_list1[[1]][6],
                          "err_L1_group7" = L1_list1[[1]][7], "err_L1_group8" = L1_list1[[1]][8],
                          "err_L1_group9" = L1_list1[[1]][9], "err_L1_group10" = L1_list1[[1]][10],
                          "err_L1_group11" = L1_list1[[1]][11], "err_L1_group12" = L1_list1[[1]][12],
                          "err_L1_group13" = L1_list2[[1]][1], "err_L1_group14" = L1_list2[[1]][2],
                          "err_L1_group15" = L1_list2[[1]][3],
                          "err_L1_mean" = (d1 * L1_list1[[2]] + d2 * L1_list2[[2]])/d,
                          "err_wass2_group1" = wass2[[1]][1], "err_wass2_group2" = wass2[[1]][2],
                          "err_wass2_group3" = wass2[[1]][3], "err_wass2_group4" = wass2[[1]][4],
                          "err_wass2_group5" = wass2[[1]][5], "err_wass2_group6" = wass2[[1]][6],
                          "err_wass2_group7" = wass2[[1]][7], "err_wass2_group8" = wass2[[1]][8],
                          "err_wass2_group9" = wass2[[1]][9], "err_wass2_group10" = wass2[[1]][10],
                          "err_wass2_group11" = wass2[[1]][11], "err_wass2_group12" = wass2[[1]][12],
                          "err_wass2_group13" = wass2[[1]][13], "err_wass2_group14" = wass2[[1]][14],
                          "err_wass2_group15" = wass2[[1]][15],
                          "err_wass2_mean" = wass2[[2]]
  )

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
    local_coclus_list[[j]] = Compute_coclust_error( real_partition[take_indices],sim_matrix[take_indices,take_indices] )
  }

  sd = sd1 # c'? un bug nella funzione, ho scritto sd invece che sigma. da sistemare!!!!
  L1_list1 = Compute_L1_dist(Pred = Pred_median[1:d1],
                             p_mix = mix_probs1, mu = mu1, sigma = sd1,
                             grid = grid)
  sd = sd2
  L1_list2 = Compute_L1_dist(Pred = Pred_median[(1+d1):d],
                             p_mix = mix_probs2, mu = mu2, sigma = sd2,
                             grid = grid)
  # wasserstein 2
  wass2 = vector("list",2)
  wass2[[1]] = rep(0,d)
  for(j in 1:d1){
    f_true = GDFMM::dmix(x = grid, w_j = mix_probs1[j,], mu_vec = mu1, sigma_vec = sd1)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f_true, p = 2)
  }
  for(j in (d1+1):d){
    f_true = GDFMM::dmix(x = grid, w_j = mix_probs2[j-d1,], mu_vec = mu2, sigma_vec = sd2)
    wass2[[1]][j] = transport::wasserstein1d(a = Pred_median[[j]], b = f_true, p = 2)
  }
  wass2[[2]] = mean(wass2[[1]])
  names(wass2) = c("wass2_per_level","wass2_average")

  save_res_temp  = tibble("type" = "VecFDP-Cond", "time" = time,
                          "ARI_est_part" = ARI_VI,
                          "ARI_group1" = ARI_j[1], "ARI_group2" = ARI_j[2],
                          "ARI_group3" = ARI_j[3], "ARI_group4" = ARI_j[4],
                          "ARI_group5" = ARI_j[5], "ARI_group6" = ARI_j[6],
                          "ARI_group7" = ARI_j[7], "ARI_group8" = ARI_j[8],
                          "ARI_group9" = ARI_j[9], "ARI_group10" = ARI_j[10],
                          "ARI_group11" = ARI_j[11], "ARI_group12" = ARI_j[12],
                          "ARI_group13" = ARI_j[13], "ARI_group14" = ARI_j[14],
                          "ARI_group15" = ARI_j[15],
                          "K_mode" = K_mod,
                          "K_mode_group1" = Kj_mod[1], "K_mode_group2" = Kj_mod[2],
                          "K_mode_group3" = Kj_mod[3], "K_mode_group4" = Kj_mod[4],
                          "K_mode_group5" = Kj_mod[5], "K_mode_group6" = Kj_mod[6],
                          "K_mode_group7" = Kj_mod[7], "K_mode_group8" = Kj_mod[8],
                          "K_mode_group9" = Kj_mod[9], "K_mode_group10" = Kj_mod[10],
                          "K_mode_group11" = Kj_mod[11], "K_mode_group12" = Kj_mod[12],
                          "K_mode_group13" = Kj_mod[13], "K_mode_group14" = Kj_mod[14],
                          "K_mode_group15" = Kj_mod[15],
                          "K_ARI" = K_ARI,
                          "K_ARI_group1" = Kj_ARI[1], "K_ARI_group2" = Kj_ARI[2],
                          "K_ARI_group3" = Kj_ARI[3], "K_ARI_group4" = Kj_ARI[4],
                          "K_ARI_group5" = Kj_ARI[5], "K_ARI_group6" = Kj_ARI[6],
                          "K_ARI_group7" = Kj_ARI[7], "K_ARI_group8" = Kj_ARI[8],
                          "K_ARI_group9" = Kj_ARI[9], "K_ARI_group10" = Kj_ARI[10],
                          "K_ARI_group11" = Kj_ARI[11], "K_ARI_group12" = Kj_ARI[12],
                          "K_ARI_group13" = Kj_ARI[13], "K_ARI_group14" = Kj_ARI[14],
                          "K_ARI_group15" = Kj_ARI[15],
                          "err_coclust" = coclus_list$coclust_err,
                          "err_coclust_group1" = local_coclus_list[[1]][[1]], "err_coclust_group2" = local_coclus_list[[2]][[1]],
                          "err_coclust_group3" = local_coclus_list[[3]][[1]], "err_coclust_group4" = local_coclus_list[[4]][[1]],
                          "err_coclust_group5" = local_coclus_list[[5]][[1]], "err_coclust_group6" = local_coclus_list[[6]][[1]],
                          "err_coclust_group7" = local_coclus_list[[7]][[1]], "err_coclust_group8" = local_coclus_list[[8]][[1]],
                          "err_coclust_group9" = local_coclus_list[[9]][[1]], "err_coclust_group10" = local_coclus_list[[10]][[1]],
                          "err_coclust_group11" = local_coclus_list[[11]][[1]], "err_coclust_group12" = local_coclus_list[[12]][[1]],
                          "err_coclust_group13" = local_coclus_list[[13]][[1]], "err_coclust_group14" = local_coclus_list[[14]][[1]],
                          "err_coclust_group15" = local_coclus_list[[15]][[1]],
                          "err_L1_group1" = L1_list1[[1]][1], "err_L1_group2" = L1_list1[[1]][2],
                          "err_L1_group3" = L1_list1[[1]][3], "err_L1_group4" = L1_list1[[1]][4],
                          "err_L1_group5" = L1_list1[[1]][5], "err_L1_group6" = L1_list1[[1]][6],
                          "err_L1_group7" = L1_list1[[1]][7], "err_L1_group8" = L1_list1[[1]][8],
                          "err_L1_group9" = L1_list1[[1]][9], "err_L1_group10" = L1_list1[[1]][10],
                          "err_L1_group11" = L1_list1[[1]][11], "err_L1_group12" = L1_list1[[1]][12],
                          "err_L1_group13" = L1_list2[[1]][1], "err_L1_group14" = L1_list2[[1]][2],
                          "err_L1_group15" = L1_list2[[1]][3],
                          "err_L1_mean" = (d1 * L1_list1[[2]] + d2 * L1_list2[[2]])/d,
                          "err_wass2_group1" = wass2[[1]][1], "err_wass2_group2" = wass2[[1]][2],
                          "err_wass2_group3" = wass2[[1]][3], "err_wass2_group4" = wass2[[1]][4],
                          "err_wass2_group5" = wass2[[1]][5], "err_wass2_group6" = wass2[[1]][6],
                          "err_wass2_group7" = wass2[[1]][7], "err_wass2_group8" = wass2[[1]][8],
                          "err_wass2_group9" = wass2[[1]][9], "err_wass2_group10" = wass2[[1]][10],
                          "err_wass2_group11" = wass2[[1]][11], "err_wass2_group12" = wass2[[1]][12],
                          "err_wass2_group13" = wass2[[1]][13], "err_wass2_group14" = wass2[[1]][14],
                          "err_wass2_group15" = wass2[[1]][15],
                          "err_wass2_mean" = wass2[[2]]
  )

  save_res = save_res %>% rbind(save_res_temp)
  rm(GDFMM, sim_matrix, Kj, Pred_all,Local_Clustering)
}


save_res$type = as.factor(save_res$type)
#save(save_data, file = "data_exp3_cond_d15_variable.Rdat")
#save(save_res,  file = "res_exp3_cond_d15_variable.Rdat")
save(save_data, file = "data_exp4_d15_newmarginal.Rdat")
save(save_res,  file = "res_exp4_d15_newmarginal.Rdat")







