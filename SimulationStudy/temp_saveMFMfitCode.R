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
xrange = c(-7,7)
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
