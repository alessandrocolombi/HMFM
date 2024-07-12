# Libraries
suppressWarnings(suppressPackageStartupMessages(library(GDFMM)))
suppressWarnings(suppressPackageStartupMessages(library(hdp)))

setwd("./SimulationStudy")
Rcpp::sourceCpp("../src/lastirling1.cpp")

# Function ----------------------------------------------------------------

SimStudy_runningtimes = function(seed, ndatas){
  suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
  suppressWarnings(suppressPackageStartupMessages(library(GDFMM)))
  suppressWarnings(suppressPackageStartupMessages(library(hdp)))
  # Initialize return object
  names_results = c("HMFM_marg","HDP","HMFM_cond")
  results = lapply(1:length(names_results), FUN = function(x){rep(-1,length(ndatas))})
  names(results) = names_results

  # Set number of iterations
  niter  <-  5000
  burnin <-  5000
  thin   <- 1

  # Main for-loop for all values of n
  for(ii in 1:length(ndatas)){

    cat("\n ----------------------------- \n")
    cat(" n =  ",ndatas[ii])
    cat("\n ----------------------------- \n")

    #1) Data Generation
    d  = 2                   # number of groups
    K  = 3                   # number of global clusters
    mu = c(-5,0,7)         # vectors of means
    sd = c(sqrt(0.5),sqrt(0.5),sqrt(0.5))       # vector of sd
    n_j=rep(ndatas[ii]/2, d)           # set cardinality of the groups
    mix_probs = matrix(c(1/3,1/3,1/3,1/3,1/3,1/3), nrow = d, ncol = K, byrow = T)


    set.seed(seed) #set the seed
    genD = simulate_data(d = d, K = K, p = mix_probs, mu = mu, sd = sd, n_j = n_j, seed = seed)
    data = genD$data
    real_partition = genD$real_partition

    data_list = vector("list",length = d)
    for(j in 1:d)
      data_list[[j]] = data[j,1:n_j[j]]


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
    Ncenters = 4
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
    results$HMFM_marg[ii] = as.numeric(time)
    rm(GDFMM)


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

    results$HDP[ii] = time
    rm(fit)


    # c) HMFM cond ------------------------------------------------------------

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
    Exp_Lambda   =  5#25
    Var_Lambda   =  5
    gamma_guess  =  0.5#1/sum(n_j)
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

    Ncenters = 4
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
    gamma0 = rep(0.25,d)#rep(0.025,d)#1/n_j # per ora il valore magico è 0.025
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
      GDFMM = ConditionalSampler(dt[1:11], niter, burnin, thin, seed = 123, option = option,
                                 FixPartition = F,
                                 P0.prior = prior)
    time = tictoc::toc()
    time = time$toc - time$tic

    results$HMFM_cond[ii] = time
    rm(GDFMM)
  }

  # Normalize wrt number of iterations
  results$HMFM_marg = results$HMFM_marg/(niter + burnin)
  results$HDP = results$HDP/(niter + burnin)
  results$HMFM_cond = results$HMFM_cond/(niter + burnin)

  # return
  return(results)
}


# Run - no precomputation ---------------------------------------------------------------------
ndatas = c(24,50,100,200,400,800,1600)
Nrep  = 10

seed0 = 1605
set.seed(seed0)
seeds = sample(1:999999, size = Nrep)
num_cores = 7

tictoc::tic()
  cluster <- parallel::makeCluster(num_cores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)
  parallel::clusterExport(cluster, list())
  res = parallel::parLapply( cl = cluster, seeds,
                             fun = SimStudy_runningtimes,
                             ndatas = ndatas)
  parallel::stopCluster(cluster)
tictoc::toc()


# Stirling numbers computation -----------------------------------------------

B = rep(1,length(ndatas)) # repetitions of Striling computations
B[1] <- B[2] <- 10000
B[3] <- B[4] <- 1000
B[5] <- B[6] <- 10
res_Stirling = lapply(1:Nrep,function(x){rep(-1,length(ndatas))})


for(rep in 1:Nrep){
  for(ii in 1:length(ndatas)){

    cat("\n ----------------------------- \n")
    cat(" repetition =  ",rep,"; ",ii)
    cat("\n ----------------------------- \n")

    tictoc::tic()
      lstr = lapply(1:B[ii], function(x){GDFMM:::lastirling1(ndatas[ii])})
    time = tictoc::toc()
    time = time$toc - time$tic
    res_Stirling[[rep]][ii] = time/B[ii]
  }
}


# Plot --------------------------------------------------------------------
col_type = c("chartreuse3","orange","darkred","cadetblue4")
add = function(x){Reduce("+", x)} # used to sum vectors/matrices stored in elements of a list

HDP_line = add( lapply(res, function(x){log(x$HDP)}) )/Nrep
HMFMmarg = add( lapply(res, function(x){log(x$HMFM_marg)}) )/Nrep
HMFMcond = add( lapply(res, function(x){log(x$HMFM_cond)}) )/Nrep
Stirling_line = add( lapply(res_Stirling, function(x){log(x)}) )/Nrep

plot_lines = list(HDP_line,HMFMmarg,HMFMcond,Stirling_line)

grid = log(ndatas)
grid_lines = seq(5.5,7.3,length.out = 4)
par( mar = c(4,4,1,1), bty = "l")
plot(0,0,main = " ", ylab = " log( Execution Times )", xlab = "log(n)" ,
     type = "n",
     xlim = log(c(20,1800)), ylim = c(-10.5,-2))
grid(lty = 1,lwd = 1, col = "gray90" )
for(m in 1:4){
  lines(x = grid, y = plot_lines[[m]], col = col_type[m], lwd = 3)
  lines(x = grid_lines, y = grid_lines - 14.5,
        col = "black", lty = 2, lwd = 2, type = "b", pch = 16)
  lines(x = grid_lines, y = 2*grid_lines - 17.4,
        col = "grey45", lty = 3, lwd = 2, type = "b", pch = 4)
  legend("topleft",
         legend = c("HDP","HMFM-marg","HMFM-cond","Stirling number"),
         col = col_type, lwd = 2, lty = 1 )
}

# Warning: the final plot as well as the values on the y-axis depends on the machine this experiment is run
