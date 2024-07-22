

#  Librerie ---------------------------------------------------------------


suppressWarnings(suppressPackageStartupMessages(library(GDFMM)))
suppressWarnings(suppressPackageStartupMessages(library(ACutils)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(library(salso)))
suppressWarnings(suppressPackageStartupMessages(library(wesanderson)))
suppressWarnings(suppressPackageStartupMessages(library(mcclust.ext)))
suppressWarnings(suppressPackageStartupMessages(library(abind)))


# Custom functions --------------------------------------------------------


my_invgamma = function(x, nu0, sigma0){
  a = nu0/2
  b = a*sigma0
  b^a/gamma(a) * (1/x)^(a+1) * exp(-b/x)
}



# Load data ---------------------------------------------------------------

data_aligned = read.csv("../Shotput_longformat_preproc_nofew50.csv", row.names=1)
data_aligned = as_tibble(data_aligned) %>% mutate(ID = as.factor(ID),
                                                  SeasonNumber = as.factor(SeasonNumber),
                                                  Gender = as.factor(Gender),
                                                  Environment = as.factor(Environment)) %>%
  select(ID,SeasonNumber,Result,
         Gender,Environment,Age,AgeEntrance,
         Days,t_ji)


# select first 15 seasons only
n = nrow(data_aligned %>% distinct(ID))
selectIDs = 1:n
d = 15
IDs  = data_aligned %>% distinct(ID) %>% pull(ID)
data_longform_input = data_aligned %>%
  filter(ID %in% IDs[selectIDs]) %>%
  filter(SeasonNumber %in% as.factor(1:d)) %>%
  filter(ID != "76011") %>%
  select(ID,SeasonNumber,Result,Gender,Environment,Age,AgeEntrance,Days,t_ji)

# Center data and covariates
data_longform_input$Result      = data_longform_input$Result
data_longform_input$Result      = data_longform_input$Result      - mean(data_longform_input$Result)
data_longform_input$Age         = data_longform_input$Age         - mean(data_longform_input$Age)
data_longform_input$AgeEntrance = data_longform_input$AgeEntrance - mean(data_longform_input$AgeEntrance)


# Hyperparameters: P0 ---------------------------------------------------------

#Res_range = range( data_longform_input$Result )
Res_range = quantile(data_longform_input$Result, probs = c(0.005,0.995)) # provo nuova
R = Res_range[2] - Res_range[1]
mu0 = mean(data_longform_input$Result) # should be 0
# Per il momento scegliamo
# k0  = 10/R^2; 1/R^2
# nu0 = 10;     4
# sigma0 = 10;  10
k0  = 1/R^2
nu0 = 4
sigma0 = 20/2

scale = sqrt( (k0 + 1)/(k0) * sigma0 )
mean_marginal = mu0
var_marginal  = nu0/(nu0-2) * scale^2





# Initial values ----------------------------------------------------------

data_med4season = data_longform_input %>% group_by(ID,SeasonNumber) %>%
  mutate(MedResult = mean(Result)) %>%
  select(ID,SeasonNumber,Result, MedResult,
         Gender,Environment,Age,AgeEntrance,
         Days,t_ji) %>%
  distinct(ID,SeasonNumber, .keep_all = TRUE) %>%
  ungroup() %>% arrange(SeasonNumber)


# Con covariate, devo prima fare una regressione
res = lm(Result ~ Gender, data = as.data.frame(data_med4season))$residuals
Ncenters = 20
Kmeans0 = kmeans(x = res, centers = Ncenters, iter.max = 50, nstart = 10 )
KmeansCl = Kmeans0$cluster
centers = Kmeans0$centers
data_med4season = data_med4season %>% cbind("Kmeans" = KmeansCl)

# Tibble con tutti i dati iniziali

data_with_init = data_longform_input %>% left_join(data_med4season %>%
                                                     select("ID","SeasonNumber","MedResult","Kmeans"),
                                                   by = c("ID","SeasonNumber"))



# Setup -------------------------------------------------------------------

# Leggo i dati nel modo in cui poi devono essere inseriti

dt = input_handle(data_with_init[,c(1:3,11,4)], intercept = FALSE)

n = dt$n
d = dt$d
r = dt$r
n_j = dt$n_j
#n;d;r;sum(n_j);sum(dt$N_ji)

# Hyperparameters: process ------------------------------------------------

Exp_Lambda   = 25
Var_Lambda   =  3
gamma_guess  =  1/sum(n_j) #0.005
Lambda_guess = Exp_Lambda

b_lambda = Exp_Lambda/Var_Lambda
a_lambda = Exp_Lambda * b_lambda

a_gamma = a_lambda/d
b_gamma = a_gamma / (gamma_guess * Lambda_guess)



# Run ---------------------------------------------------------------------

niter  <-  20000 #number of saved iteration
burnin <-  50000
thin   <-     10

#total number of iterations is burnin + thin*niter


# initial values
beta0 = rep(-2,dt$r)
Sigma0 = diag(dt$r)

Lambda0 = 5
gamma0 = rep(0.01,d)
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
                      "IncludeCovariates" = TRUE,
                      "UpdateU" = T, "UpdateM" = T, "UpdateGamma" = T,
                      "UpdateS" = T, "UpdateTau" = T, "UpdateLambda" = T,
                      "UpdateBeta" = T )

prior = "Normal-InvGamma"

GDFMM = ConditionalSampler(dt[1:11], niter, burnin, thin, seed = 123, option = option, FixPartition = F,
                           P0.prior = prior, algorithm = "Neal2")


# Clustering --------------------------------------------------------------


# Get labels for each iterations for each data point
part_matrix <- GDFMM$Partition[(niter/2):niter,] #GDFMM$Partition is a (n_iter x n_data) matrix

# Compute similarity matrix
sim_matrix <- psm(part_matrix)

#heatmap(sim_matrix)
VI_sara = minVI(sim_matrix)


dt$finalPartition = vector("list", length = d)
dt$finalPartition = lapply(1:d, FUN = function(s){dt$finalPartition[[s]] = vector("list", length = n) })
dt$Clustering = vector("list", length = d)
dt$Clustering = lapply(1:d, FUN = function(s){dt$finalPartition[[s]] = vector("list", length = n) })

counter_obs = 1
for(j in 1:d){
  for(i in 1:n){
    if(dt$N_ji[j,i] > 0){
      dt$Clustering[[j]][[i]] = GDFMM$Partition[,counter_obs]+1 # parte da 0 quella salvata
      dt$finalPartition[[j]][[i]] = VI_sara$cl[counter_obs]
      counter_obs = counter_obs + 1
    }
  }
}





# Level dependent clustering ----------------------------------------------

Local_Clustering = list("list", length = d)
idx_start = c(1,cumsum(n_j)[1:(d-1)]+1)
idx_end = cumsum(n_j)
Kj = matrix(0,nrow = niter, ncol = d)
for(j in 1:d){
  Local_Clustering[[j]] = GDFMM$Partition[ , idx_start[j]:idx_end[j] ]
  Kj[,j] = apply(Local_Clustering[[j]], 1, FUN = function(Part_it){length(table(Part_it))})
}


# save --------------------------------------------------------------------

application_result = list("GDFMM" = GDFMM,
                          "sim_matrix" = sim_matrix,
                          "VI_sara" = VI_sara,
                          "dt" = dt,
                          "Local_Clustering" = Local_Clustering,
                          "Kj" = Kj)

save(application_result, file = "application_result.Rdat")


# heatmap -----------------------------------------------------------------

tictoc::tic()
heatmap(sim_matrix) #34 minuti ed ? veramente brutta
time = tictoc::toc()
time = time$toc - time$toc



