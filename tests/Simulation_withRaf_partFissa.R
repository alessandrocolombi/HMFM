library(GDFMM)
library(ACutils)
# data generation ---------------------------------------------------------

d = 10               # number of groups
K = 3               # number of global clusters
mu = c(-20,0,20)   # vectors of means
sd = c(1,1,1)      # vector of sd
n_j = rep(200, d)  # set cardinality of the groups
p = matrix(0, nrow = d, ncol = K) # matrix with components weights

set.seed(124123)
Kgruppo = c()
componenti_gruppo = NULL
data = matrix(NA, nrow = d, ncol = max(n_j))     # d x max(n_j) matrix
cluster = matrix(NA, nrow = d, ncol = max(n_j))  # d x max(n_j) matrix
real_partition = c()      # real_partition is a vector of length sum(n_j), it collects all the group membership.
                          # values are collected level by level, so first all the values in level 1, the all values in level 2 and so on
                          # cluster label must always start from 0!
for(j in 1:d){
  Kgruppo[j] = sample(1:K,1) # number of clusters in each level
  componenti_gruppo[[j]] = sample(1:K,Kgruppo[j], replace = F) # choose the components
  p[j,1:Kgruppo[j]] = rep(1/Kgruppo[j], Kgruppo[j]) # set the weights all equals
  appoggio = genera_mix_gas(n = n_j[j], pro = p[j,1:Kgruppo[j]], means = mu[ componenti_gruppo[[j]] ],
                            sds = sd[ componenti_gruppo[[j]] ] )

  data[j, 1:n_j[j]] = appoggio$y
  #cluster[j, 1:n_j[j]] = appoggio$clu, #errore, genera_mix_gas usa sempre indici che partono da 1!
  cluster[j, 1:n_j[j]] = unlist(lapply(1:n_j[j], function(h){componenti_gruppo[[j]][appoggio$clu[h]]}))
  real_partition = c(real_partition, cluster[j, 1:n_j[j]])
}
# In real partition devo avere valore da 0 a K senza buchi. Per esempio, se ho solo 0 e 2 non va bene!
# quella che viene modificata dentro il sampler è
# partion_within_sampler = arrange_partition(real_partition)
N_m = table(real_partition)

#x11();hist(data[1,])
#x11();hist(data[2,])
#x11();hist(data[3,])


data_level1 = data[cluster==1]
data_level2 = data[cluster==2]
data_level3 = data[cluster==3]

#mean1 = mean(data_level1); var1 = var(data_level1); x11(); hist(data_level1); N_m1 = length(data_level1)
#mean2 = mean(data_level2); var2 = var(data_level2); x11(); hist(data_level2); N_m2 = length(data_level2)
#mean3 = mean(data_level3); var3 = var(data_level3); x11(); hist(data_level3); N_m3 = length(data_level3)
#c(N_m1, N_m2, N_m3)
#c(mean1, mean2, mean3)
#c(var1, var2, var3)

# Run  --------------------------------------------------------------------


niter  <- 1000
burnin <- 1000
thin   <- 1

option<-list("nu" = 0.7, "Mstar0" = 3, "Lambda0" = 3, "mu0" = 0,"sigma0"= 1, "gamma0" = 10,
             "Adapt_MH_hyp1"= 0.7,"Adapt_MH_hyp2"= 0.234, "Adapt_MH_power_lim"=10, "Adapt_MH_var0"=1,
             "k0"= 1/10, "nu0"=10, "alpha_gamma"=1,
             "beta_gamma"=1, "alpha_lambda"=1, "beta_lambda"=1,
             "UpdateU" = T, "UpdateM" = T, "UpdateGamma" = T, "UpdateS" = T,
             "UpdateTau" = T, "UpdateLambda" = T, "partition" = real_partition
)

#GDFMM = GDFMM_sampler(data, niter, burnin, thin, seed = 123, option = option)
GDFMM = GDFMM_sampler(data, niter, burnin, thin, seed = 123, FixPartition = T, option = option)

View(GDFMM)

# - Mstar0 è il numero iniziale di componenti non allocate, decidono loro a priori che K0 = 1, quindi si parte
#   da M = Mstar0 + 1. Bisogna fare in modo che l'utente possa decidere il numero iniziale di componenti allocate
#   e partire con Mstar0 = 0
# - Non va bene non entrare nemmeno nella funzione FC_M::update, perché quando si aggiorna la partizione viene
#   modificato K ma non viene aggiornato Mstar e nemmeno M. Vengono aggiornate in FC_M::update dopo aver estratto
#   il valore nuovo. Quindi se non entro lì, niente ha senso dopo l'update della partizione-
#   --> ho messo un update manuale
# - Ho dovuto fare il nuovo allocate di tau e S. è strano perché distrugge il contenuto di quello che c'era prima,
#   però ai fini dell'update non sembra essere un problema. è un problema invece se voglio tenere S e gamma fissi,
#   a quel punto cosi mi dà errore.
# - C'era un errore nella varianza di tau, dividevano per N_k[m] e non per N_k[m]-1. Cambia poco

# - Assicurati che possa andare fissando S e gamma. anche tau forse viene inizializzato random
# Analisi output ----------------------------------------------------------


#K
summary(GDFMM$K)
x11();plot(GDFMM$K, type = 'l', main = "K")

#Mstar
summary(GDFMM$Mstar)
x11();plot(GDFMM$Mstar, type = 'l', main = "Mstar")


#log_sum = sum(  gamma_j * log( (u_j + 1) )   )
summary(GDFMM$log_sum)
x11();plot(GDFMM$log_sum, type = 'l', main = "log_sum")

#Parametro Poisson: lambda * exp(-log_sum)
summary(GDFMM$lambda*exp(-GDFMM$log_sum))
x11();plot(GDFMM$lambda*exp(-GDFMM$log_sum), type = 'l', main = "Poisson parameter")

#gammas
gamma = GDFMM$gamma
post_mean_gamma = rowMeans(gamma)
x11();matplot(x = 1:ncol(gamma), y = t(gamma), type = 'l')
# Predictive --------------------------------------------------------------

l_grid = 200
grid = seq(-25,25,length.out = l_grid)

# Predictive in all groups
Pred_all = predictive_all_groups(grid = grid, fit = GDFMM)

for(j in 1:d){
  x11();hist(data[j,], freq = F, breaks = l_grid/10, col = ACutils::t_col("darkred", 70), xlim = range(grid),
             main = paste0("level = ",j))
  matplot(x = grid, y = t(Pred_all[[j]]), type = 'l', col = 'black', lty = 1, lwd = 2, add = T)
  points(x = data[j,], y = rep(0, length(data[j,])), pch = 16, col = ACutils::t_col("darkred", 10))

}

#Pred = pred_uninorm(idx_group = 1, grid = grid, fit = GDFMM)
Pred = predictive(idx_group = 1, grid = grid, fit = GDFMM)
x11();hist(data[1,], freq = F, breaks = l_grid/10, col = ACutils::t_col("darkred", 70), xlim = range(grid) )
matplot(x = grid, y = t(Pred), type = 'l', col = 'black', lty = 1, lwd = 2, add = T)
points(x = data[1,], y = rep(0, length(data[1,])), pch = 16, col = ACutils::t_col("darkred", 10))








x11();par(mfrow = c(d,1))
lapply(Pred_all, function(pred){
  matplot(t(pred), type = 'l', col = 'red', main = 'Predictive')
})

# questo plot è sbagliatissimo. trova praticamente una sola componenti
# - guarda i valori che escono in (mu,sigma):
# --> Ci sono tante mu vicino a -10 ma quella componente praticamente non si vede.
## --> Le sigma sono orrende, c'è sempre un valore enorme e uno vicino a 0.



# Test predictive ---------------------------------------------------------
# Test solo per vedere se le funzioni pred_uniform e predictive sono corrette
niter  <- 10
GDFMM$Mstar = rep(0,niter)
GDFMM$K     = rep(3,niter)
for(i in 1:niter){
  GDFMM$mu[[i]] = c(-20,0,20) + rnorm(n=1,sd=0.01)
  GDFMM$sigma[[i]] = c(1,1,1)
  GDFMM$S[[i]] = matrix(1/d,nrow = d, ncol = GDFMM$K[i] + GDFMM$Mstar[i], byrow = T)
}


l_grid = 100
grid = seq(-25,25,length.out = l_grid)

Pred = pred_uninorm(idx_group = 1, grid = grid, fit = GDFMM)
Pred = predictive(idx_group = 1, grid = grid, fit = GDFMM)
x11();matplot(x = grid, y = t(Pred), type = 'l', col = 'red')


