library(GDFMM)
# data generation ---------------------------------------------------------

d = 3 #number of groups
K = 3  # number of global clusters
mu = c(-20,0,20)
sd = c(1,1,1)
n_j = rep(200, d)
p = matrix(0, nrow = d, ncol = K)

set.seed(1234)
Kgruppo = c()
componenti_gruppo = NULL
data = matrix(NA, nrow = d, ncol = max(n_j))
cluster = matrix(NA, nrow = d, ncol = max(n_j))
for(j in 1:d){
  Kgruppo[j] = sample(1:3,1)
  componenti_gruppo[[j]] = sample(1:3,Kgruppo[j], replace = F)
  p[j,1:Kgruppo[j]] = rep(1/Kgruppo[j], Kgruppo[j])
  appoggio = genera_mix_gas(n = n_j[j], pro = p[j,1:Kgruppo[j]], means = mu[ componenti_gruppo[[j]] ],
                            sds = sd[ componenti_gruppo[[j]] ] )

  data[j, 1:n_j[j]] = appoggio$y
  cluster[j, 1:n_j[j]] = appoggio$clu
}


#x11();hist(data[1,])
# Run  --------------------------------------------------------------------


niter  <- 2000
burnin <- 1000
thin   <- 1

option<-list("Mstar0" = 5,"Lambda0" = 5,"mu0" = 0,"nu0"=10,"sigma0"= 1, "gamma0" = 1,
             "Adapt_MH_hyp1"= 0.7,"Adapt_MH_hyp2"= 0.234, "Adapt_MH_power_lim"=10, "Adapt_MH_var0"=1,
             "k0"= 1/10, "alpha_gamma"=1,
             "beta_gamma"=1, "alpha_lambda"=1, "beta_lambda"=1,
             "UpdateU" = T, "UpdateM" = T, "UpdateGamma" = T, "UpdateS" = T,
             "UpdateTau" = T, "UpdateLambda" = T
)


GDFMM = GDFMM_sampler(data, niter, burnin, thin, seed = 123, option = option)




# Analisi output ----------------------------------------------------------


#K
summary(GDFMM$K)
x11();plot(GDFMM$K, type = 'l')

#Mstar
summary(GDFMM$Mstar)
x11();plot(GDFMM$Mstar, type = 'l')




# Predictive --------------------------------------------------------------

l_grid = 10000
grid = seq(-25,25,length.out = l_grid)

Pred = pred_uninorm(idx_group = 1, grid = grid, fit = GDFMM)
Pred = predictive(idx_group = 1, grid = grid, fit = GDFMM)
x11();matplot(x = grid, y = t(Pred), type = 'l', col = 'red')

Pred_all = predictive_all_groups(grid = grid, fit = GDFMM)

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



