library(GDFMM)
library(ACutils)
library(tidyverse)
library(RColorBrewer)

# color palette -----------------------------------------------------------
mycol = hcl.colors(n=3,palette = "Zissou1")

# data generation ---------------------------------------------------------

d = 3               # number of groups
K = 4               # number of global clusters
mu = c(-20,-10,0, 10)   # vectors of means
sd = c(1,1,1,1)      # vector of sd
n_j = rep(200, d)  # set cardinality of the groups
seed = 20051131

genD = generate_data(d=d, K=K, mu = mu, sd = sd, n_j = n_j, seed = seed)
data = genD$data
real_partition = genD$real_partition


mycol_cluster = brewer.pal(n=K, name = "Dark2")
x11();hist(data[1,], breaks = 100, freq = F)
x11();hist(data[2,])
x11();hist(data[3,])
# Run  --------------------------------------------------------------------


niter  <- 1000
burnin <- 3000
thin   <- 1

option<-list("nu" = 1, "Mstar0" = 10, "Lambda0" = 3, "mu0" = 0,"sigma0"= 1, "gamma0" = 1,
             "Adapt_MH_hyp1"= 0.7,"Adapt_MH_hyp2"= 0.234, "Adapt_MH_power_lim"=10, "Adapt_MH_var0"=1,
             "k0"= 1/10, "nu0"=10, "alpha_gamma"=1,
             "beta_gamma"=1, "alpha_lambda"=1, "beta_lambda"=1,
             "UpdateU" = T, "UpdateM" = F, "UpdateGamma" = T, "UpdateS" = T,
             "UpdateTau" = T, "UpdateLambda" = T
)

#GDFMM = GDFMM_sampler(data, niter, burnin, thin, seed = 123, option = option)
GDFMM = GDFMM_sampler(data, niter, burnin, thin, seed = 123, FixPartition = F, option = option)

View(GDFMM)

# Analisi output ----------------------------------------------------------

#K
summary(GDFMM$K)
x11();plot(GDFMM$K, type = 'l', main = "K")
x11();hist(GDFMM$K)

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
x11();print(
  tibble(value = data[j,]) %>%
    ggplot(aes(x=value)) +
    geom_histogram(aes(y = ..density..), color = 'black', fill = 'darkred', alpha = 0.3,
                   binwidth = 0.5) +
    geom_point(y = rep(0,length(data[j,])), col = 'darkred') +
    geom_path(data = as_tibble(t(Pred_all[[j]])), color = 'black', aes(x=grid,y=`50%`), size = 1.2 ) +
    geom_ribbon(data = as_tibble(t(Pred_all[[j]])), fill = mycol[1],
                aes(x = grid, ymin=`2.5%`, ymax=`97.5%`, y=`50%`, fill = "band"),
                alpha = 0.5) +
    labs(title=paste0("level = ",j), x = "x-axis", y = "y-axis") + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))  # set title in center and set the
)
}


# Nel plot sotto colori i pallini secondo il vero cluster.
# Assurdo ma gli istogrammi non funzionano più, evidentemente raggruppa in modo diverso
idx_start = 1
idx_end = n_j[1]
for(j in 1:d){

  x11();print(
  tibble(value = data[j,],
         true_clus = as.factor(real_partition[idx_start:idx_end])
         ) %>%
    ggplot(aes(x=value, col = true_clus, fill = true_clus)) +
    geom_point(y = rep(0.005,length(data[j,])), size = 2) +
    geom_histogram(data = tibble(value = data[j,]),
                   aes(x=value, y = ..density..),
                   color = 'black', alpha = 0.3,
                   binwidth = 0.5, inherit.aes = F) +
    scale_color_manual(values = mycol_cluster[ sort(unique(real_partition[idx_start:idx_end])) ] ) +
    scale_fill_manual(values = mycol_cluster[ sort(unique(real_partition[idx_start:idx_end])) ] ) +
    geom_path(data = as_tibble(t(Pred_all[[j]])), color = 'black', aes(x=grid,y=`50%`), size = 1.2,
              inherit.aes = F) +
    labs(title=paste0("level = ",j), x = "x-axis", y = "y-axis") + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) + # set title in center and set the
    geom_ribbon(data = as_tibble(t(Pred_all[[j]])), fill = mycol[1],
                aes(x = grid, ymin=`2.5%`, ymax=`97.5%`, y=`50%`, fill = "band"),
                alpha = 0.5, inherit.aes = F)
  )

  idx_start = idx_end + 1
  idx_end = idx_end + n_j[j]
}


#Pred = pred_uninorm(idx_group = 1, grid = grid, fit = GDFMM)
Pred = predictive(idx_group = 1, grid = grid, fit = GDFMM)

idx_group = 1
x11();tibble(value = data[idx_group,]) %>%
    ggplot(aes(x=value)) +
    geom_histogram(aes(y = ..density..), color = 'black', fill = 'darkred', alpha = 0.3,
                   binwidth = 0.5) +
    geom_point(y = rep(0,length(data[idx_group,])), col = 'darkred') +
    geom_path(data = as_tibble(t(Pred)), color = 'black', aes(x=grid,y=`50%`), size = 1.2 ) +
    geom_ribbon(data = as_tibble(t(Pred)), fill = mycol[1],
                aes(x = grid, ymin=`2.5%`, ymax=`97.5%`, y=`50%`, fill = "band"),
                alpha = 0.5) +
    labs(title=paste0("level = ",j), x = "x-axis", y = "y-axis") + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))  # set title in center and set the


tibble(value = data[idx_group,], true_clus = as.factor(real_partition[1:200]))

tibble(value = data[idx_group,], true_clus = as.factor(real_partition[1:200])) %>%
  ggplot(aes(x=value, col = true_clus)) +
  geom_point(y = rep(0,length(data[idx_group,]))) +
  geom_histogram(aes(y = ..density.., fill = true_clus),  alpha = 0.3,
                 binwidth = 0.5) +
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2")



mycol_cluster

  geom_histogram(aes(y = ..density..), color = 'black', fill = 'darkred', alpha = 0.3,
                 binwidth = 0.5) +
  geom_point(y = rep(0,length(data[idx_group,])), col = 'darkred') +
  geom_path(data = as_tibble(t(Pred)), color = 'black', aes(x=grid,y=`50%`), size = 1.2 ) +
  geom_ribbon(data = as_tibble(t(Pred)), fill = mycol[1],
              aes(x = grid, ymin=`2.5%`, ymax=`97.5%`, y=`50%`, fill = "band"),
              alpha = 0.5) +
  labs(title=paste0("level = ",j), x = "x-axis", y = "y-axis") + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))  # set title in center and set the


# Vecchii plot
x11();hist(data[1,], freq = F, breaks = l_grid/10, col = ACutils::t_col("darkred", 70), xlim = range(grid) )
matplot(x = grid, y = t(Pred), type = 'l', col = 'black', lty = 1, lwd = 2, add = T)
points(x = data[1,], y = rep(0, length(data[1,])), pch = 16, col = ACutils::t_col("darkred", 10))


for(j in 1:d){
  x11();hist(data[j,], freq = F, breaks = l_grid/10, col = ACutils::t_col("darkred", 70), xlim = range(grid),
             main = paste0("level = ",j))
  matplot(x = grid, y = t(Pred_all[[j]]), type = 'l', col = 'black', lty = 1, lwd = 2, add = T)
  points(x = data[j,], y = rep(0, length(data[j,])), pch = 16, col = ACutils::t_col("darkred", 10))

}






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


