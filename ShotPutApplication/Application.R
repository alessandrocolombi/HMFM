#  Libraries ---------------------------------------------------------------

suppressWarnings(suppressPackageStartupMessages(library(GDFMM)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(salso)))
suppressWarnings(suppressPackageStartupMessages(library(mcclust.ext)))
suppressWarnings(suppressPackageStartupMessages(library(abind)))
suppressWarnings(suppressPackageStartupMessages(library(knitr)))
# suppressWarnings(suppressPackageStartupMessages(library(ACutils)))
# suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
# suppressWarnings(suppressPackageStartupMessages(library(wesanderson)))

setwd(here::here())
data("ShotPutData")
setwd("./ShotPutApplication")
seed0 = 271296
set.seed(seed0)
# Load data ---------------------------------------------------------------
data_longform_input = ShotPutData

# Center data and covariates
data_longform_input$Result      = data_longform_input$Result
data_longform_input$Result      = data_longform_input$Result      - mean(data_longform_input$Result)
data_longform_input$Age         = data_longform_input$Age         - mean(data_longform_input$Age)
data_longform_input$AgeEntrance = data_longform_input$AgeEntrance - mean(data_longform_input$AgeEntrance)


# Hyperparameters: P0 ---------------------------------------------------------

Res_range = range( data_longform_input$Result )
# Res_range = quantile(data_longform_input$Result, probs = c(0.005,0.995))
R = Res_range[2] - Res_range[1]
mu0 = mean(data_longform_input$Result) # should be 0
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


# covariates

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


dt = input_handle(data_with_init[,c(1:3,11,4)], intercept = FALSE)

n = dt$n
d = dt$d
r = dt$r
n_j = dt$n_j

# Hyperparameters: process ------------------------------------------------

Exp_Lambda   = 25
Var_Lambda   =  3
gamma_guess  =  1/sum(n_j)
Lambda_guess = Exp_Lambda

b_lambda = Exp_Lambda/Var_Lambda
a_lambda = Exp_Lambda * b_lambda

a_gamma = a_lambda/d
b_gamma = a_gamma / (gamma_guess * Lambda_guess)



# Run ---------------------------------------------------------------------

niter  <-   20000 #number of saved iteration
burnin <-   50000
thin   <-      10

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

beepr::beep()
# Clustering --------------------------------------------------------------


# Get labels for each iterations for each data point
part_matrix <- GDFMM$Partition[(niter/2):niter,] #GDFMM$Partition is a (n_iter x n_data) matrix

# Compute similarity matrix
sim_matrix <- psm(part_matrix)

VI_sara = minVI(sim_matrix)
table(VI_sara$cl)
beepr::beep()


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




# Analysis: betas -------------------------------------------------------------------

Beta = GDFMM$beta[15000:20000] #get sampled values
Beta_table = abind(Beta, along = 3) # transform into array  (d x r x niter)

Beta_post_mean = apply(Beta_table, c(1,2), mean) # posterior mean
Beta_post_quant = apply(Beta_table, c(1,2), quantile, prob=c(0.025,0.5,0.975)) # calcolo quantili


Beta_tibble = tibble("coefficient" = as.numeric(), "season" = as.integer())
for(j in 1:d){
  temp = tibble("coefficient" = Beta_table[j,1,], "season" = rep(j,length(Beta_table[j,1,]))  )
  Beta_tibble = rbind(Beta_tibble, temp)
}
Beta_tibble = Beta_tibble %>% mutate(season = as.factor(season))


col_season = "grey47"
name = "coefficient"

lower  = Beta_tibble %>% group_by(season) %>% summarise(lowerCL  = quantile(coefficient, probs = c(0.025)))
median = Beta_tibble %>% group_by(season) %>% summarise(medianCL = quantile(coefficient, probs = c(0.500)))
upper  = Beta_tibble %>% group_by(season) %>% summarise(upperCL  = quantile(coefficient, probs = c(0.975)))
Beta_CI = lower %>% left_join(median) %>% left_join(upper)

Beta_CI %>%
  ggplot(aes(y = medianCL, x = season, fill = season)) +
  geom_point(size = 4, color = col_season) +
  geom_segment(aes(y = lowerCL, x = season, yend = upperCL, xend = season),
               size = 2, color = col_season, lineend = "round") +
  labs(y=" ") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        text = element_text(size = 10))


# Analysis: clustering --------------------------------------------------------------

data_med4season = data_longform_input %>% group_by(ID,SeasonNumber) %>%
  mutate(MedResult = mean(Result)) %>%
  select(ID,SeasonNumber,Result, MedResult,
         Gender,Environment,Age,AgeEntrance,
         Days,t_ji) %>%
  distinct(ID,SeasonNumber, .keep_all = TRUE) %>%
  ungroup() %>% arrange(SeasonNumber)

data_med4season$Clustering = rep(0,nrow(data_med4season))
for(idx in 1:nrow(data_med4season)){
  id = data_med4season$ID[idx]
  season_num = data_med4season$SeasonNumber[idx]
  nobs = which(dt$ID_i == id)
  data_med4season$Clustering[idx] = dt$finalPartition[[season_num]][[nobs]]
}


data_with_clustering = data_longform_input %>%
  left_join(data_med4season %>%
              select("ID","SeasonNumber","MedResult","Clustering"),
            by = c("ID","SeasonNumber"))


# Global statistics on pooled data

NMale   = data_with_clustering %>% distinct(SeasonNumber,ID, .keep_all = T) %>%
            filter(Gender == "M") %>% summarise(count = n() ) %>% pull(count)
NFemale = data_with_clustering %>% distinct(SeasonNumber,ID, .keep_all = T) %>%
            filter(Gender == "W") %>% summarise(count = n() ) %>% pull(count)
MeanMale   = data_with_clustering %>% filter(Gender == "M") %>%
              summarise(mean = mean(Result) ) %>% pull(mean)
MeanFemale = data_with_clustering %>% filter(Gender == "W") %>%
              summarise(mean = mean(Result) ) %>% pull(mean)
VarMale   = data_with_clustering %>% filter(Gender == "M") %>% summarise(var = var(Result) ) %>% pull(var)
VarFemale = data_with_clustering %>% filter(Gender == "W") %>% summarise(var = var(Result) ) %>% pull(var)


AgeMale   = data_with_clustering %>% filter(Gender == "M") %>% summarise(age = mean(Age) ) %>% pull(age)
AgeFemale = data_with_clustering %>% filter(Gender == "W") %>% summarise(age = mean(Age) ) %>% pull(age)

AoeMale   = data_with_clustering %>% filter(Gender == "M") %>% summarise(aoe = mean(AgeEntrance) ) %>% pull(aoe)
AoeFemale = data_with_clustering %>% filter(Gender == "W") %>% summarise(aoe = mean(AgeEntrance) ) %>% pull(aoe)


Nclus = length(table(data_with_clustering$Clustering))
cluster_summary = tibble(Cluster = as.integer(),
                         Nputs = as.integer(),
                         NMales = as.integer(),NFemales = as.integer(),
                         MeanMale = as.numeric(), MeanFemale = as.numeric(),
                         MeanMaleFemalecl = as.numeric(),
                         VarMale = as.numeric(), VarFemale = as.numeric(),
                         AgeMale = as.numeric(), AgeFemale = as.numeric(),
                         AoeMale = as.numeric(), AoeFemale = as.numeric())

for(cl in 1:Nclus){
  # Statistics within each cluster
  temp = data_with_clustering %>% filter(Clustering == cl)
  size = table(data_with_clustering$Clustering)[cl]
  NMalecl = temp %>% distinct(SeasonNumber,ID, .keep_all = T) %>%
              filter(Gender == "M") %>% summarise(count = n() ) %>% pull(count)
  NFemalecl = temp %>% distinct(SeasonNumber,ID, .keep_all = T) %>%
                filter(Gender == "W") %>% summarise(count = n() ) %>% pull(count)
  MeanMalecl   = temp %>% filter(Gender == "M") %>% summarise(mean = mean(Result) ) %>% pull(mean)
  MeanFemalecl = temp %>% filter(Gender == "W") %>% summarise(mean = mean(Result) ) %>% pull(mean)
  MeanMaleFemalecl = mean(c(MeanMalecl,MeanFemalecl))

  VarMalecl   = temp %>% filter(Gender == "M") %>% summarise(var = var(Result) ) %>% pull(var)
  VarFemalecl = temp %>% filter(Gender == "W") %>% summarise(var = var(Result) ) %>% pull(var)

  AgeMalecl   = temp %>% filter(Gender == "M") %>% summarise(age = mean(Age) ) %>% pull(age)
  AgeFemalecl = temp %>% filter(Gender == "W") %>% summarise(age = mean(Age) ) %>% pull(age)

  AoeMalecl   = temp %>% filter(Gender == "M") %>% summarise(aoe = mean(AgeEntrance) ) %>% pull(aoe)
  AoeFemalecl = temp %>% filter(Gender == "W") %>% summarise(aoe = mean(AgeEntrance) ) %>% pull(aoe)


  cluster_summary = cluster_summary %>%
    rbind( tibble(Cluster = cl, Nputs = size,
                  NMales = NMalecl, NFemales = NFemalecl,
                  MeanMale = MeanMalecl, MeanFemale = MeanFemalecl, MeanMaleFemalecl = MeanMaleFemalecl,
                  VarMale = VarMalecl, VarFemale = VarFemalecl,
                  AgeMale = AgeMalecl, AgeFemale = AgeFemalecl,
                  AgeEntranceMale = AoeMalecl, AgeEntranceFemale = AoeFemalecl  ) )
}

cluster_summary$Cluster = as.factor(cluster_summary$Cluster)

# Define cluster ordering according to decreasing mean male
cluster_summary = cluster_summary %>% arrange( desc(MeanMaleFemalecl))
# Define cluster ordering according to decreasing mean male/female
# cluster_summary %>% arrange( desc(MeanMaleFemalecl))
# --> only minor changes happen and only due to noisy clusters
old_labels = cluster_summary$Cluster
new_labels = as.factor(1:nlevels(cluster_summary$Cluster))
cluster_summary = cluster_summary %>% mutate(Cluster = new_labels)

# new_labels <- recode(cluster_summary$Cluster,
#                      `1` = 10, `2` = 9, `3` = 3, `4` = 8, `5` = 11,
#                      `6` = 5 , `7` = 2 , `8` = 7, `9` = 12, `10` = 4,
#                      `11`= 1 , `12` = 6, `13` = 15, `14` = 13, `15` = 14)

cluster_summary =  rbind( tibble(Cluster = "Pooled", Nputs = sum(n_j),
                                 NMales = NMale, NFemales = NFemale,
                                 MeanMale = MeanMale, MeanFemale = MeanFemale,
                                 MeanMaleFemalecl = MeanMaleFemalecl,
                                 VarMale = VarMale, VarFemale = VarFemale,
                                 AgeMale = AgeMale, AgeFemale = AgeFemale,
                                 AgeEntranceMale = AoeMale, AgeEntranceFemale = AoeFemale),
                          cluster_summary)

cluster_summary = cluster_summary %>% mutate(Size = NMales+NFemales) %>% mutate(across(5:14, ~ round(., digits = 2))) %>% select(Cluster,Size,everything())

cluster_summary_means = cluster_summary[,c(1:10)]
cluster_summary_ages = cluster_summary[,c(1,10:13)]

kable(cluster_summary_means, caption = "Cluster Interpretation - Means and Variances")
# kable(cluster_summary_ages, caption = "Cluster Interpretation - Ages")




# Global clusters' sizes

# Compute similarity matrix
recode_map <- setNames(new_labels, old_labels)
VI_sara$cl <- recode(VI_sara$cl, !!!recode_map)
table(VI_sara$cl)

# Local number of clusters and associated traceplots
Kj_VI = vector(length = d)
idx_start = c(1,cumsum(n_j)[1:(d-1)]+1)
idx_end = cumsum(n_j)
for(j in 1:d){
  Tab_j = table(VI_sara$cl[idx_start[j]:idx_end[j]])
  Kj_VI[j] = length( which(Tab_j>0) )
}

cat("Number of local cluster: \n")
Kj_VI


# Now, given the 15 global clusters ordered as above, I want to print their sizes in each season


Nclus = length(table(data_with_clustering$Clustering))
Local_sizes = matrix(0, nrow = Nclus, ncol = d)

idx_start = c(1,cumsum(n_j)[1:(d-1)]+1)
idx_end = cumsum(n_j)
for(j in 1:d){
  Local_table = table(VI_sara$cl[idx_start[j]:idx_end[j]])
  Local_sizes[as.numeric(names(Local_table)), j] = Local_table
}

rownames(Local_sizes) = unlist(lapply(list(1:Nclus), function(k){paste0("Cluster ",k)}))
colnames(Local_sizes) = unlist(lapply(list(1:d), function(j){paste0("S",j)}))

kable(Local_sizes, caption = "Cluster sizes across different seasons. Each row represents a cluster, each column represents a season")


# Local sizes plots ----------------------------------------------------

n_j_seasons = apply(Local_sizes,2,sum)
Nseason = 15
Nclus = length(table(data_with_clustering$Clustering))
mycol_clus = hcl.colors(n = Nclus, palette = "Temps")
Local_sizes_perc = Local_sizes
for(jj in 1:ncol(Local_sizes)){
  Local_sizes_perc[,jj] = Local_sizes_perc[,jj]/n_j_seasons[jj]
}


par(mar = c(2,2,2,1), mfrow = c(1,1), bty = "l")
matplot(t(Local_sizes[1:4,]), type = "b", pch = 16, lty = 1, col = mycol_clus[1:4])
par(mar = c(2,2,2,1), mfrow = c(1,1), bty = "l")
matplot(t(Local_sizes[5:8,]), type = "b", pch = 16, lty = 1, col = mycol_clus[5:8])
par(mar = c(2,2,2,1), mfrow = c(1,1), bty = "l")
matplot(t(Local_sizes[9:12,]), type = "b", pch = 16, lty = 1, col = mycol_clus[9:12])
par(mar = c(2,2,2,1), mfrow = c(1,1), bty = "l")
matplot(t(Local_sizes[13:15,]), type = "b", pch = 16, lty = 1, col = mycol_clus[13:15])

par(mar = c(2,2,2,1), mfrow = c(1,1), bty = "l")
matplot(t(Local_sizes_perc[1:4,]), type = "b", pch = 16, lty = 1, col = mycol_clus[1:4])
par(mar = c(2,2,2,1), mfrow = c(1,1), bty = "l")
matplot(t(Local_sizes_perc[5:8,]), type = "b", pch = 16, lty = 1, col = mycol_clus[5:8])
par(mar = c(2,2,2,1), mfrow = c(1,1), bty = "l")
matplot(t(Local_sizes_perc[9:12,]), type = "b", pch = 16, lty = 1, col = mycol_clus[9:12])
par(mar = c(2,2,2,1), mfrow = c(1,1), bty = "l")
matplot(t(Local_sizes_perc[13:15,]), type = "b", pch = 16, lty = 1, col = mycol_clus[13:15])

# Final clustering  Visualization -----------------------------------------

counter_obs = 1
for(j in 1:d){
  for(i in 1:n){
    if(dt$N_ji[j,i] > 0){
      dt$finalPartition[[j]][[i]] = VI_sara$cl[counter_obs]
      counter_obs = counter_obs + 1
    }
  }
}


old2 = factor(data_with_clustering$Clustering)
old2 <- recode(old2, !!!recode_map, .default = as.factor("16"))
levels(old2)
any(is.na(old2))
data_with_clustering$Clustering = old2



Nclus = length(table(data_with_clustering$Clustering))
mycol_clus = hcl.colors(n = Nclus, palette = "Temps")

seasons = 1:d
cl_plots = 1:Nclus #c(1:6,8:10,12:13)

par(mar = c(2,2,2,1), mfrow = c(1,1), bty = "l")
plot( x = 0, y = 0, type = "n",
      ylab = "Result", xlab = "Season",
      main = "Male athletes",
      xlim = c(0,d), ylim = c(-6,4.5))
abline(v = seasons, lty = 2, col = "grey45", lwd = 1)
for( cl in cl_plots ){
  temp = data_with_clustering %>% filter(Clustering == cl)%>% filter(Gender == "M")
  #par(mar = c(4,4,2,1))
  points( x = temp$t_ji,
          y = temp$Result,
          pch = 16, cex = 0.3, col = mycol_clus[cl])
}



# Then, I repeat the previous plot **female athletes**, which  is coherent with the one for male.

Nclus = length(table(data_with_clustering$Clustering))
mycol_clus = hcl.colors(n = Nclus, palette = "Temps")

seasons = 1:d
cl_plots = 1:Nclus #c(1:6,8:10,12:13)

par(mar = c(2,2,2,1), mfrow = c(1,1), bty = "l")
plot( x = 0, y = 0, type = "n",
      ylab = "Result", xlab = "Season",
      main = "Female athletes",
      xlim = c(0,d), ylim = c(-6,4.5))
abline(v = seasons, lty = 2, col = "grey45", lwd = 1)
for( cl in cl_plots ){
  temp = data_with_clustering %>% filter(Clustering == cl)%>% filter(Gender == "W")
  #par(mar = c(4,4,2,1))
  points( x = temp$t_ji,
          y = temp$Result,
          pch = 16, cex = 0.3, col = mycol_clus[cl])
}
# highlight women in cluster 14
cl = 14
temp = data_with_clustering %>% filter(Clustering == cl)%>% filter(Gender == "W")
points( x = temp$t_ji,
        y = temp$Result,
        pch = 8, cex = 0.4, col = mycol_clus[15])



# Extra - Cluster specific plots ------------------------------------------

# Male

seasons = 1:d
cl_plots = 1:Nclus #c(1:6,8:10,12:13)

par(mar = c(2,2,2,1), mfrow = c(4,4), bty = "l")
for( cl in cl_plots ){
plot( x = 0, y = 0, type = "n",
      ylab = "Result", xlab = "Season",
      main = "Male athletes",
      xlim = c(0,d), ylim = c(-6,4.5))
  abline(v = seasons, lty = 2, col = "grey45", lwd = 1)
  temp = data_with_clustering %>% filter(Clustering == cl)%>% filter(Gender == "M")
  #par(mar = c(4,4,2,1))
  points( x = temp$t_ji,
          y = temp$Result,
          pch = 16, cex = 0.3, col = mycol_clus[cl])
}

# Female

seasons = 1:d
cl_plots = 1:Nclus #c(1:6,8:10,12:13)

par(mar = c(2,2,2,1), mfrow = c(4,4), bty = "l")
for( cl in cl_plots ){
  plot( x = 0, y = 0, type = "n",
        ylab = "Result", xlab = "Season",
        main = "Female athletes",
        xlim = c(0,d), ylim = c(-6,4.5))
  abline(v = seasons, lty = 2, col = "grey45", lwd = 1)
  temp = data_with_clustering %>% filter(Clustering == cl)%>% filter(Gender == "W")
  #par(mar = c(4,4,2,1))
  points( x = temp$t_ji,
          y = temp$Result,
          pch = 16, cex = 0.3, col = mycol_clus[cl])
}

# Extra - Season specific plots ------------------------------------------

# Male

seasons = 1:d
cl_plots = 1:Nclus #c(1:6,8:10,12:13)

par(mar = c(2,2,2,1), mfrow = c(4,4), bty = "l")
for( j in seasons ){
  plot( x = 0, y = 0, type = "n",
        ylab = "Result", xlab = "Season",
        main = "Male athletes",
        xlim = c(0,d), ylim = c(-6,4.5))
  abline(v = seasons, lty = 2, col = "grey45", lwd = 1)
  temp = data_with_clustering %>% filter(SeasonNumber == j)%>% filter(Gender == "M")
  #par(mar = c(4,4,2,1))
  points( x = temp$t_ji,
          y = temp$Result,
          pch = 16, cex = 0.3, col = mycol_clus[temp$Clustering])
}

# Female

seasons = 1:d
cl_plots = 1:Nclus #c(1:6,8:10,12:13)

par(mar = c(2,2,2,1), mfrow = c(4,4), bty = "l")
for( j in seasons ){
  plot( x = 0, y = 0, type = "n",
        ylab = "Result", xlab = "Season",
        main = "Female athletes",
        xlim = c(0,d), ylim = c(-6,4.5))
  abline(v = seasons, lty = 2, col = "grey45", lwd = 1)
  temp = data_with_clustering %>% filter(SeasonNumber == j)%>% filter(Gender == "W")
  #par(mar = c(4,4,2,1))
  points( x = temp$t_ji,
          y = temp$Result,
          pch = 16, cex = 0.3, col = mycol_clus[temp$Clustering])
}

# Trajectories ------------------------------------------------------------

# plot some trajectories c(32,110,214,254)
seasons = 1:d
ID_plyrs = data_longform_input %>% distinct(ID) %>% pull(ID)
yrange = c(-4.7,3.5)
pt_size = 1 #1.3 (for slide/poster)

mycol = hcl.colors(n=3,palette = "Zissou1")
par(mar = c(4,4,3,1), mfrow = c(1,1))

ID_ply = ID_plyrs[32]
temp = data_longform_input %>% filter(ID == ID_ply) %>% arrange(t_ji)
plot( x = temp$t_ji, y = temp$Result,
      ylab = "Result", xlab = "Season",
      main = paste0("Athlete - ",ID_ply),
      xlim = c(0,15),
      ylim = yrange,
      pch = 16, cex = pt_size, col = "black")
abline(v = seasons, lty = 2, col = "grey45", lwd = 1)


ID_ply = ID_plyrs[110]
temp = data_longform_input %>% filter(ID == ID_ply) %>% arrange(t_ji)
par(mar = c(4,4,3,1), mfrow = c(1,1))
plot( x = temp$t_ji, y = temp$Result,
      ylab = "Result", xlab = "Season",
      main = paste0("Athlete - ",ID_ply),
      xlim = c(0,15),
      ylim = yrange,
      pch = 16, cex = pt_size, col = "black")
abline(v = seasons, lty = 2, col = "grey45", lwd = 1)


ID_ply = ID_plyrs[214]
temp = data_longform_input %>% filter(ID == ID_ply) %>% arrange(t_ji)
par(mar = c(4,4,3,1), mfrow = c(1,1))
plot( x = temp$t_ji, y = temp$Result,
      ylab = "Result", xlab = "Season",
      main = paste0("Athlete - ",ID_ply),
      xlim = c(0,15),
      ylim = yrange,
      pch = 16, cex = pt_size, col = "black")
abline(v = seasons, lty = 2, col = "grey45", lwd = 1)


ID_ply = ID_plyrs[254]
temp = data_longform_input %>% filter(ID == ID_ply) %>% arrange(t_ji)
par(mar = c(4,4,3,1), mfrow = c(1,1))
plot( x = temp$t_ji, y = temp$Result,
      ylab = "Result", xlab = "Season",
      main = paste0("Athlete - ",ID_ply),
      xlim = c(0,15),
      ylim = yrange,
      pch = 16, cex = pt_size, col = "black")
abline(v = seasons, lty = 2, col = "grey45", lwd = 1)




# Analysis: trajectories --------------------------------------------------




# posterior plot of trajectories with bands
ID_plyrs = data_with_ordered_clustering %>% distinct(ID) %>% pull(ID)
raf_ply = c(32,110,214,254)
mycol_clus = hcl.colors(n = Nclus, palette = "Temps")
mycol_clus[12] = mycol_clus[15]

for(ii in raf_ply){
  ID_ply = ID_plyrs[ii]
  curva_ply = predictive_players(ID_ply = ID_ply, dt = application_result$dt,
                                 fit = application_result$GDFMM,
                                 burnin = niter/2 )
  n_ji = nrow(curva_ply)
  temp = data_with_ordered_clustering %>% filter(ID == ID_ply) %>% arrange(t_ji)

  mycol = hcl.colors(n=3,palette = "Zissou1")
  par(mar = c(4,4,2,1), mfrow = c(1,1))
  plot( x = 0, y = 0,
        ylab = "Result", xlab = "Season",
        main = paste0("Athlete - ",ID_ply),
        xlim = range(data_with_ordered_clustering$t_ji),
        ylim = yrange,
        pch = 16, cex = pt_size, type = "n")
  abline(v = seasons, lty = 2, col = "grey45", lwd = 1)
  for(j in 1:n_ji){
    x_ji = temp %>% filter(SeasonNumber == j) %>% pull(t_ji)
    y_ji = temp %>% filter(SeasonNumber == j) %>% pull(Result)
    cl_ji = temp %>% filter(SeasonNumber == j) %>% pull(Clustering)
    points( x = x_ji, y = y_ji, col = mycol_clus[cl_ji[1]],
            pch = 16, cex = pt_size )
    segments(x0 = j-1, x1 = j,
             y0 = curva_ply[j,2], y1 = curva_ply[j,2],
             mycol[1], lwd = 2 )
    rect(xleft = j-1, xright = j,
         ybottom = curva_ply[j,1], ytop = curva_ply[j,3],
         col = ACutils::t_col(mycol[1], percent = 85),border = NA )
  }


}










