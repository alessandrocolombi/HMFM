#  Libraries ---------------------------------------------------------------

suppressWarnings(suppressPackageStartupMessages(library(GDFMM)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(salso)))
suppressWarnings(suppressPackageStartupMessages(library(mcclust.ext)))
suppressWarnings(suppressPackageStartupMessages(library(abind)))
suppressWarnings(suppressPackageStartupMessages(library(knitr)))

setwd(here::here())
data("ShotPutData")
setwd("./ShotPutApplication")
seed0 = 29071996
set.seed(seed0)
# Load data ---------------------------------------------------------------
data_longform_input = ShotPutData

# Center data and covariates
data_longform_input$Result      = data_longform_input$Result      - mean(data_longform_input$Result)

# Create a pooled ID
data_longform_input$ID_ji = as.factor(paste0(as.character(data_longform_input$ID),"_",data_longform_input$SeasonNumber))
data_longform_input$pooled = as.factor("1")
data_longform_input$Partition0 = 1


# Define dt object --------------------------------------------------------


dt = input_handle(data_longform_input[,c(10,11,3,12,4)], intercept = FALSE)


n = dt$n
d = dt$d
r = dt$r
n_j = dt$n_j

# Hyperparameters: P0 ---------------------------------------------------------

# Res_range = range(data_longform_input$Result)
Res_range = quantile(data_longform_input$Result, probs = c(0.005,0.995))
R = Res_range[2] - Res_range[1]
mu0 = mean(data_longform_input$Result) # should be 0
k0  = 1/R^2
nu0 = 4
sigma0 = 20/2

scale = sqrt( (k0 + 1)/(k0) * sigma0 )
mean_marginal = mu0
var_marginal  = nu0/(nu0-2) * scale^2

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

cluster_mean = 0
cluster_var  = data_longform_input %>% summarise(ClusterVar = var(Result)) %>% pull(ClusterVar)
initial_partition = unlist(unlist(dt$initialPartition))


option = set_options( "mu0" = mu0,"sigma0"= sigma0, "k0"= k0, "nu0"=nu0,
                      "Adapt_MH_hyp2" = 0.234, "Adapt_MH_var0"=0.1,
                      "proposal_Mstar" = 1,
                      "Lambda0" = Lambda0, "gamma0" = gamma0, "Mstar0" = Mstar0,
                      "beta0" = beta0, "Sigma0" = Sigma0,
                      "alpha_gamma" = a_gamma, "beta_gamma" = b_gamma,
                      "alpha_lambda" = a_lambda, "beta_lambda" = b_lambda,
                      "init_mean_cluster" = c(cluster_mean, rep(0,Mstar0)),
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

VI_sara = minVI(sim_matrix)
SeasonNumber_all = data_longform_input %>% distinct(ID_ji, SeasonNumber = SeasonNumber) %>% pull(SeasonNumber)

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

data_med4season = data_longform_input %>% group_by(ID_ji) %>% mutate(MedResult = mean(Result)) %>%
  select(ID,SeasonNumber,ID_ji,Result, MedResult,Gender,Environment,Age,AgeEntrance,Days,t_ji) %>%
  distinct(ID_ji, .keep_all = TRUE) %>%
  ungroup() %>% arrange(SeasonNumber)


data_med4season$Clustering = rep(0,nrow(data_med4season))
for(idx in 1:nrow(data_med4season)){
  id_ji = data_med4season$ID_ji[idx]
  # id    = data_med4season$ID[idx]
  # season_num = data_med4season$SeasonNumber[idx]
  nobs = which(dt$ID_i == id_ji)
  data_med4season$Clustering[idx] = dt$finalPartition[[1]][[nobs]]
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
  if(is.na(MeanMalecl))
    MeanMalecl = 0
  if(is.na(MeanFemalecl))
    MeanFemalecl = 0

  MeanMaleFemalecl = mean(temp$Result)
  # MeanMaleFemalecl = (NMalecl*MeanMalecl + NFemalecl*MeanFemalecl )/(NMalecl+NFemalecl)

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
# cluster_summary %>% arrange( desc(MeanMale))
# --> they coincide

old_labels = cluster_summary$Cluster
new_labels = as.factor(1:nlevels(cluster_summary$Cluster))
cluster_summary = cluster_summary %>% mutate(Cluster = new_labels)

cluster_summary =  rbind( tibble(Cluster = "Pooled", Nputs = sum(n_j),
                                 NMales = NMale, NFemales = NFemale,
                                 MeanMale = MeanMale, MeanFemale = MeanFemale,
                                 MeanMaleFemalecl = (NMale*MeanMale + NFemale*MeanFemale )/(NMale+NFemale),
                                 VarMale = VarMale, VarFemale = VarFemale,
                                 AgeMale = AgeMale, AgeFemale = AgeFemale,
                                 AgeEntranceMale = AoeMale, AgeEntranceFemale = AoeFemale),
                          cluster_summary)

cluster_summary = cluster_summary %>% mutate(Size = NMales+NFemales) %>% mutate(across(5:13, ~ round(., digits = 2))) %>% select(Cluster,Size,everything())

cluster_summary_means = cluster_summary[,c(1:10)]
cluster_summary_ages = cluster_summary[,c(1,10:13)]

kable(cluster_summary_means, caption = "Cluster Interpretation - Means and Variances")
# kable(cluster_summary_ages, caption = "Cluster Interpretation - Ages")


paper_summary =  cbind(cluster_summary_means[,c(1,8,4:7,9:10)],
                       cluster_summary_ages[,c(3,4)])
kable(paper_summary, caption = "Cluster Interpretation - Paper table")


# Recode lables
recode_map <- setNames(new_labels, old_labels)
VI_sara$cl <- recode(VI_sara$cl, !!!recode_map)
table(VI_sara$cl)



# Now, given the 15 global clusters ordered as above, I want to print their sizes in each season


Nseason = 15
Nclus = length(table(data_with_clustering$Clustering))
Local_sizes = matrix(0, nrow = Nclus, ncol = Nseason)
for(j in 1:Nseason){
  idx_season_j = which(SeasonNumber_all == j)
  Local_table = table(VI_sara$cl[idx_season_j])
  Local_sizes[as.numeric(names(Local_table)), j] = Local_table
}

rownames(Local_sizes) = unlist(lapply(list(1:Nclus), function(k){paste0("Cluster ",k)}))
colnames(Local_sizes) = unlist(lapply(list(1:Nseason), function(j){paste0("S",j)}))

kable(Local_sizes, caption = "Cluster sizes across different seasons. Each row represents a cluster, each column represents a season")



Fake_Kj = apply(Local_sizes,2,function(x){length(which(x > 0))})
cat("Number of local cluster: \n")
Fake_Kj

# Plot cluster sizes evolution
library(plot.matrix)
Local_sizes_plot = Local_sizes
rownames(Local_sizes_plot) = as.character(1:Nclus)
colnames(Local_sizes_plot) = as.character(1:Nseason)
my_breaks = c(0,1,10,30,50,70,
              120,140,253 )
n_breaks  = length(my_breaks)
my_col_mat = hcl.colors(n = n_breaks - 2, palette = "Heat 2", rev = TRUE)
my_col_mat = c("white",my_col_mat)

# par(mar = c(2.75,2.75,0.5,2), mgp = c(1.75,0.75,0))
# plot(Local_sizes_plot,
#      breaks = my_breaks, col = my_col_mat,
#      xlab = "Season", ylab = "Cluster", main = "",
#      cex.axis = 0.7, axes = FALSE,
#      key=list(side=4, cex.axis=0.75), axis.key=NULL, spacing.key=0.75, fmt.key="%.0f",
# )

par(mar = c(3.5,3.5,0.5,3), mgp = c(2,0.55,0))
plot(Local_sizes_plot,
     breaks = my_breaks, col = my_col_mat,
     xlab = "Season", ylab = "Cluster", main = "",
     cex.axis = 1, cex.lab = 1.75, axes = FALSE,
     key=list(side=4, cex.axis=1.5), axis.key=NULL, spacing.key=0.75, fmt.key="%.0f",
)


# Final clustering  Visualization -----------------------------------------



# I must also recode labels in data_with_clustering
old2 = factor(data_with_clustering$Clustering)
old2 <- recode(old2, !!!recode_map, .default = as.factor("16"))
levels(old2)
any(is.na(old2))
data_with_clustering$Clustering = old2


Nclus = length(table(data_with_clustering$Clustering))
mycol_clus = hcl.colors(n = Nclus, palette = "Temps")

seasons = 1:15 #1:d
cl_plots = 1:Nclus

par(mar = c(2,2,2,1), mfrow = c(1,1), bty = "l")
plot( x = 0, y = 0, type = "n",
      ylab = "Result", xlab = "Season",
      main = "Male athletes",
      xlim = c(0,15), ylim = c(-6,4.5))
abline(v = seasons, lty = 2, col = "grey45", lwd = 1)
for( cl in cl_plots ){
  temp = data_with_clustering %>% filter(Clustering == cl) %>% filter(Gender == "M")
  #par(mar = c(4,4,2,1))
  points( x = temp$t_ji,
          y = temp$Result,
          pch = 16, cex = 0.3, col = mycol_clus[cl])
}



# Then, I repeat the previous plot **female athletes**, which  is coherent with the one for male.

Nclus = length(table(data_with_clustering$Clustering))
mycol_clus = hcl.colors(n = Nclus, palette = "Temps")

seasons = 1:15
cl_plots = 1:Nclus #c(1:6,8:10,12:13)

par(mar = c(2,2,2,1), mfrow = c(1,1), bty = "l")
plot( x = 0, y = 0, type = "n",
      ylab = "Result", xlab = "Season",
      main = "Female athletes",
      xlim = c(0,15), ylim = c(-6,4.5))
abline(v = seasons, lty = 2, col = "grey45", lwd = 1)
for( cl in cl_plots ){
  temp = data_with_clustering %>% filter(Clustering == cl)%>% filter(Gender == "W")
  #par(mar = c(4,4,2,1))
  points( x = temp$t_ji,
          y = temp$Result,
          pch = 16, cex = 0.3, col = mycol_clus[cl])
}



# Peak histogram --------------------------------------------------------------------

ID_plyrs = data_with_clustering %>% distinct(ID) %>% pull(ID)
peak = c()
for(i in seq_along(ID_plyrs)){
  ID_i = ID_plyrs[i]
  temp_i = data_with_clustering %>% filter(ID == ID_i)
  temp_i = temp_i %>% filter(Clustering == min(as.integer(temp_i$Clustering)))
  seasons_i = as.integer(temp_i %>% pull(SeasonNumber))
  peak = c(peak,min(seasons_i))
}

val_peak = rep(0,Nseason)
val_peak[as.numeric(names(table(peak)))] = table(peak)
peak_tibble = tibble( "Season" = as.factor(1:Nseason), "val" = as.integer(val_peak) )

col_season = "grey47"
ggplot(peak_tibble, aes(x = Season, y = val)) +
  geom_bar(stat = "identity", fill = col_season) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none",
        text = element_text(size = 10)) +
  labs(y=" ", x = "Season")



