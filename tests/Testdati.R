# LIBRARIES #
library(bbricks)
library(dplyr)
library(salso)
library(LaplacesDemon)

#-------------------------------------------------------------------------------
#
#k0=1 / (max(dat, na.rm = T) - min(dat, na.rm = T)) ^ 2 --> paper
#interpret hyper ->hoff non abbiamo dati precedenti media 0 su campione 200 dati
#per far coincidere con i dati dele ragazze nu0 deve essere 1/2 e sigma0 0.5^2

diagplot<-function (param, color, title = "title", it = 0){

  if(it == 0){
    y_lab = paste(title)
  }else{
    y_lab = paste(title, "[", it,"]")
  }

  plot(param,
       type = "l",
       col = color,
       ylab= y_lab,
       main="")
  acf(param, lag = 100,main="")
  mtext(title,
        side = 3,
        line = - 1,
        outer = TRUE,
        font = 4 )
}


# Data generation -----------------------------------------------------

y1_m1 = rnorm(20,-3, 1/2) # 1st level, 1st comp
y1_m2 = rnorm(20, 0, 1/2) # 1st level, 2nd comp
y2_m1 = rnorm(200,-3, 1/2) # 2nd level, 1st comp
y2_m2 = rnorm(200, 3, 1/2) # 2nd level, 2nd comp
y3_m1 = rnorm(200,-3, 1/2) # 3nd level, 1st comp
y3_m2 = rnorm(50, 3, 1/2) # 3nd level, 2nd comp
y3_m3 = rnorm(200, 0, 1/2) # 3nd level, 3rd comp

real_partition = c(rep(0, 20), rep(1, 20), rep(0, 200), rep(2, 200), rep(0, 200),
                   rep(2, 50), rep(1, 200) )

data_level1 <- c(y1_m1, y1_m2)
data_level2 <- c(y2_m1, y2_m2)
data_level3 <- c(y3_m1, y3_m2, y3_m3)

comp_m3 <- c(y1_m1, y2_m1, y3_m1)
comp_0 <- c(y1_m2, y3_m3)
comp_3 <- c(y2_m2, y3_m2)

data_all <- c(data_level1, data_level2,data_level3)


# Plots -------------------------------------------------------------------


plot = F
if(plot){
    x11()
    par(mfrow = c(2,2))
    plot(density(data_level1),main="Density - Level 1", col = 'salmon', lwd = 3,
         xlim = c(-20,20))
    plot(density(data_level2),main="Density - Level 2", col = 'aquamarine2', lwd = 3,
         xlim = c(-20,20))
    plot(density(data_level3),main="Density - Level 3", col = 'aquamarine2', lwd = 3,
         xlim = c(-20,20))
    plot(density(data_all),main="Density - all levels", col = 'purple', lwd = 3,
         xlim = c(-20,20))
}

# Data pre-processing -------------------------------------------------------
d=3
ncol_data <- max(length(data_level1), length(data_level2),length(data_level3))
dat <- matrix(NA, nrow = d, ncol = ncol_data)
dat[1, 1:length(data_level1)] <- data_level1
dat[2, 1:length(data_level2)] <- data_level2
dat[3, 1:length(data_level3)] <- data_level3


# Gibbs Sampler 1st run ---------------------------------------------------

niter <-2000
burnin <- 1000
thin <- 2

option<-list("Mstar0" =10,"Lambda0"=2,"mu0"=0,"nu0"=10,"sigma0"=1^2,
             "Adapt_MH_hyp1"=0.7,"Adapt_MH_hyp2"=0.234, "Adapt_MH_power_lim"=10, "Adapt_MH_var0"=1,
             "k0"= 1/14, "alpha_gamma"=1,
             "beta_gamma"=1, "alpha_lambda"=1, "beta_lambda"=1)

GDFMM = GDFMM_sampler(dat, niter, burnin, thin, seed = 123, option = option)

# COMPUTE BINDER LOSS FUNCTION TO SELECT BEST PARTITION -------------------

part_matrix <- GDFMM$Partition

sim_matrix <- psm(part_matrix)
# VI_dahl <- dlso(matr, loss = 'VI', estimate=NULL)
binder_dahl <- dlso(part_matrix, loss = 'binder', estimate = sim_matrix)

estimate_partition = as.vector(binder_dahl)

# misclassified = sum(estimate_partition != real_partition)
# cat("Number of misclassified object = ", misclassified, " ; accuracy = ",
#     misclassified/length(real_partition) )

# Gibbs Sampler second run ------------------------------------------------

option_fixed <- list("Mstar0" = 100, "Lambda0"= 2, "mu0"=0, "nu0"=10, "sigma0"= 1,
                     "Adapt_MH_hyp1"=0.7,"Adapt_MH_hyp2"=0.234, "Adapt_MH_power_lim"=10, "Adapt_MH_var0"=1,
                     "k0"= 1/14, "alpha_gamma"=1,
                     "beta_gamma"=1, "alpha_lambda"=1, "beta_lambda"=1, "partition" = estimate_partition)

GDFMM_fixed = GDFMM_sampler(dat, niter, burnin, 2, seed = 123, FixPartition = T, option = option_fixed)


# Inference on parameters -------------------------------------------------


table(GDFMM$K)

plot_GDFMM = T
if(plot_GDFMM){
  x11()
  par(mfrow = c(1,3))
  hist(GDFMM$Mstar)
  hist(GDFMM$K)
  hist(GDFMM$lambda)
}

mu<-lapply(GDFMM_fixed$mu, mean)
sigma<-lapply(GDFMM_fixed$mu, var)


# Diagnostic plot ---------------------------------------------------------


#K
x11()
par(mfrow = c(2,1))
diagplot(GDFMM$K, "black", "K")

#M
x11()
par(mfrow = c(2,1))
diagplot(GDFMM$M, "black", "M")

#lambda
x11()
par(mfrow = c(2,1))
diagplot(GDFMM$lambda, "black")

#gamma
x11()
par(mfcol = c(2, 3))
for (j in 1:dim(GDFMM_fixed$gamma)[1]){
  diagplot(GDFMM_fixed$gamma[j,], "black", "Gamma", j)
}

#mu
x11()
par(mfcol = c(2, 3))
for (j in 1:length(GDFMM_fixed$mu)[1]){
  diagplot(GDFMM_fixed$mu[[j]], "black", "mu", j)
}
#sigma
x11()
par(mfcol = c(2, 3))
for (j in 1:length(GDFMM_fixed$sigma)[1]){
  diagplot(GDFMM_fixed$sigma[[j]], "black", "sigma", j)
}


# Density plot ------------------------------------------------------------

#Numerosità per gruppo 0->249 1->221 ->420

#Group Numerosity
tab<-table(estimate_partition)

tab[1]

tab1<-table(estimate_partition[1:40])

tab2<-table(estimate_partition[41:640])

tab3<-table(estimate_partition[641:890])

liv<-list()
simul<-c()

for (i in 1:length(mu)){
  liv[[i]] <-rnorm(tab[i],mu[[i]],sigma[[i]])
  simul<-c(simul,rnorm(tab[i],mu[[i]],sigma[[i]]))
}



liv11<-rnorm(tab1[[1]],mu[[1]],sigma[[1]])
liv13<-rnorm(tab1[[2]],mu[[3]],sigma[[3]])

simul1<-c(liv11,liv13)

liv22<-rnorm(tab2[[1]],mu[[2]],sigma[[2]])
liv23<-rnorm(tab2[[2]],mu[[3]],sigma[[3]])
simul2<-c(liv22,liv23)

liv31<-rnorm(tab3[[1]],mu[[1]],sigma[[1]])
liv32<-rnorm(tab3[[2]],mu[[2]],sigma[[2]])
liv33<-rnorm(tab3[[3]],mu[[3]],sigma[[3]])
simul3<-c(liv31,liv32,liv33)

x11()
plot(density(simul))

plot = T
if(plot){
  x11()
  par(mfrow = c(2,2))
  plot(density(data_level1),main="Density - Level 1", col = 'salmon', lwd = 3,
       xlim = c(-20,20))
  lines(density(simul1), lty="dashed", lwd=2)
  plot(density(data_level2),main="Density - Level 2", col = 'aquamarine2', lwd = 3,
       xlim = c(-20,20))
  lines(density(simul2), lty="dashed", lwd=2)
  plot(density(data_level3),main="Density - Level 3", col = 'aquamarine2', lwd = 3,
       xlim = c(-20,20))
  lines(density(simul3), lty="dashed", lwd=2)
  plot(density(data_all),main="Density - all levels", col = 'purple', lwd = 3,
       xlim = c(-20,20))
    lines(density(simul), lty="dashed", lwd=2)
}

# Density plot V2 --------------------------------------------

grid = seq(min(data_all)-1, max(data_all)+1, by = 0.001)

w_ji = GDFMM_fixed$w_ji
w_1 = rowMeans(w_ji[[1]])
w_2 = rowMeans(w_ji[[2]])
w_3 = rowMeans(w_ji[[3]])

mu = GDFMM_fixed$mu
mu_1 = mean(mu[[1]])
mu_2 = mean(mu[[2]])
mu_3 = mean(mu[[3]])
mu_vec = c(mu_1, mu_2, mu_3)

sigma = GDFMM_fixed$sigma
sigma1 = sqrt(mean(sigma[[1]]))
sigma2 = sqrt(mean(sigma[[2]]))
sigma3 = sqrt(mean(sigma[[3]]))
sigma_vec = c(sigma1, sigma2, sigma3)

mix_1 = dmix(grid, w_1, mu_vec, sigma_vec)
mix_2 = dmix(grid, w_2, mu_vec, sigma_vec)
mix_3 = dmix(grid, w_3, mu_vec, sigma_vec)
# mix_all = dmix(grid, 1/3*(w_1 + w_2 + w_3), mu_vec, sigma_vec)

# PLOT OF MIXTURES IN THE GROUPS

x11()
par(mfrow = c(2,2))

plot(grid, dmix(grid,c(1/2, 1/2), c(0, -3), c(1/2, 1/2)),
     main = "Density comparison for group 1", ylab = "Group 1",
     type = "l", lty = 2, col = "gray", lwd = 2)
lines(density(data_level1), lty = 2, col = "red", lwd = 2)
lines(grid, mix_1, lwd = 2, col = "blue")
legend("topright", legend = c("real", "sample", "estimate"),
       col = c("gray", "red", "blue"), lty = c(2, 2, 1), lwd = 2)

plot(grid, dmix(grid,c(1/2, 1/2), c(3, -3), c(1/2, 1/2)),
     main = "Density comparison for group 2", ylab = "Group 2",
     type = "l", lty = 2, col = "gray", lwd = 2)
lines(density(data_level2), lty = 2, col = "red", lwd = 2)
lines(grid, mix_2, lwd = 2, col = "blue")
legend("top", legend = c("real", "sample", "estimate"),
       col = c("gray", "red", "blue"), lty = c(2, 2, 1), lwd = 2)

plot(grid, dmix(grid, c(4/9, 1/9, 4/9), c(-3, 3, 0), c(1/2, 1/2, 1/2)),
     main = "Density comparison for group 3", ylab = "Group 3",
     type = "l", lty = 2, col = "gray", lwd = 2)
lines(density(data_level3), lty = 2, col = "red", lwd = 2)
lines(grid, mix_3, lwd = 2, col = "blue")
legend("topright", legend = c("real", "sample", "estimate"),
       col = c("gray", "red", "blue"), lty = c(2, 2, 1), lwd = 2)


# PLOT OF COMPONENTS
x11()
par(mfrow = c(2,2))

d_comp1 = GDFMM::dnorm_est(grid, GDFMM_fixed$mu[[1]], sqrt(GDFMM_fixed$sigma[[1]]))

plot(density(comp_m3), main = "Estimate vs Sample density - mu=-3", lty = 2, col = "gray")
plot(grid, dnorm(grid, -3, 1/2), main = "Estimate vs Real density - mu=3", type = "l", lty = 2, col = "gray")
lines(grid, d_comp1$Inf. , lty = 4, col = "red")
lines(grid, d_comp1$Sup. , lty = 4, col = "red")
lines(grid, d_comp1$Est. ,  col = "red", lwd = 2)


d_comp2 = GDFMM::dnorm_est(grid, GDFMM_fixed$mu[[2]], sqrt(GDFMM_fixed$sigma[[2]]))

plot(density(comp_0), main = "Estimate vs Sample density - mu=0", lty = 2, col = "gray")
plot(grid, dnorm(grid, 0, 1/2), main = "Estimate vs Real density - mu=3",
     type = "l", lty = 2, col = "gray")
lines(grid, d_comp2$Inf. , lty = 4, col = "red")
lines(grid, d_comp2$Sup. , lty = 4, col = "red")
lines(grid, d_comp2$Est. ,  col = "red", lwd = 2)

d_comp3 = GDFMM::dnorm_est(grid, GDFMM_fixed$mu[[3]], sqrt(GDFMM_fixed$sigma[[3]]))

plot(density(comp_3), main = "Estimate vs Sample density - mu=3", lty = 2, col = "gray")
plot(grid, dnorm(grid, 3, 1/2), main = "Estimate vs Real density - mu=3", type = "l", lty = 2, col = "gray")
lines(grid, d_comp3$Inf. , lty = 4, col = "red")
lines(grid, d_comp3$Sup. , lty = 4, col = "red")
lines(grid, d_comp3$Est. ,  col = "red", lwd = 2)


# Simulated Data testing --------------------------------------------------

LaplacesDemon::rnormm()



