# LIBRARIES #
library(bbricks)
library(dplyr)
library(salso)

#
#k0=1 / (max(dat, na.rm = T) - min(dat, na.rm = T)) ^ 2 --> paper
#interpret hyper ->hoff non abbiamo dati precedenti media 0 su campione 200 dati
#per far coincidere con i dati dele ragazze nu0 deve essere 1/2 e sigma0 0.5^2

diagplot<-function (param, color, title, it){
  plot(param,
       type = "l",
       col = color,
       ylab= paste(title, " ", j),
       main="")
  acf(param, lag = 100,main="")
  mtext(title,
        side = 3,
        line = - 1,
        outer = TRUE,
        font = 4 )
}


# Data pre processing -----------------------------------------------------


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

d=3
ncol_data <- max(length(data_level1), length(data_level2),length(data_level3))
dat <- matrix(NA, nrow = d, ncol = ncol_data)
dat[1, 1:length(data_level1)] <- data_level1
dat[2, 1:length(data_level2)] <- data_level2
dat[3, 1:length(data_level3)] <- data_level3




# Gibbs Sampler 1st run ---------------------------------------------------


niter<-2500
burnin<-5000


option<-list("Mstar0" =10,"Lambda0"=2,"mu0"=0,"nu0"=200,"sigma0"=1^2,
             "Adapt_MH_hyp1"=0.7,"Adapt_MH_hyp2"=0.234, "Adapt_MH_power_lim"=10, "Adapt_MH_var0"=1,
             "k0"= 1/sqrt(200), "alpha_gamma"=1,
             "beta_gamma"=1, "alpha_lambda"=1, "beta_lambda"=1)
GDFMM = GDFMM_sampler(dat, niter, burnin, 2, seed = 123, option = option)


# COMPUTE BINDER LOSS FUNCTION TO SELECT BEST PARTITION -------------------

part_matrix <- GDFMM$Partition

sim_matrix <- psm(part_matrix)
# VI_dahl <- dlso(matr, loss = 'VI', estimate=NULL)
binder_dahl <- dlso(part_matrix, loss = 'binder', estimate = sim_matrix)

estimate_partition = as.vector(binder_dahl)


# Gibbs Sampler second run ------------------------------------------------



option_fixed <- list("Mstar0" =0,"Lambda0"=2,"mu0"=0,"nu0"=200,"sigma0"=1^2,
                     "Adapt_MH_hyp1"=0.7,"Adapt_MH_hyp2"=0.234, "Adapt_MH_power_lim"=10, "Adapt_MH_var0"=1,
                     "k0"= 1/sqrt(200), "alpha_gamma"=1,
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
diagplot(GDFMM$K, "black")
#M
diagplot(GDFMM$M, "black")
#lambda
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
  diagplot(GDFMM_fixed$mu[[j]], "black")
}
#sigma
x11()
par(mfcol = c(2, 3))
for (j in 1:length(GDFMM_fixed$sigma)[1]){
  diagplot(GDFMM_fixed$sigma[[j]], "black")
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



table(GDFMM$K)




# Simulated Data testing --------------------------------------------------





