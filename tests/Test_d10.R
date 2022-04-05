
y1_m1 = rnorm(20,-3, 1/2) # 1st level, 1st comp
y1_m2 = rnorm(20, 0, 1/2) # 1st level, 2nd comp

y2_m1 = rnorm(200,-3, 1/2) # 2nd level, 1st comp
y2_m2 = rnorm(200, 3, 1/2) # 2nd level, 2nd comp

y3_m1 = rnorm(200,-3, 1/2) # 3nd level, 1st comp
y3_m2 = rnorm(20, 3, 1/2) # 3nd level, 2nd comp
y3_m3 = rnorm(200, 0, 1/2) # 3nd level, 3rd comp

y4_m1 = rnorm(20,-3, 1/2) # 1st level, 1st comp
y4_m2 = rnorm(20, 0, 1/2) # 1st level, 2nd comp

y5_m1 = rnorm(200,-3, 1/2) # 2nd level, 1st comp
y5_m2 = rnorm(200, 3, 1/2) # 2nd level, 2nd comp

y6_m1 = rnorm(200,-3, 1/2) # 3nd level, 1st comp
y6_m2 = rnorm(20, 3, 1/2) # 3nd level, 2nd comp
y6_m3 = rnorm(200, 0, 1/2) # 3nd level, 3rd comp

y7_m1 = rnorm(20,-3, 1/2) # 1st level, 1st comp
y7_m2 = rnorm(20, 0, 1/2) # 1st level, 2nd comp

y8_m1 = rnorm(200,-3, 1/2) # 2nd level, 1st comp
y8_m2 = rnorm(200, 3, 1/2) # 2nd level, 2nd comp

y9_m1 = rnorm(200,-3, 1/2) # 3nd level, 1st comp
y9_m2 = rnorm(20, 3, 1/2) # 3nd level, 2nd comp
y9_m3 = rnorm(200, 0, 1/2) # 3nd level, 3rd comp

y10_m1 = rnorm(20,-3, 1/2) # 1st level, 1st comp
y10_m2 = rnorm(20, 0, 1/2) # 1st level, 2nd comp

y11_m1 = rnorm(200,-3, 1/2) # 2nd level, 1st comp
y11_m2 = rnorm(200, 3, 1/2) # 2nd level, 2nd comp

y12_m1 = rnorm(200,-3, 1/2) # 3nd level, 1st comp
y12_m2 = rnorm(20, 3, 1/2) # 3nd level, 2nd comp
y12_m3 = rnorm(200, 0, 1/2) # 3nd level, 3rd comp



real_partition = c(rep(0, 20), rep(1, 20), rep(0, 200), rep(2, 200), rep(0, 200),
                   rep(2, 50), rep(1, 200) )

data_level1 <- c(y1_m1, y1_m2)
data_level2 <- c(y2_m1, y2_m2)
data_level3 <- c(y3_m1, y3_m2, y3_m3)

data_level4 <- c(y4_m1, y4_m2)
data_level5 <- c(y5_m1, y5_m2)
data_level6 <- c(y6_m1, y6_m2, y6_m3)

data_level7 <- c(y7_m1, y7_m2)
data_level8 <- c(y8_m1, y8_m2)
data_level9 <- c(y9_m1, y9_m2, y9_m3)

data_level10 <- c(y10_m1, y10_m2)
data_level11 <- c(y11_m1, y11_m2)
data_level12 <- c(y12_m1, y12_m2, y12_m3)



data_all <- c(data_level1, data_level2,data_level3,
              data_level4, data_level5,data_level6,
              data_level7, data_level7,data_level9,
              data_level8, data_level11,data_level12)


# Data pre-processing -------------------------------------------------------
d=12
ncol_data <- max(length(data_level1), length(data_level2),length(data_level3),
                 length(data_level4), length(data_level5),length(data_level6),
                 length(data_level7), length(data_level8),length(data_level9),
                 length(data_level10), length(data_level11),length(data_level12))

dat <- matrix(NA, nrow = d, ncol = ncol_data)
dat[1, 1:length(data_level1)] <- data_level1
dat[2, 1:length(data_level2)] <- data_level2
dat[3, 1:length(data_level3)] <- data_level3
dat[4, 1:length(data_level4)] <- data_level4
dat[5, 1:length(data_level5)] <- data_level5
dat[6, 1:length(data_level6)] <- data_level6
dat[7, 1:length(data_level7)] <- data_level7
dat[8, 1:length(data_level8)] <- data_level8
dat[9, 1:length(data_level9)] <- data_level9
dat[10, 1:length(data_level10)] <- data_level10
dat[11, 1:length(data_level11)] <- data_level11
dat[12, 1:length(data_level12)] <- data_level12

# Gibbs Sampler 1st run ---------------------------------------------------

niter <-5000
burnin <- 1
thin <- 1

option<-list("Mstar0" = 50,"Lambda0" = 5,"mu0" = 0,"nu0"=10,"sigma0"= 1,
             "Adapt_MH_hyp1"= 0.7,"Adapt_MH_hyp2"= 0.234, "Adapt_MH_power_lim"=10, "Adapt_MH_var0"=1,
             "k0"= 1/10, "alpha_gamma"=1,
             "beta_gamma"=1, "alpha_lambda"=1, "beta_lambda"=1)

GDFMM = GDFMM_sampler(dat, niter, burnin, thin, seed = 123, option = option)

# Interessante questo caso con Mstar0 = 25.
# Inizialmente mixa bene, è che quando arriva Mstar=0 poi non si muove più. In particolare, interessante
# il traceplot di lambda


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
diagplot(GDFMM$lambda, "black", "lambda")

#gamma
x11()
par(mfcol = c(2, 3))
for (j in 1:dim(GDFMM$gamma)[1]){
  diagplot(GDFMM$gamma[j,], "black", "Gamma", j)
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



# U
U = GDFMM$U
x11();matplot(t(U), type = 'l')

log_sum_u = apply(log(U), 2, sum)

x11();plot(log_sum_u, type = 'l')
x11();plot(exp(-log_sum_u), type = 'l')

min(exp(log_sum_u))
