
# Librerie
suppressWarnings(suppressPackageStartupMessages(library(GDFMM)))
suppressWarnings(suppressPackageStartupMessages(library(ACutils)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(library(salso)))
suppressWarnings(suppressPackageStartupMessages(library(wesanderson)))
suppressWarnings(suppressPackageStartupMessages(library(mcclust.ext)))
suppressWarnings(suppressPackageStartupMessages(library(abind)))

# Custom functions

my_invgamma = function(x, nu0, sigma0){
  a = nu0/2
  b = a*sigma0
  b^a/gamma(a) * (1/x)^(a+1) * exp(-b/x)
}

# Esperimento

#Qua carico tutti i dati a disposizione
data_aligned = read.csv("Shotput_longformat_preproc_all.csv", row.names=1)
data_aligned = as_tibble(data_aligned) %>% mutate(ID = as.factor(ID),
                                                  SeasonNumber = as.factor(SeasonNumber),
                                                  Gender = as.factor(Gender),
                                                  Environment = as.factor(Environment)) %>%
  select(ID,SeasonNumber,Result,
         Gender,Environment,Age,AgeEntrance,
         Days,t_ji)
data_aligned

n = nrow(data_aligned %>% distinct(ID))



#Seleziono solo pochi dati e solo un numero ristretto di stagioni
#L'atleta con ID 76011 è quello che nella stagione 3 ha fatto un punteggio ridicolo (1.48). Però ha comunque giocato molte partite e molte stagioni, posso decidere se togliere solo quella stagione (però poi ho un buco) oppure togliere lui completamente
selectIDs = 1:n#c(1:10,800:810)#
IDs  = data_aligned %>% distinct(ID) %>% pull(ID)
data_longform_input = data_aligned %>%
                      filter(ID %in% IDs[selectIDs]) %>%
                      filter(SeasonNumber %in% as.factor(1:4)) %>%
                      filter(ID != "76011") %>%
                      select(ID,SeasonNumber,Result,Gender,Environment,Age,AgeEntrance,Days,t_ji)

#Trasformo i dati - Centro i dati e le covariate quantitative

data_longform_input$Result      = data_longform_input$Result * 100
data_longform_input$Result      = data_longform_input$Result      - mean(data_longform_input$Result)
data_longform_input$Age         = data_longform_input$Age         - mean(data_longform_input$Age)
data_longform_input$AgeEntrance = data_longform_input$AgeEntrance - mean(data_longform_input$AgeEntrance)


# Leggo i dati nel modo in cui poi devono essere inseriti
dt = input_handle(data_longform_input[,c(1:4)], intercept = FALSE)

n = dt$n
d = dt$d
r = dt$r
n_j = dt$n_j

#save(dt, file = "dt_d19_data51415.Rdat")

#View(dt)
n;d;r;sum(n_j);sum(dt$N_ji)

## Plot

mycol_gender = c("deepskyblue3", "palevioletred1")
cols = data_longform_input$Gender
levels(cols) = mycol_gender
cols = as.character(cols)

seasons = 1:d

par(mar = c(4,4,2,1))
plot( x = data_longform_input$t_ji,
      y = data_longform_input$Result,
      ylab = "Result", xlab = "Season",
      pch = 16, cex = 0.4, col = cols)
abline(v = seasons, lty = 2, col = "grey45", lwd = 1)


Res_range = range( data_longform_input$Result )
mu0 = mean(data_longform_input$Result) # should be 0
k0  = 1/(Res_range[2] - Res_range[1])^2
nu0 = 4
sigma0 = 0.5

c("mu0"=mu0,"k0"=k0,"nu0"=nu0,"sigma0"=sigma0)

scale = sqrt( (k0 + 1)/(k0) * sigma0 )
mean_marginal = mu0
var_marginal  = nu0/(nu0-2) * scale^2

mean_marginal
var_marginal

xrange = c(-10000,10000)
l_grid = 10000
grid = seq(xrange[1],xrange[2],length.out = l_grid)

#par(mfrow = c(1,3), mar = c(2,2,2,1), bty = "l")
# Qua plotto la marginale
Prior_grid = GDFMM:::dnct(x = grid, n0 = nu0, mu0 = mu0, gamma0 = scale)
par(mar = c(2,2,2,1), bty = "l")
plot(1,1, xlim = xrange, ylim = c(0,max(Prior_grid)),
     type = "n",
     main = "Marginal data Density")
grid(lty = 1,lwd = 2,col = "gray90" )
points(x = na.omit(as.numeric(data_longform_input$Result)),
       y = rep(0, length(na.omit(as.numeric(data_longform_input$Result)))),
       ylim = c(0,0.1),
       pch = 16)
points(grid, Prior_grid, col = "black", lwd = 2, type = "l")

# Qua plotto la inverse-gamma
empricial_var = sum(dt$N_ji * dt$var_ji)/sum(dt$N_ji)
sigma_grid = seq(1e-5, 10, length.out = 1000)
par(mar = c(2,2,2,1), bty = "l")
plot(x = sigma_grid, y = my_invgamma(sigma_grid, nu0 = nu0, sigma0 = sigma0),
     main = "Inverse-gamma prior Density", type = "l", lwd = 2, lty = 1)
abline(v = empricial_var, lwd = 2, lty = 2, col = "grey45")
grid(lty = 1,lwd = 2,col = "gray90" )
legend("topright", legend = c( "empricial var", paste0("InvGamma(",nu0/2,", ",sigma0*nu0/2,")")  ),
       lwd = 2, lty = c(2,1), col = c("grey45","black"))

# Qua invece plotto la normale, condizionata a dei valori di sigma estratti
Invgamma_grid = 1/rgamma(n=10000, shape = nu0/2, rate = sigma0*nu0/2)
media = rnorm(n=10000,mu0, sqrt(Invgamma_grid/k0))
par(mar = c(2,2,2,1), bty = "l")
plot(density(media),main = "Normal prior Density", lwd = 2)
grid(lty = 1,lwd = 2,col = "gray90" )



## Run
#Kobj = 4
niter  <-  20000#0
burnin <-  10000#0
thin   <-       1

mu_gamma  =  0.1
var_gamma =  1
a_gamma = (mu_gamma*mu_gamma)/var_gamma
b_gamma = mu_gamma/(var_gamma)
#a_gamma <- 2 ; b_gamma  <- 0.005

mu_lambda  = 10
var_lambda = 3
a_lambda = (mu_lambda*mu_lambda)/var_lambda
b_lambda = mu_lambda/(var_lambda)

# initial values
beta0 = rep(0,dt$r)
Sigma0 = 100*diag(dt$r)


# load the 'mlrMBO' package
library(mlrMBO)

# STEP 1: definition of the objective function
obj.fun <- makeSingleObjectiveFunction(

  # give a name to the objective function
  name = "GDFMM-Nclus",

  # implement the objective function
  fn = function(x) {
    Lambda0 = x[1]
    gamma0  = x[2:5]

    option = set_options( "mu0" = mu0,"sigma0"= sigma0, "k0"= k0, "nu0"=nu0,
                          "Adapt_MH_hyp2" = 0.234, "Adapt_MH_var0"=0.001,
                          "proposal_Mstar" = 1,
                          "Lambda0" = Lambda0, "gamma0" = gamma0, "Mstar0" = 1,
                          "beta0" = beta0, "Sigma0" = Sigma0,
                          "alpha_gamma" = a_gamma, "beta_gamma" = b_gamma,
                          "alpha_lambda" = a_lambda, "beta_lambda" = b_lambda,
                          "init_mean_cluster" = NULL,
                          "init_var_cluster" = NULL,
                          "partition" = rep(0,sum(n_j)),#NULL, #seq(1:sum(n_j)),
                          "IncludeCovariates" = TRUE,
                          "UpdateU" = T, "UpdateM" = T, "UpdateGamma" = T,
                          "UpdateS" = T, "UpdateTau" = T, "UpdateLambda" = T,
                          "UpdateBeta" = T )

    prior = "Normal"
    #prior = "Normal-InvGamma"
    GDFMM = ConditionalSampler(dt, niter, burnin, thin, seed = 123, option = option, FixPartition = F,
                               P0.prior = prior)



    ### Clustering
    # Get labels for each iterations for each data point
    part_matrix <- GDFMM$Partition #GDFMM$Partition is a (n_iter x n_data) matrix

    # Compute similarity matrix
    sim_matrix <- psm(part_matrix)
    VI_sara = minVI(sim_matrix)
    Nclus = length(table(VI_sara$cl))
    #return( abs(Nclus - Kobj))
    Nclus
  },

  # define if the objective function has to be minimized or maximized (i.e., accuracy must be maximized)
  minimize = F,

  # define the search space
  par.set = makeParamSet(
    makeNumericParam( "Lambda", lower=0.05, upper=10 ),
    makeNumericVectorParam("gamma",4,lower = 0.0001,upper = 1)
  )

)


# STEP 2: generation of the initial design
set.seed(123452)
des = generateDesign( n=100, getParamSet(obj.fun), fun=lhs::randomLHS )
des$y = apply( des, 1, obj.fun )

des$y
#
# # STEP 3: define the Probabilistic Surrogate Model
#
# # Gaussian Process (regression)
# #psm = makeLearner( "regr.bgp", predict.type="se", covtype="exp", control=list(trace=F))
# # covtype: "gauss", "exp", "powexp", "matern3_2", "matern5_2
#
#
# # STEP 4: Sequential process and acquisition
# control = makeMBOControl()
# control = setMBOControlTermination( control, iters=10 )
# # acquisition function
# #control = setMBOControlInfill( control, crit=makeMBOInfillCritEI() ) # EI
#  control = setMBOControlInfill( control, crit=makeMBOInfillCritCB() ) # CB
#
# # start sequential BO
# run = mbo( obj.fun, design=des, learner=NULL, control=control, show.info=T )
#
#
# # retrieving and plotting results
# best.seen <- getOptPathY(run$opt.path)
# best.seen <- c( max(best.seen[1:5]), best.seen[6:15] )
# # we are maximizing accuracy
# plot( cummax(best.seen), type="o", lwd=3, col="blue", ylim=c(min(best.seen),max(best.seen)),
#       ylab="best seen", xlab="trials")
# lines( best.seen, type="o", lty=2, col="green", lwd=3 )
# legend( "bottomright", legend=c("best seen","Accuracy[i]"), col=c("blue","green"), lty=1:2, lwd=3, pch=1 )
#
#
#
#
#
