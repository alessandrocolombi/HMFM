
# Real Dataset ------------------------------------------------------------


# Libraries ---------------------------------------------------------------

library(bbricks)
library(dplyr)
library(salso)
library(ggplot2)


# Preprocessing of the data -----------------------------------------------


data2<-data(hlrData)
data1<-data(hlrData)
df <- data.frame(hlrData$mathScore, hlrData$socioeconomicStatus, hlrData$schoolID)
names(df) <- c('math_score', 'soc_status', 'school_id')
count <- df %>% count(school_id)
count_ordered <- count[order(count$n, decreasing = T),]
data1 <- matrix(NA, nrow = 100, ncol = max(count$n))
dim(data1)
d = 100
n_j <- rep(0, d)
for (j in 1:d) {
  data_level <- subset(df, school_id == j)
  data_level <- data_level$math_score
  n_j[j] <- length(data_level)
  data1[j, 1:n_j[j]] <- data_level
}
data2<-data1[c(67,51,79,92,100,94,5,72,17),]



# Gibbs sampler 1st run ---------------------------------------------------

option2 <-list("Mstar0"= 5,"Lambda0"= 2,"mu0" = mean(data2, na.rm = T), "nu0"= 100,"sigma0"= 400,
               "Adapt_MH_hyp1"=0.7,"Adapt_MH_hyp2"=0.234, "Adapt_MH_power_lim"=10,  "Adapt_MH_var0"=1,
               "k0"= 100, "alpha_gamma"=1, "beta_gamma"=1, "alpha_lambda"=1, "beta_lambda"=1)

GS = GDFMM_sampler(data1,2500,2500,2,seed=123, option=option2)


# Get partition by binder loss function minimization ----------------------


part_matrix <- GS$Partition

sim_matrix <- psm(part_matrix)
# VI_dahl <- dlso(matr, loss = 'VI', estimate=NULL)
binder_dahl <- dlso(part_matrix, loss = 'binder', estimate = sim_matrix)

estimate_partition = as.vector(binder_dahl)


# Gibbs Sampler 2nd run ---------------------------------------------------

option2_fixed <-list("Mstar0"= 0,"Lambda0"=2,"mu0"=mean(data2, na.rm = T),"nu0"= 150,"sigma0"=40,
               "Adapt_MH_hyp1"=0.7,"Adapt_MH_hyp2"=0.234, "Adapt_MH_power_lim"=10,  "Adapt_MH_var0"=1,
               "k0"= 1/sqrt(2500), "alpha_gamma"=1, "beta_gamma"=1, "alpha_lambda"=1, "beta_lambda"=1,
               "partition" = estimate_partition)

GS_fixed = GDFMM_sampler(data1, 500, 500, 2, seed=123, FixPartition = T, option=option2_fixed)



# Plot for Pres -----------------------------------------------------------
k<-na.omit(as.vector(data1[c(67,51,79,92,100,94,5,72,17),]))

df2 <- df1[which(df1$school_id %in% c(67,51,79,92,100,94,5,72,17)) , ]

df_plot <- df %>%
  group_by(school_id) %>%
  mutate(avgmathscore = mean(math_score), maxms=max(math_score), minms=min(math_score))%>%as.data.frame()
df_plot <- df_plot[order(df_plot$avgmathscore),]
df_plot$rank <- as.integer(as.factor(df_plot$avgmathscore))

clust = as.factor(estimate_partition +1)

ggplot(df_plot, aes(x=rank,y=math_score,group=school_id, ymin=minms,
                    ymax=maxms, group = clust)) +
  geom_point( aes(color = clust, shape = clust), size =3 ) +
  geom_linerange() +
  labs(x = "Rank - Math score average", y = "Math Score", title ="2002 US Math Score")


l[order(l$math_score),]


