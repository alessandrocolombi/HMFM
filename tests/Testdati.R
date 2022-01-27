#datapreprocessing
library(bbricks)
library(dplyr)
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
option2 <-list("Mstar0"=3,"Lambda0"=2,"mu0"=mean(data, na.rm = T),"nu0"=2.5,"sigma0"=var(data, na.rm=T),
             "Adapt_MH_hyp1"=0.7,"Adapt_MH_hyp2"=0.234, "Adapt_MH_power_lim"=10, "k0"=1, "alpha_gamma"=1,
             "beta_gamma"=1, "alpha_lambda"=1, "beta_lambda"=1)
example_GDFMM_sampler(data2,6,6,2,123, option2)

y1_m1 = rnorm(20,-3, 1/2) # 1st level, 1st comp
y1_m2 = rnorm(20, 0, 1/2) # 1st level, 2nd comp
y2_m1 = rnorm(200,-3, 1/2) # 2nd level, 1st comp
y2_m2 = rnorm(200, 3, 1/2) # 2nd level, 2nd comp
y3_m1 = rnorm(200,-3, 1/2) # 3nd level, 1st comp
y3_m2 = rnorm(50, 3, 1/2) # 3nd level, 2nd comp
y3_m3 = rnorm(200, 0, 1/2) # 3nd level, 3rd comp
data_level1 <- c(y1_m1, y1_m2)
data_level2 <- c(y2_m1, y2_m2)
data_level3 <- c(y3_m1, y3_m2, y3_m3)
data_all <- c(data_level1, data_level2,data_level3)

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

d=3
ncol_data <- max(length(data_level1), length(data_level2),length(data_level3))
dat <- matrix(NA, nrow = d, ncol = ncol_data)
dat[1, 1:length(data_level1)] <- data_level1
dat[2, 1:length(data_level2)] <- data_level2
dat[3, 1:length(data_level3)] <- data_level3
option<-list("Mstar0"=2,"Lambda0"=2,"mu0"=mean(dat, na.rm = T),"nu0"=1,"sigma0"=0.5,
             "Adapt_MH_hyp1"=0.7,"Adapt_MH_hyp2"=0.234, "Adapt_MH_power_lim"=10,
             "k0"= 1 / (max(dat, na.rm = T) - min(dat, na.rm = T)) ^ 2, "alpha_gamma"=1,
             "beta_gamma"=1, "alpha_lambda"=1, "beta_lambda"=1, "Adapt_MH_var0"=1)
example_GDFMM_sampler(dat,10,10,2,123,option = option)

df1<-df %>%
  group_by(school_id) %>%
  mutate(avgmathscore = mean(math_score), maxms=max(math_score), minms=min(math_score))%>%as.data.frame()
df1<-df1[order(df1$avgmathscore),]
df1$rank<-as.integer(as.factor(df1$avgmathscore))

ggplot(df1, aes(x=rank,y=math_score,group=school_id, ymin=minms, ymax=maxms)) +
  geom_point()+geom_linerange()+labs(x = "Rank - Math score average", y = "Math Score",
                                       title ="2002 US Math Score")

library(dplyr)






l[order(l$math_score),]

