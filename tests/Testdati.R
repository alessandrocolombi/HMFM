#datapreprocessing
library(bbricks)
library(dplyr)
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


example_GDFMM_sampler(data1,6000,2000,2)



