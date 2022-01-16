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

example_GDFMM_sampler(data1,10,10,2,123)

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

