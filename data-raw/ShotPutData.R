## code to prepare `ShotPutData` dataset goes here

suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
data_aligned = read.csv("./data-raw/ShotPutData.csv", row.names=1)
data_aligned = as_tibble(data_aligned) %>%
               mutate(ID = as.factor(ID),
                      SeasonNumber = as.factor(SeasonNumber),
                      Gender = as.factor(Gender),
                      Environment = as.factor(Environment)) %>%
               select(ID,SeasonNumber,Result,
                      Gender,Environment,Age,AgeEntrance,
                      Days,t_ji)


# select first 15 seasons only
n = nrow(data_aligned %>% distinct(ID))
selectIDs = 1:n
d = 15
IDs  = data_aligned %>% distinct(ID) %>% pull(ID)
ShotPutData = data_aligned %>%
              filter(ID %in% IDs[selectIDs]) %>%
              filter(SeasonNumber %in% as.factor(1:d)) %>%
              filter(ID != "76011") %>%
              select(ID,SeasonNumber,Result,Gender,Environment,Age,AgeEntrance,Days,t_ji)

usethis::use_data(ShotPutData, overwrite = TRUE)
