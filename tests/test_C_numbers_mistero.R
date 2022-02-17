my_logC(n = 4, scale = -exp(2), location = 0)
lgfactorial_ns(n = 4, sigma = -exp(2))
my_logC2(n = 4, scale = -exp(2), location = 0)
my_logC2_central(n = 4, scale = -exp(2))
compute_logC(n = 4, scale = -exp(2), location = 0)


start_time1 <- Sys.time()
for( i in 1:1000){
  A = my_logC(n = 500, scale = -1, location = 0)
}
end_time1 <- Sys.time()
end_time1 - start_time1

start_time2 <- Sys.time()
for(i in 1:1000){

  a = calcola_stirling_ricor_log(n = 500, gamma = 1, r = 0)

}
end_time2 <- Sys.time()
end_time2 - start_time2

start_time3 <- Sys.time()
for(i in 1:1000){

  r = lgfactorial_ns(n = 500, sigma = -1)

}
end_time3 <- Sys.time()
end_time3 - start_time3



start_time4 <- Sys.time()
for(i in 1:1000){

  ac2 = my_logC2(n = 500, scale = -1, location = 0)

}
end_time4 <- Sys.time()
end_time4 - start_time4


start_time5 <- Sys.time()
for(i in 1:1000){

  ac_c = my_logC2_central(n = 500, scale = -1)

}
end_time5 <- Sys.time()
end_time5 - start_time5



start_time6 <- Sys.time()
for(i in 1:1000){

  ac_c3 = compute_logC(n = 500, scale = -1, location = 0)

}
end_time6 <- Sys.time()
end_time6 - start_time6


