
# Test one group
m = 1
n = 10
k = 3
tot = 0
for(i in 0:m){
  pp = p_distinct_posterior(r=i,k=k, m_j = m, n_j = n, gamma = 2, prior = "Poisson", lambda = 2, Max_iter = 100)
  print(pp)
  tot = tot + pp
}
tot


# Test 2 groups
m1 = 1;m2 = 1
n1 = 10;n2 = 10
n = c(n1,n2)
m = c(m1,m2)
k = 3
tot = 0
for(i in 0:sum(m)){
  pp = p_distinct_posterior(r=i,k=k, m_j = m, n_j = n, gamma = rep(2,length(n)), prior = "Poisson", lambda = 2, Max_iter = 100)
  print(pp)
  tot = tot + pp
}
tot


# Test 3 groups
m1 = 1;m2 = 1;m3=1
n1 = 10;n2 = 10;n3=10
n = c(n1,n2,n3)
m = c(m1,m2,m3)
k = 3
tot = 0
for(i in 0:sum(m)){
  pp = p_distinct_posterior(r=i,k=k, m_j = m, n_j = n, gamma = rep(2,length(n)), prior = "Poisson", lambda = 2, Max_iter = 100)
  print(pp)
  tot = tot + pp
}
tot

# Test d groups
d = 4
n_max = 10
n = sample(1:n_max,size = d,replace = T)
m = rep(1,d)
k = 3
n
tot = 0
for(i in 0:sum(m)){
  pp = p_distinct_posterior(r=i,k=k, m_j = m, n_j = n, gamma = rep(2,length(n)), prior = "Poisson", lambda = 2, Max_iter = 100)
  print(pp)
  tot = tot + pp
}
tot



library(GDFMM)
p_distinct_posterior(r=0,k=2, m_j = c(1,1), n_j = c(2,2), gamma = c(2,2), prior = "Poisson", lambda = 2, Max_iter = 100)

0.7816969 + 0.1995299 +0.01877326
compute_logC(  1, -2, - ( 1*2 + 1 )  )



Vprior = exp(-0.693147)
sum_M = 0
r = 1
for(M in r:100){
  sum_M = sum_M + ((factorial(M+1))/factorial(M-r)) * ((2^M)*exp(-2))/(factorial(M)*( 2*(M+1)+1 )*(2*(M+1)))
}
val = sum_M/Vprior
# val(r=0) = 0.1700015
# val(r=1) = 0.2449977
# --> log_Vprior = -0.693147
# se r=0 ---> sum_M = 0.08500075, log(sum_M) = -2.465095, Vpost = 0.1700015, log_Vpost = -1.771948
# se r=1 ---> sum_M = 0.1224989,  log(sum_M) = -2.099653, Vpost = 0.2449977, log_Vpost = -1.406506


prob_0 = val*exp(1.0986123)
# prob(r=0) = 0.5100044
prob_1 = val*exp(0.6931472)
# prob(r=1) = 0.4899954




sum_M = 0
r = 0
for(M in r:100){
  sum_M = sum_M + ( (2^M)*exp(-2))/(factorial(M) )
}
sum_M
log(sum_M)
log(sum_M) - log(Vprior)





p_distinct_posterior(r=0,k=2, m_j = c(1,1,0), n_j = c(2,2,0), gamma = c(2,2,2), prior = "Poisson", lambda = 2, Max_iter = 100)

0.9160297 + 0.07794079 + 0.005827313 + 0.000202199




