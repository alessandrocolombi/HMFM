
# Test one group
n = 104
tot = 0
for(i in 0:n){
  pp = p_distinct_prior(k=i, n_groups = n, gamma = 2, prior = "Poisson", lambda = 2, Max_iter = 1000)
  print(pp)
  tot = tot + pp
}
tot
# Esplode per n>=104


# Test 2 groups
n1 = 2;n2 = 2
n = c(n1,n2)
tot = 0
for(i in 0:sum(n)){
  pp = p_distinct_prior(k=i, n_j = n, gamma = c(2,2), prior = "Poisson", lambda = 2, Max_iter = 1000)
  print(pp)
  tot = tot + pp
}
tot

# Se una delle n è > 100 esplode
0.2073875 + 0.4756042 + 0.2754024 + 0.04160585




library(GDFMM)
p_distinct_prior(k=1, n_j = c(1,1), gamma = c(2,2), prior = "Poisson", lambda = 2, Max_iter = 1000)
p_distinct_prior(k=1, n_j = c(1), gamma = c(2), prior = "Poisson", lambda = 2, Max_iter = 1000)

exp(1.38629)*exp(0.693147)*exp(-3.68111)
1.38629+0.693147
p_distinct_prior(k=3, n_j = c(1,1,1), gamma = c(2,2,2), prior = "Poisson", lambda = 2, Max_iter = 1000)

0.2492789 + 0.5491604 + 0.4031214



exp( log( choose(1,3-2) ) + my_log_falling_factorial(3-1,0) )

