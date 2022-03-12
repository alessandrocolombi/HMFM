
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
n1 = 100;n2 = 100
n = c(n1,n2)
tot = 0
for(i in 0:sum(n)){
  pp = p_distinct_prior(k=i, n_groups = n, gamma = c(2,2), prior = "Poisson", lambda = 2, Max_iter = 1000)
  #print(pp)
  tot = tot + pp
}
tot

# Se una delle n è > 100 esplode







Test_multiple_groups(k=2, n_groups = c(1,2,3), gamma = c(2,2,2), prior = "Poisson", lambda = 2, Max_iter = 1000)









