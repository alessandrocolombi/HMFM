
# Test one group
n = 2
tot = 0
for(i in 0:n){
  pp = p_distinct_prior(k=i, n_j = n, gamma = 2, prior = "Poisson", lambda = 2, Max_iter = 1000)
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

# Test 3 groups
n1 = 2;n2 = 2;n3=0
n = c(n1,n2,n3)
tot = 0
for(i in 0:sum(n)){
  pp = p_distinct_prior(k=i, n_j = n, gamma = rep(2,length(n)), prior = "Poisson", lambda = 2, Max_iter = 1000)
  print(pp)
  tot = tot + pp
}
tot

# Test d groups
d = 4
n_max =20
n = sample(1:n_max,size = d,replace = T)
n
tot = 0
for(i in 0:sum(n)){
  pp = p_distinct_prior(k=i, n_j = n, gamma = rep(2,length(n)), prior = "Poisson", lambda = 2, Max_iter = 1000)
  print(pp)
  tot = tot + pp
}
tot




library(GDFMM)
p_distinct_prior(k=3, n_j = c(1,1,0), gamma = c(2,2,2), prior = "Poisson", lambda = 2, Max_iter = 1000)

