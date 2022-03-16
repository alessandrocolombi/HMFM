
# Test one group
n = 25
tot = 0
for(i in 0:n){
  pp = p_distinct_prior(k=i, n_j = n, gamma = 2, prior = "Poisson", lambda = 2, Max_iter = 1000)
  print(pp)
  tot = tot + pp
}
tot


# Test 2 groups
n1 = 7;n2 = 9
n = c(n1,n2)
tot = 0
for(i in 0:sum(n)){
  pp = p_distinct_prior(k=i, n_j = n, gamma = c(2,2), prior = "Poisson", lambda = 2, Max_iter = 1000)
  print(pp)
  tot = tot + pp
}
tot


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



# Test Shared species

# d=2
n1 = 50;n2 = 30
n = c(n1,n2)
tot = 0
for(i in 0:min(n)){
  pp = p_shared_prior(s=i, n_j = n, gamma = c(2,2), prior = "Poisson", lambda = 2, Max_iter = 1000)
  print(pp)
  tot = tot + pp
}
tot

# d=3
n1 = 5;n2 = 7;n3=10
n = c(n1,n2,n3)
tot = 0
for(i in 0:min(n)){
  pp = p_shared_prior(s=i, n_j = n, gamma = rep(2,length(n)), prior = "Poisson", lambda = 2, Max_iter = 100)
  print(pp)
  tot = tot + pp
}
tot

# Test d groups
d = 4
n_max =10
n = sample(1:n_max,size = d,replace = T)
n
tot = 0
for(i in 0:min(n)){
  pp = p_shared_prior(s=i, n_j = n, gamma = rep(2,length(n)), prior = "Poisson", lambda = 2, Max_iter = 100)
  print(pp)
  tot = tot + pp
}
tot







library(GDFMM)
p_distinct_prior(k=1, n_j = c(0,0,1), gamma = c(2,2,2), prior = "Poisson", lambda = 2, Max_iter = 1000)
p_distinct_prior(k=10, n_j = c(7,7), gamma = c(2,2), prior = "Poisson", lambda = 2, Max_iter = 100)
p_distinct_prior(k=0, n_j = 0, gamma = 2, prior = "Poisson", lambda = 2, Max_iter = 10)


library(GDFMM)
p_shared_prior(s=0,n_j = c(2,2,2), gamma = c(2,2,2), prior = "Poisson", lambda = 2, Max_iter = 1000)


0.4848428 + 0.4662059 + 0.0248566
