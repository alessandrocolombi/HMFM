# Test one group
p_distinct_posterior(r=0,k=1, m_j = c(0,2), n_j = c(2,2), gamma = c(2,2), prior = "Poisson", lambda = 2, Max_iter = 100)
p_distinct_prior(k=1, n_j = 3, gamma = 2, prior = "Poisson", lambda = 2, Max_iter = 1000)

n = 250
tot = 0
for(i in 0:n){
  pp = p_distinct_prior(k=i, n_j = n, gamma = 2, prior = "Poisson", lambda = 2, Max_iter = 1000)
  print(pp)
  tot = tot + pp
}
tot
