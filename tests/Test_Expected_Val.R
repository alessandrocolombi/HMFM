# Test Expected values

# Prior distinct
n1 = 5;n2 = 3
n = c(n1,n2)
Expected_prior(n_j = n, gamma = rep(2, length(n)), type = "distinct", prior = "Poisson",
               lambda = 2, Max_iter = 100, tol = 10^-12)

# Prior shared
n1 = 5;n2 = 3
n = c(n1,n2)
Expected_prior(n_j = n, gamma = rep(2, length(n)), type = "shared", prior = "Poisson",
               lambda = 2, Max_iter = 100, tol = 10^-12)

# Posterior distinct
m1 = 5;m2 = 3
n1 = 10;n2 = 10
n = c(n1,n2)
m = c(m1,m2)
k = 18
Expected_posterior(k = k, m_j = m, n_j = n, gamma = rep(2, length(n)), type = "distinct", prior = "Poisson",
                   lambda = 2, Max_iter = 100, tol = 10^-12)


# Posterior shared
m1 = 15;m2 = 13
n1 = 10;n2 = 10
n = c(n1,n2)
m = c(m1,m2)
k = 2
Expected_posterior(k = k, m_j = m, n_j = n, gamma = rep(2, length(n)), type = "shared", prior = "Poisson",
                   lambda = 2, Max_iter = 100, tol = 10^-12)















