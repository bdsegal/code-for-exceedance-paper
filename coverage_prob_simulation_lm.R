library(exceedProb)
library(tidyr)

# simulation of lm slope ------------------------------------------------------
set.seed(12)

K <- 1e4 # number of iterations
n_vec <- c(20, 40, 60, 80, 100)  # sample size for observed data
d <- 2

beta <- c(1, 2)
nu <- 5  # sd of generating distribution
alpha <- 0.05

cutoff <- seq(from = beta[2] - 1, to = beta[2] + 1, by = 0.1)

coverage_long <- NULL

for (n in n_vec) {

  m <- n

  # generate common covariate values x points for all simulations
  x <- runif(n = n, min = 0, max = 10)
  X <- cbind(1, x)
  
  sd_true <- sqrt(n * nu^2 * solve(t(X) %*% X)[2, 2])

  ep_list <- vector(mode = "list", length = K)

  for (k in 1:K) {


    if(k%%100 == 0) {print(paste0("n = ", n, ", k = ", k))}

    # generate data of size n
    y <- rnorm(n = n, mean = beta[1] + x * beta[2], sd = nu)
    fit <- lm(y ~ x)
    theta_hat <- coef(fit)[2]
    sd_hat <- sqrt(n * vcov(fit)[2, 2])

    
    ep_list[[k]] <- exceedProb(cutoff = cutoff, 
          theta_hat = theta_hat, 
          sd_hat = sd_hat, 
          d = d,
          alpha = alpha, 
          n = n,
          m = n)

    ep_list[[k]]$k <- k

  }

  # get true exceedance probabilities for data of size m
  exceed_true <- exceedProb(cutoff = cutoff, 
                            theta_hat = beta[2], 
                            sd_hat = sd_true, 
                            d = d,
                            alpha = alpha, 
                            n = n,
                            m = n)$point

  coverage_prob <- data.frame(cutoff = cutoff)

  coverage_list <- lapply(ep_list,
                          function(ep) {
                            ep$lower <= exceed_true & 
                            exceed_true <= ep$upper
                          })

  coverage <- do.call(cbind, coverage_list)
  coverage_prob$cp <- apply(coverage, 1, mean)
  coverage_prob$n <- n

  coverage_long <- rbind(coverage_long, coverage_prob)
}

coverage_long %>%
  spread(n, cp) %>%
  signif(., 3)
