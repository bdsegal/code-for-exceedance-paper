library(exceedProb)
library(tidyr)

# simulation of sample mean ---------------------------------------------------
set.seed(12)

K <- 1e4  # number of iterations
n_vec <- c(20, 40, 60, 80, 100)  # sample size for observed data
d <- 1

mu <- 0  # mean of generating distribution
sd_true <- 1  # sd of generating distribution
alpha <- 0.05

cutoff <- seq(from = -0.5, to = 0.5, by = 0.1)

coverage_long <- NULL

for (n in n_vec) {

  m <- n

  ep_list <- vector(mode = "list", length = K)

  for (k in 1:K) {

    if(k%%100 == 0) {print(paste0("n = ", n, ", k = ", k))}
    
    # simulate data of size n
    x <- rnorm(n = n, mean = mu, sd = sd_true)
    theta_hat <- mean(x)
    sd_hat <- sd(x)

    ep_list[[k]] <- exceedProb(cutoff = cutoff, 
          theta_hat = theta_hat, 
          sd_hat = sd_hat, 
          d = d,
          alpha = alpha, 
          n = n,
          m = n)

    ep_list[[k]]$k <- k

  }

  # get true exceedance probability for data of size m
  exceed_true <- exceedProb(cutoff = cutoff, 
                            theta_hat = mu, 
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
