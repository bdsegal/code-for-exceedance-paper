library(ggplot2)
library(reshape2)
library(dplyr)

source("functions.R")
paper.path <- "../paper"
present.path <- "../presentation"

# simulation of lm slope ------------------------------------------------------

K <- 1e4  # number of iterations
n <- 100  # sample size for observed data
d <- 2
m.vec <- c(50, 100, 150)  # sample size for future data

beta <- c(1, 2)
nu <- 5  # sd of generating distribution
alpha <- 0.05

cutoff.vec <- seq(from = beta[2] - 1.5, to = beta[2] + 1.5, by = 0.2)

coverage.long <- NULL

# generate common covariate values x points for all simulations
x <- runif(n = n, min = 0, max = 10)
X <- cbind(1, x)
sd.beta <- sqrt(n * nu^2 * solve(t(X) %*% X)[2, 2])
  
for (m in m.vec) {

  all.data.list <- vector(mode = "list", length = K)

  for (k in 1:K) {

    if(k%%100 == 0) {print(k)}

    # generate data of size n
    y <- rnorm(n = n, mean = beta[1] + x * beta[2], sd = nu)
    fit <- lm(y ~ x)
    theta.hat <- coef(fit)[2]
    sd.hat <- sqrt(n * vcov(fit)[2, 2])

    delta.ci <- getDeltaCI(cutoff = cutoff.vec, 
                           theta.hat = theta.hat, 
                           sd.hat = sd.hat, 
                           n = n,
                           d = d, 
                           alpha = alpha)

    lower <- pnorm(sqrt(m/n) * delta.ci["upper", ], lower.tail = FALSE)
    upper <- pnorm(sqrt(m/n) * delta.ci["lower", ], lower.tail = FALSE)

    all.data.list[[k]] <- data.frame(cutoff = cutoff.vec,
                                     lower = lower, 
                                     upper = upper,
                                     k = k)
  }

  # get true exceedance probabilities for data of size m
  exceed.true <- exceedProb(cutoff.vec = cutoff.vec, theta.hat = beta[2], sd.hat = sd.beta, m = m)

  coverage.prob <- data.frame(cutoff = cutoff.vec)

  coverage.list <- lapply(all.data.list,
                          function(dat) {
                            dat$lower <= exceed.true & 
                            exceed.true <= dat$upper
                          })

  coverage <- do.call(cbind, coverage.list)
  coverage.prob$cp <- apply(coverage, 1, mean)
  coverage.prob$m <- paste0("m = ", prettyNum(m, big.mark = ","))

  coverage.long <- rbind(coverage.long, coverage.prob)
}

coverage.long$m <- factor(coverage.long$m, 
                          levels = paste0("m = ", prettyNum(m.vec, big.mark = ",")))

coverage.long <- coverage.long %>%
                   mutate(low = pmax(0, cp - 1.96 * sqrt(cp * (1 - cp) / K)),
                          high = pmin(1, cp + 1.96 * sqrt(cp * (1 - cp) / K)))

# save(coverage.long, file = "coverage_long_slr.Rdata")
load("coverage_long_slr.Rdata")

# plot coverage probability
dev.new(width = 9, height = 3.5)
ggplot(aes(x = cutoff, y = cp), data = coverage.long) + 
  geom_point() +
  # geom_errorbar(aes(ymin = low, ymax = high)) +
  facet_wrap(~ m) +
  theme_bw(17) + 
  labs(x = "Cutoff, c",
       y = "Coverage probability") +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(file.path(paper.path, "lm_cov_prob_100.png"))


coverage.long.present <- coverage.long
coverage.long.present$m <- as.character(coverage.long.present$m)
coverage.long.present$m <- gsub("m", "Future sample size", coverage.long.present$m)
coverage.long.present$m <- factor(coverage.long.present$m,
                                  levels = paste0("Future sample size = ", c(50, 100, 150)))

# plot coverage probability
dev.new(width = 9, height = 3.5)
ggplot(aes(x = cutoff, y = cp), data = coverage.long.present) + 
  # geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  theme_bw(18) + 
  facet_wrap(~ m) +
  labs(x = "Cutoff c",
       y = "Coverage probability") +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed", color = "black") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(file.path(present.path, "lm_cov_prob_100.png"))
