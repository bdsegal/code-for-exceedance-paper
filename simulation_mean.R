library(ggplot2)
library(reshape2)
library(dplyr)

source("functions.R")
paper.path <- "/home/bsegal/Dropbox/Research/exceedance/exceedance_prob/paper/segal_exceedance_TAS"
present.path <- "/home/bsegal/Dropbox/Research/exceedance/exceedance_prob/presentation"

# simulation of sample mean ---------------------------------------------------

K <- 1e4  # number of iterations
n <- 100  # sample size for observed data
d <- 1
m.vec <- c(50, 100, 150)  # sample size for future data

mu <- 0  # mean of generating distribution
sd <- 1  # sd of generating distribution
alpha <- 0.05

cutoff.vec <- seq(from = mu - sd, to = mu + sd, by = 0.1)

coverage.long <- NULL

for (m in m.vec) {

  all.data.list <- vector(mode = "list", length = K)

  for (k in 1:K) {

    if(k%%100 == 0) {print(k)}
    
    # simulate data of size n
    x <- rnorm(n = n, mean = mu, sd = sd)
    theta.hat <- mean(x)
    sd.hat <- sd(x)

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

  # get true exceedance probability for data of size m
  exceed.true <- exceedProb(cutoff.vec = cutoff.vec, 
                            theta.hat = mu, 
                            sd.hat = sd, 
                            m = m)

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

# save(coverage.long, file = "coverage_long_mean.Rdata")
load("coverage_long_mean.Rdata")

coverage.long$m <- factor(coverage.long$m, 
                          levels = paste0("m = ", prettyNum(m.vec, big.mark = ",")))

# adjusted Wald confidence intervals
coverage.long <- coverage.long %>%
                   mutate(cp.tilde = (cp * K + 2) / (K + 4),
                          low = pmax(0, cp.tilde - 1.96 * sqrt(cp.tilde * (1 - cp.tilde) / (K + 4))),
                          high = pmin(1, cp.tilde + 1.96 * sqrt(cp.tilde * (1 - cp.tilde) / (K + 4))))


# plot coverage probability
dev.new(width = 9, height = 3.5)
ggplot(aes(x = cutoff, y = cp), data = coverage.long) + 
  # geom_line() +
  geom_point() +
  # geom_errorbar(aes(ymin = low, ymax = high)) +
  theme_bw(17) + 
  facet_wrap(~ m) +
  labs(x = "Cutoff c",
       y = "Coverage probability") +
  geom_hline(yintercept = 1 - alpha, linetype = "dashed", color = "black") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(file.path(paper.path, "sample_mean_cov_prob_100.png"))

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
ggsave(file.path(present.path, "sample_mean_cov_prob_100.png"))
