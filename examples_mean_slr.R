library(ggplot2)
library(BayesFactor)
library(dplyr)
source("functions.R")
paper.path <- "../paper"
present.path <- "../presentation"

alpha <- 0.05

# sample mean -----------------------------------------------------------------
set.seed(12345)

n <- 100
d <- 1
m.vec <- c(50, 100, 150)

sd <- 1

x <- rnorm(n = n, mean = 0, sd = sd)
theta.hat <- mean(x)
sd.hat <- sd(x)

c(theta.hat, sd.hat)
# [1] 0.2451972 1.1147308

cutoff.vec <- seq(from = theta.hat - 1*sd, to = theta.hat + 1*sd, by = 0.01)

exceed.point.list <- list()
exceed.ci.list <- list()
lower.zero <- data.frame(m = m.vec, lower = NA)

for (j in 1:length(m.vec)) {

  m <- m.vec[j]

  point.est <- exceedProb(cutoff.vec, 
                          theta.hat = theta.hat, 
                          sd.hat = sd.hat, 
                          m = m)

  delta.ci <- getDeltaCI(cutoff = cutoff.vec, 
                         theta.hat = theta.hat, 
                         sd.hat = sd.hat, 
                         n = n, 
                         d = d,
                         alpha = alpha)

  lower <- pnorm(sqrt(m/n) * delta.ci["upper", ], lower.tail = FALSE)
  upper <- pnorm(sqrt(m/n) * delta.ci["lower", ], lower.tail = FALSE)

  # point estimates for exceedance prob
  exceed.point.list[[j]] <- data.frame(c = cutoff.vec, 
                                       val = point.est, 
                                       m = paste0("m = ", prettyNum(m, big.mark = ",")))

  # confidence intervals for exceedance prob
  exceed.ci.list[[j]] <- data.frame(c = c(cutoff.vec, rev(cutoff.vec)), 
                                    val = c(upper, rev(lower)),
                                    m = paste0("m = ", prettyNum(m, big.mark = ",")))

  delta.ci.zero <- getDeltaCI(cutoff = 0, 
                         theta.hat = theta.hat, 
                         sd.hat = sd.hat, 
                         n = n, 
                         d = d,
                         alpha = alpha)

  lower.zero$lower[j] <- pnorm(sqrt(m/n) * delta.ci.zero["upper", ], lower.tail = FALSE)
}

exceed.point <- do.call(rbind, exceed.point.list)
exceed.ci <- do.call(rbind, exceed.ci.list)

dev.new(width = 5, height = 4)
ggplot(aes(x = x), data = data.frame(x = x)) +
  geom_histogram(binwidth = 0.5, boundary = 0, fill = "grey", color = "black") +
  theme_bw(17) +
  labs(y = "Count", x = "y") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave(file.path(paper.path, "sample_mean_100.png"))

dev.new(width = 5, height = 4)
ggplot(aes(x = x), data = data.frame(x = x)) +
  geom_histogram(binwidth = 0.5, boundary = 0, fill = "grey", color = "black") +
  theme_bw(16) +
  labs(y = "Count", x = "y",
       title = bquote(paste(bar(y), " = ", .(signif(theta.hat, 2)), ", ",
                            "sd = ", .(signif(sd.hat, 2)), ", ",
                            "n = ", .(n), sep = "")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave(file.path(present.path, "sample_mean_100.png"))

dev.new(width = 9, height = 3.5)
ggplot() +
  geom_polygon(aes(x = c, y = val), data = exceed.ci, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = c, y = val), data = exceed.point) + 
  theme_bw(17) +
  facet_wrap(~ m) +
  labs(x = "Cutoff c", 
       y = expression(paste("Pr", ""[hat(theta)], ","[hat(sigma)], "(", 
                            hat(theta)^"rep", ">c)", sep = ""))) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave(file.path(paper.path, "S_100.png"))

signif(lower.zero, 2)
#     m lower
# 1  50  0.56
# 2 100  0.58
# 3 150  0.60

exceed.ci.present <- exceed.ci
exceed.ci.present$m <- gsub("m", "Future sample size", exceed.ci.present$m)
exceed.ci.present$m <- factor(exceed.ci.present$m,
                              levels = paste0("Future sample size = ", c(50, 100, 150)))

exceed.point.present <- exceed.point
exceed.point.present$m <- gsub("m", "Future sample size", exceed.point.present$m)
exceed.point.present$m <- factor(exceed.point.present$m,
                              levels = paste0("Future sample size = ", c(50, 100, 150)))

dev.new(width = 9, height = 3.5)
ggplot() +
  geom_polygon(aes(x = c, y = val), data = exceed.ci.present, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = c, y = val), data = exceed.point.present) + 
  theme_bw(18) +
  facet_wrap(~ m) +
  labs(x = "Cutoff c", 
       y = expression(paste("Pr", ""[hat(theta)], ","[hat(sigma)], "(", 
                            hat(theta)^"new", ">c)", sep = ""))) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave(file.path(present.path, "S_100.png"))


# t-tests and Bayes factors

t.test(x, alternative = "greater")
#   One Sample t-test
# data:  x
# t = 2.1996, df = 99, p-value = 0.01508
# alternative hypothesis: true mean is greater than 0
# 95 percent confidence interval:
#  0.06010828        Inf
# sample estimates:
# mean of x 
# 0.2451972 

bf.int <- ttestBF(x, nullInterval = c(-Inf, 0))
bf <- bf.int[1] / bf.int[2]
bf
# Bayes factor analysis
# --------------
# [1] Alt., r=0.707 -Inf<d<0 : 0.01648142 ±0.04%

# Against denominator:
#   Alternative, r = 0.707106781186548, mu =/= 0 !(-Inf<d<0) 
# ---
# Bayes factor type: BFoneSample, JZS

1/exp(bf@bayesFactor$bf)
# [1] 60.67438

# two-sided interval
theta.hat + c(-1, 1) * qt(p = 1 - alpha/2, df = n - d) * sd.hat / sqrt(n)
# [1] 0.02401042 0.46638398

t.test(x)

bf.zero <- ttestBF(x, mu = 0)
1 / exp(bf.zero@bayesFactor$bf)

# linear regression -----------------------------------------------------------
set.seed(12345)

n <- 100
d <- 2
m.vec <- c(50, 100, 150)

beta <- c(1, 2)
nu <- 5

x <- runif(n = n, min = 0, max = 10)
y <- rnorm(n = n, mean = beta[1] + x * beta[2], sd = nu)

dev.new(width = 5, height = 4)
qplot(x = x, y = y) +
  theme_bw(17) +
  geom_abline(intercept = beta[1], slope = beta[2], linetype = "dashed")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(file.path(paper.path, "lin_reg_data_100.png"))

fit <- lm(y ~ x)
theta.hat <- coef(fit)[2]
Sigma <- vcov(fit)
sd.hat <- sqrt(Sigma[2, 2]) * sqrt(n)

c(theta.hat, sd.hat)
# 1.856129 1.906534

cutoff.vec <- seq(from = theta.hat - 2, to = theta.hat + 2, by = 0.01)

exceed.point.list <- list()
exceed.ci.list <- list()

for (j in 1:length(m.vec)) {
  m <- m.vec[j]

  point.est <- exceedProb(cutoff.vec, 
                          theta.hat = theta.hat, 
                          sd.hat = sd.hat, 
                          m = m)

  delta.ci <- getDeltaCI(cutoff = cutoff.vec, 
                         theta.hat = theta.hat, 
                         sd.hat = sd.hat, 
                         n = n, 
                         d = d,
                         alpha = alpha)

  lower <- pnorm(sqrt(m/n) * delta.ci["upper", ], lower.tail = FALSE)
  upper <- pnorm(sqrt(m/n) * delta.ci["lower", ], lower.tail = FALSE)

  # point estimates for exceedance prob
  exceed.point.list[[j]] <- data.frame(c = cutoff.vec, 
                                       val = point.est, 
                                       m = paste0("m = ", prettyNum(m, big.mark = ",")))

  # confidence intervals for exceedance prob
  exceed.ci.list[[j]] <- data.frame(c = c(cutoff.vec, rev(cutoff.vec)), 
                                    val = c(upper, rev(lower)),
                                    m = paste0("m = ", prettyNum(m, big.mark = ",")))
}

exceed.point <- do.call(rbind, exceed.point.list)
exceed.ci <- do.call(rbind, exceed.ci.list)

dev.new(width = 9, height = 3.5)
ggplot(aes(x = c, y = val), data = exceed.point) +
  geom_polygon(data = exceed.ci, fill = "grey", color = "grey", linetype = "solid") +
  geom_line() + 
  theme_bw(17) +
  facet_wrap(~ m) +
  labs(x = "Cutoff c", 
       y = expression(paste("Pr",""[hat(theta)],","[hat(sigma)],"(", tilde(theta)^m, ">c)", sep = ""))) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 1.5, linetype = "dashed")
ggsave(file.path(paper.path, "S_lin_reg_100.png"))

# t-tests and confidence intervals

# one-sided
t <- sqrt(n) * (theta.hat - 1.5) / sd.hat
pt(t, lower.tail = FALSE, df = n - d)
# 0.0323793 

# lower bound of one-sided confidence interval
theta.hat - qt(p = 1-alpha, df = n - d) * sd.hat / sqrt(n)
# 1.539539 

# two-sided confidence interval
theta.hat + c(-1, 1) * qt(p = 1-alpha/2, df = n - d) * sd.hat / sqrt(n)
# [1] 1.477783 2.234474
