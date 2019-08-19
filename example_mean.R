library(exceedProb)
library(ggplot2)
library(BayesFactor)
library(dplyr)

paper.path <- "../paper"

alpha <- 0.05

# sample mean -----------------------------------------------------------------
set.seed(12345)

n <- 100
d <- 1

sd_true <- 1

x <- rnorm(n = n, mean = 0, sd = sd_true)
theta_hat <- mean(x)
sd_hat <- sd(x)

c(theta_hat, sd_hat)
# [1] 0.2451972 1.1147308

cutoff <- seq(from = theta_hat - 1*sd_true, to = theta_hat + 1*sd_true, by = 0.01)

ep <- exceedProb(cutoff = cutoff, 
      theta_hat = theta_hat, 
      sd_hat = sd_hat, 
      d = d,
      alpha = alpha, 
      n = n,
      m = n)

# confidence intervals for exceedance prob
exceed_ci <- with(ep, data.frame(cutoff = c(cutoff, rev(cutoff)), 
                        val = c(upper, rev(lower))))

ep0 <- exceedProb(cutoff = 0, 
      theta_hat = theta_hat, 
      sd_hat = sd_hat, 
      d = d,
      alpha = alpha, 
      n = n,
      m = n)

dev.new(width = 5, height = 4)
ggplot(aes(x = x), data = data.frame(x = x)) +
  geom_histogram(binwidth = 0.5, boundary = 0, fill = "grey", color = "black") +
  theme_bw(17) +
  labs(y = "Count", x = "y") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave(file.path(paper.path, "sample_mean_100.png"))

dev.new(width = 5.4, height = 3.5)
ggplot() +
  geom_polygon(aes(x = cutoff, y = val), data = exceed_ci, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = cutoff, y = point), data = ep) + 
  theme_bw(17) +
  labs(x = "Cutoff c", 
       y = expression(paste("Pr", ""[hat(theta)], ","[hat(sigma)], "(", 
                            hat(theta)^"rep", ">c)", sep = ""))) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave(file.path(paper.path, "S_100.png"))

signif(ep0, 2)
#   cutoff point lower upper
# 1      0  0.99  0.58     1

# t-tests and Bayes factors ---------------------------------------------------

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

bf_int <- ttestBF(x, nullInterval = c(-Inf, 0))
bf <- bf_int[1] / bf_int[2]
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
t.test(x)
#   One Sample t-test

# data:  x
# t = 2.1996, df = 99, p-value = 0.03016
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#  0.02401042 0.46638398
# sample estimates:
# mean of x 
# 0.2451972 

bf_zero <- ttestBF(x, mu = 0)
1 / exp(bf_zero@bayesFactor$bf)
# [1] 0.9056418
