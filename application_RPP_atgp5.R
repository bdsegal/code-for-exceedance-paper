# Heading from R script at https://osf.io/9ivaj -------------------------------

##@@ GENERAL INFORMATION @@##
#@ A single-system account of the relationship between priming, recognition, and fluency.
#@ Coder name: John Hodsoll
#@ Coder e-mail: john.hodsoll@kcl.ac.uk
#@ Type of statistic: t test
#@ Type of effect-size: Cohen's D
#@ OSF link project: https://osf.io/atgp5/
#@ OSF link replication report: https://osf.io/yc2fe/
#
##@@ REQUIRED PACKAGES @@##
library("httr")
library("RCurl")
source("http://sachaepskamp.com/files/OSF/getOSFfile.R") # the getOSFfile function
  
##@@ DATA LOADING @@##
#@ NOTE: Data must be loaded from OSF directly and NOT rely on any local files.

#Option to download SPSS data file to tab delimited text (.dat). 
#Download dat data file as no need for file conversion for reading into R
file <- getOSFfile("https://osf.io/btr65/")
berry.data <- read.table(file, sep="\t", header=TRUE) 

# Newly added code ------------------------------------------------------------
library(ggplot2)
library(exceedProb)
library(BayesFactor)
paper.path <- "../paper"

alpha <- 0.05

x <- with(berry.data, meanRT_cr - meanRT_miss)
n <- length(x)
m.vec <- n
d <- 1

dev.new(width = 5, height = 4)
ggplot(aes(x = x), data = data.frame(x = x)) +
  geom_histogram(binwidth = 50, boundary = 0, fill = "grey", color = "black") +
  theme_bw(17) +
  labs(y = "Count", x = "Within-volunteer difference in RT (ms)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave(file.path(paper.path, "application_data_hist.png"))

# exceedance probability ----------------------------------
theta_hat <- mean(x)
sd_hat <- sd(x)

alpha <- 0.05
t_alpha <- qt(p = 1 - alpha / 2, df = n - d, lower.tail = TRUE)

cutoff <- seq(from = theta_hat - sd_hat, to = theta_hat + sd_hat, by = 1)

ep <- exceedProb(cutoff = cutoff, 
      theta_hat = theta_hat, 
      sd_hat = sd_hat, 
      d = d,
      alpha = alpha, 
      n = n,
      m = n)

exceed_ci <- with(ep, data.frame(cutoff = c(cutoff, rev(cutoff)), 
                        val = c(upper, rev(lower))))

ep0 <- exceedProb(cutoff = 0, 
      theta_hat = theta_hat, 
      sd_hat = sd_hat, 
      d = d,
      alpha = alpha, 
      n = n,
      m = n)

signif(ep0, 3)
#   cutoff point lower upper
# 1      0 0.992  0.63     1

dev.new(width = 5.4, height = 3.5)
ggplot() +
  geom_polygon(aes(x = cutoff, y = val), data = exceed_ci, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = cutoff, y = point), data = ep) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(17) +
  # facet_wrap(~ m) +
  labs(x = "Within-volunteer difference in RT (c in ms)", 
       y = expression(paste("Pr",""[hat(theta)],","[hat(sigma)],"(", hat(theta)^"rep", "> c)", sep = ""))) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(file.path(paper.path, "application_exceedance.png"))

wilcox.test(x)
#   Wilcoxon signed rank test

# data:  x
# V = 386, p-value = 0.02162
# alternative hypothesis: true location is not equal to 0

t.test(x)
#   One Sample t-test

# data:  x
# t = 2.3977, df = 31, p-value = 0.02271
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#    8.647711 107.128454
# sample estimates:
# mean of x 
#  57.88808

t.test(x, alternative = "greater")
#   One Sample t-test

# data:  x
# t = 2.3977, df = 31, p-value = 0.01135
# alternative hypothesis: true mean is greater than 0
# 95 percent confidence interval:
#  16.95284      Inf
# sample estimates:
# mean of x 
#  57.88808 

# Bayes factor
bf_int <- ttestBF(x, nullInterval = c(-Inf, 0))
bf <- bf_int[1] / bf_int[2]
bf
# Bayes factor analysis
# --------------
# [1] Alt., r=0.707 -Inf<d<0 : 0.01410145 ±0.1%

# Against denominator:
#   Alternative, r = 0.707106781186548, mu =/= 0 !(-Inf<d<0) 
# ---
# Bayes factor type: BFoneSample, JZS

1/exp(bf@bayesFactor$bf)
# [1] 70.91471

bf_zero <- ttestBF(x, mu = 0)
1 / exp(bf_zero@bayesFactor$bf)
# [1] 0.4512185
