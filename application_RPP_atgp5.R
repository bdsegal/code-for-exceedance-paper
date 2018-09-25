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
# save(berry.data, file = "berry_data.Rdata")
# load("berry_data.Rdata")

##@@ DATA MANIPULATION @@##
#@ NOTE: Include here ALL difference between OSF data and data used in analysis
#@ TIP: You will want to learn all about dplyr for manipulating data.

##@@ DATA ANLAYSIS @@##
#@ NOTE: Include a print or sumarry call on the resulting object

# exceedance probability ------------------------------------------------------
library(ggplot2)
library(BayesFactor)
# library(dplyr)
source("functions.R")
paper.path <- "../paper"
present.path <- "../presentation"

alpha <- 0.05

x <- with(berry.data, meanRT_cr - meanRT_miss)
n <- length(x)
m.vec <- c(25, 32, 50)
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

# qq-plots
qqnorm(x)
qqline(x)

theta.hat <- mean(x)
sd.hat <- sd(x)

alpha <- 0.05
t.alpha <- qt(p = 1 - alpha / 2, df = n - d, lower.tail = TRUE)

cutoff.vec <- seq(from = theta.hat - sd.hat, to = theta.hat + sd.hat, by = 10)

exceed.point.list <- list()
exceed.ci.list <- list()
lower.zero <- data.frame(m = m.vec, lower = NA)

for (j in 1:length(m.vec)) {
 m <- m.vec[j]

  point.est <- exceedProb(cutoff.vec, theta.hat = theta.hat, sd.hat = sd.hat, m = m)

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

# standard confidence interval
stand.ci <- data.frame(lower = theta.hat - t.alpha * sd.hat / sqrt(n),
                       upper = theta.hat + t.alpha * sd.hat / sqrt(n),
                       point = theta.hat,
                       val = 0.5)

ci <- theta.hat + c(-1, 1) * t.alpha * sd.hat / sqrt(n)
signif(ci, 3)

signif(lower.zero, 2)
#    m lower
# 1 25  0.62
# 2 32  0.63
# 3 50  0.66

dev.new(width = 9, height = 3.5)
ggplot() +
  geom_polygon(aes(x = c, y = val), data = exceed.ci, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = c, y = val), data = exceed.point) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(17) +
  facet_wrap(~ m) +
  labs(x = "Within-volunteer difference in RT (cutoff c in ms)", 
       y = expression(paste("Pr",""[hat(theta)],","[hat(sigma)],"(", hat(theta)^"rep", ">c)", sep = ""))) +
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
bf.int <- ttestBF(x, nullInterval = c(-Inf, 0))
bf <- bf.int[1] / bf.int[2]
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

bf.zero <- ttestBF(x, mu = 0)
1 / exp(bf.zero@bayesFactor$bf)
# [1] 0.4512185
