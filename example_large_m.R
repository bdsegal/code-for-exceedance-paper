library(ggplot2)
library(survival)
source("functions.R")
paper.path <- "/home/bsegal/Dropbox/Research/exceedance/exceedance_prob/paper"

# sample mean -----------------------------------------------------------------
set.seed(123)

n <- 100
d <- 1
m.vec <- c(100, 1000, 10000)

sd <- 1
x <- rnorm(n = n, mean = 0, sd = sd)

theta.hat <- mean(x)
sd.hat <- sd(x)

alpha <- 0.05
t.alpha <- qt(p = 1 - alpha / 2, df = n - d, lower.tail = TRUE)

cutoff.vec <- seq(from = theta.hat - 1*sd, to = theta.hat + 1*sd, by = 0.01)

exceed.point.list <- list()
exceed.ci.list <- list()

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
}

exceed.point <- do.call(rbind, exceed.point.list)
exceed.ci <- do.call(rbind, exceed.ci.list)

# standard confidence interval
stand.ci <- data.frame(lower = theta.hat - t.alpha * sd.hat / sqrt(n),
                       upper = theta.hat + t.alpha * sd.hat / sqrt(n),
                       point = theta.hat,
                       val = 0.5)

dev.new(width = 9, height = 3.5)
ggplot() +
  geom_polygon(aes(x = c, y = val), data = exceed.ci, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = c, y = val), data = exceed.point) + 
  geom_errorbarh(aes(x = point, xmin = lower, xmax = upper, y = val), data = stand.ci,
                 height = 0.1, linetype = "solid") +
  geom_point(aes(x = point, y = val), data = stand.ci) +
  theme_bw(14) +
  facet_wrap(~ m) +
  labs(x = "Cutoff c", 
       y = expression(paste("Pr",""[hat(theta)],","[hat(sigma)],"(", tilde(theta)^m, ">c)", sep = ""))) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(file.path(paper.path, "S_big_m.png"))
