library(ggplot2)
library(dplyr)
library(reshape2)
source("functions.R")
paper.path <- "../paper"
present.path <- "../exceedance_prob/presentation"

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

ci <- theta.hat + c(-1, 1) * t.alpha * sd.hat / sqrt(n)

dev.new(width = 9, height = 3.5)
ggplot() +
  geom_polygon(aes(x = c, y = val), data = exceed.ci, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = c, y = val), data = exceed.point) + 
  geom_errorbarh(aes(x = point, xmin = lower, xmax = upper, y = val), data = stand.ci,
                 height = 0.1, linetype = "solid") +
  geom_vline(xintercept = ci, linetype = "dashed") +
  geom_point(aes(x = point, y = val), data = stand.ci) +
  theme_bw(17) +
  facet_wrap(~ m) +
  labs(x = "Cutoff c", 
       y = expression(paste("Pr",""[hat(theta)],","[hat(sigma)],"(", hat(theta)^"rep", ">c)", sep = ""))) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
 scale_x_continuous(breaks = c(-1, ci, 1),
                    labels = c(-1, 
                               expression(theta[L]),
                               expression(theta[U]),
                               1))
ggsave(file.path(paper.path, "S_big_m.png"))

exceed.ci.present <- exceed.ci
exceed.ci.present$m <- gsub("m", "Future sample size", exceed.ci.present$m)
exceed.ci.present$m <- factor(exceed.ci.present$m,
                              levels = paste0("Future sample size = ", c(100, "1,000","10,000")))

exceed.point.present <- exceed.point
exceed.point.present$m <- gsub("m", "Future sample size", exceed.point.present$m)
exceed.point.present$m <- factor(exceed.point.present$m,
                              levels = paste0("Future sample size = ", c(100, "1,000", "10,000")))

dev.new(width = 9, height = 3.5)
ggplot() +
  geom_polygon(aes(x = c, y = val), data = exceed.ci.present, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = c, y = val), data = exceed.point.present) + 
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = val), data = stand.ci,
                 height = 0.1, linetype = "solid") +
  geom_point(aes(x = point, y = val), data = stand.ci) +
  theme_bw(16) +
  facet_wrap(~ m) +
  labs(x = "Cutoff c", 
       y = expression(paste("Pr",""[hat(theta)],","[hat(sigma)],"(", hat(theta)^"new", ">c)", sep = ""))) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(file.path(present.path, "S_big_m.png"))


exceed.ci.present.100 <- exceed.ci.present[grep(100, exceed.ci.present$m), ]
exceed.point.present.100 <- exceed.point.present[grep(100, exceed.point.present$m), ]

dev.new(width = 3.5, height = 3.5)
ggplot() +
  geom_polygon(aes(x = c, y = val), data = exceed.ci.present.100, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = c, y = val), data = exceed.point.present.100) + 
  geom_errorbarh(aes(x = point, xmin = lower, xmax = upper, y = val), data = stand.ci,
                 height = 0.1, linetype = "solid") +
  geom_point(aes(x = point, y = val), data = stand.ci) +
  theme_bw(16) +
  facet_wrap(~m) +
  labs(x = "Cutoff c", 
       y = expression(paste("Pr",""[hat(theta)],","[hat(sigma)],"(", hat(theta)^"new", ">c)", sep = ""))) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(file.path(present.path, "S_big_m_100.png"))


# comparison to p-value -------------------------------------------------------
m <- 100
n <- 100

delta.ci <- getDeltaCI(cutoff = cutoff.vec, 
                       theta.hat = theta.hat, 
                       sd.hat = sd.hat, 
                       n = n, 
                       d = d,
                       alpha = alpha)

point.est <- exceedProb(cutoff.vec, theta.hat = theta.hat, sd.hat = sd.hat, m = m)
lower <- pnorm(sqrt(m/n) * delta.ci["upper", ], lower.tail = FALSE)
upper <- pnorm(sqrt(m/n) * delta.ci["lower", ], lower.tail = FALSE)

p.two.sided <- 2 * (1 - pt(sqrt(n) * abs(theta.hat - cutoff.vec) / sd.hat, df = n - 1))

exceed.p <- data.frame(cutoff.vec, pr.lt = 1-lower, p.two.sided)

exceed.ci <- data.frame(c = c(cutoff.vec, rev(cutoff.vec)), 
                        val = c(upper, rev(lower)))

exceed.p.val <- data.frame(cutoff.vec,
                           exceed = 1 - point.est,
                           p.two.sided = p.two.sided)
exceed.p.val.m <- melt(exceed.p.val, id = "cutoff.vec")

dev.new(width = 7, height = 4.5)
ggplot() +
  geom_polygon(aes(x = c, y = 1 - val), data = exceed.ci, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = cutoff.vec, y = value, color = variable),
            data = exceed.p.val.m) +
  theme_bw(17) +
  labs(x = "Cutoff c", 
       y = "") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = ci, linetype = "dashed") +
  geom_hline(yintercept = c(0.05, 0.5), linetype = "dashed") +
  scale_y_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = c(-1, ci, 1),
                    labels = c(-1, 
                               expression(theta[L]),
                               expression(theta[U]),
                               1)) +
  scale_color_manual("", values = c("grey", "black"),
                     labels = c(expression(paste("Pr",""[hat(theta)],","[hat(sigma)],"(", hat(theta)^"rep" <="c)", sep = "")),
                                "p-value")) +
  guides(color = guide_legend(override.aes = list(size = c(5, 0.6))))
ggsave(file.path(paper.path, "ci_p_val_overlay.png"))
