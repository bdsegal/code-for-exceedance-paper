library(exceedProb)
library(ggplot2)
library(dplyr)
library(reshape2)

paper_path <- "../paper"

# sample mean -----------------------------------------------------------------
set.seed(123)

n <- 100
d <- 1
m_vec <- c(100, 1000, 10000)

sd_true <- 1
x <- rnorm(n = n, mean = 0, sd = sd_true)

theta_hat <- mean(x)
sd_hat <- sd(x)

alpha <- 0.05
t.alpha <- qt(p = 1 - alpha / 2, df = n - d, lower.tail = TRUE)

ci <- theta_hat + c(-1, 1) * t.alpha * sd_hat / sqrt(n)

cutoff_vec <- seq(from = theta_hat - 1*sd_true, to = theta_hat + 1*sd_true, by = 0.01)

exceed_point_list <- list()
exceed_ci_list <- list()
ci_thetaL <- data.frame(m = m_vec, cutoff = NA, point = NA, lower = NA, upper = NA)

for (j in 1:length(m_vec)) {

  m <- m_vec[j]

  ep <- exceedProb(cutoff = cutoff_vec, 
                   theta_hat = theta_hat, 
                   sd_hat = sd_hat, 
                   alpha = alpha, 
                   d = d,
                   n = n,
                   m = m)

  # point estimates for exceedance prob
  exceed_point_list[[j]] <- data.frame(c = ep$cutoff, 
                                       val = ep$point, 
                                       m = paste0("m = ", prettyNum(m, big.mark = ",")))

  # confidence intervals for exceedance prob
  exceed_ci_list[[j]] <- data.frame(c = c(cutoff_vec, rev(cutoff_vec)), 
                                    val = c(ep$upper, rev(ep$lower)),
                                    m = paste0("m = ", prettyNum(m, big.mark = ",")))

  ci_thetaL[j, ] <- exceedProb(cutoff = theta_hat - t.alpha * sd_hat / sqrt(n),
                               theta_hat = theta_hat, 
                               sd_hat = sd_hat, 
                               alpha = alpha, 
                               d = d,
                               n = n,
                               m = m)
}

ci_thetaL

exceed_point <- do.call(rbind, exceed_point_list)
exceed_ci <- do.call(rbind, exceed_ci_list)

# standard confidence interval
stand_ci <- data.frame(lower = theta_hat - t.alpha * sd_hat / sqrt(n),
                       upper = theta_hat + t.alpha * sd_hat / sqrt(n),
                       point = theta_hat,
                       val = 0.5)

dev.new(width = 9, height = 3.5)
ggplot() +
  geom_polygon(aes(x = c, y = val), data = exceed_ci, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = c, y = val), data = exceed_point) + 
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = val), data = stand_ci,
                 height = 0.1, linetype = "solid") +
  geom_vline(xintercept = ci, linetype = "dashed") +
  geom_point(aes(x = point, y = val), data = stand_ci) +
  theme_bw(17) +
  facet_wrap(~ m) +
  labs(x = "Cutoff c", 
       y = expression(paste("Pr",""[hat(theta)],","[hat(sigma)],"(", hat(theta)^"rep", ">c)", sep = ""))) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(-1, ci, 1),
                     labels = c(-1, expression(hat(theta)[L]), expression(hat(theta)[U]), 1))
ggsave(file.path(paper_path, "S_big_m.png"))

# comparison to p-value -------------------------------------------------------

# two-sided -----------------------------------------------
m <- 100
n <- 100

ep <- exceedProb(cutoff = cutoff_vec, 
                   theta_hat = theta_hat, 
                   sd_hat = sd_hat, 
                   alpha = alpha, 
                   d = d,
                   n = n,
                   m = m,
                   lower_tail = TRUE)

p_two_sided <- 2 * (1 - pt(sqrt(n) * abs(theta_hat - cutoff_vec) / sd_hat, df = n - 1))

exceed_ci <- data.frame(c = c(cutoff_vec, rev(cutoff_vec)), 
                        val = c(ep$upper, rev(ep$lower)))

exceed_p_val <- data.frame(cutoff_vec,
                           exceed = ep$point,
                           p_two_sided = p_two_sided)

exceed_p_val_m <- melt(exceed_p_val, id = "cutoff_vec")

dev.new(width = 7, height = 4.5)
ggplot() +
  geom_polygon(aes(x = c, y = val), data = exceed_ci, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = cutoff_vec, y = value, color = variable),
            data = exceed_p_val_m) +
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
                    labels = c(-1, expression(hat(theta)[L]), expression(hat(theta)[U]), 1)) +
  scale_color_manual("", values = c("grey", "black"),
                     labels = c(expression(paste("Pr",""[hat(theta)],","[hat(sigma)],"(", hat(theta)^"rep" <="c)", sep = "")),
                                "p-value")) +
  guides(color = guide_legend(override.aes = list(size = c(5, 0.6))))
ggsave(file.path(paper_path, "ci_p_val_overlay.png"))

# one-sided -----------------------------------------------
m <- 100
n <- 100

t_alpha_one_sided <- qt(p = 1 - alpha, df = n - d, lower.tail = TRUE)
ci_one_sided <- theta_hat - t_alpha_one_sided * sd_hat / sqrt(n)

ep <- exceedProb(cutoff = cutoff_vec, 
                   theta_hat = theta_hat, 
                   sd_hat = sd_hat, 
                   alpha = 2*alpha, 
                   d = d,
                   n = n,
                   m = m,
                   lower_tail = TRUE)

p_one_sided <- pt(sqrt(n) * (cutoff_vec - theta_hat) / sd_hat, df = n - 1)

exceed_ci <- data.frame(c = c(cutoff_vec, rev(cutoff_vec)), 
                        val = c(ep$upper, rev(ep$lower)))

exceed_p_val <- data.frame(cutoff_vec,
                           exceed = ep$point,
                           p_one_sided = p_one_sided)
exceed_p_val_m <- melt(exceed_p_val, id = "cutoff_vec")

dev.new(width = 7, height = 4.5)
ggplot() +
  geom_polygon(aes(x = c, y = val), data = exceed_ci, 
               fill = "grey", color = "grey", linetype = "solid") +
  geom_line(aes(x = cutoff_vec, y = value, color = variable),
            data = exceed_p_val_m) +
  theme_bw(17) +
  labs(x = "Cutoff c", 
       y = "") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = ci_one_sided, linetype = "dashed") +
  geom_hline(yintercept = c(0.05, 0.5), linetype = "dashed") +
  scale_y_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = c(-1, ci_one_sided, 1),
                    labels = c(-1, 
                               expression(hat(theta)[L]),
                               1)) +
  scale_color_manual("", values = c("grey", "black"),
                     labels = c(expression(paste("Pr",""[hat(theta)],","[hat(sigma)],"(", hat(theta)^"rep" <="c)", sep = "")),
                                "p-value")) +
  guides(color = guide_legend(override.aes = list(size = c(5, 0.6))))
ggsave(file.path(paper_path, "ci_p_val_overlay_one_sided.png"))
