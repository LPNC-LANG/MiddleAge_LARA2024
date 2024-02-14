################################################################################
# Written by Clément Guichet, PhD Student
# LPNC - CNRS UMR 5105
# 2024

################################################################################

library(tidyverse)
library(readxl)
library(data.table)
library(ggpubr)
library(viridis)

rm(list = ls())

# DATA WRANGLING ----------------------------------------------------------

age <- rio::import("E:/Research_Projects/MiddleAge_LARA/STATS/1_PreliminaryAnalysis/data_wrangling/AGE_COG_data_155Subj.csv")$age
gender <- rio::import("E:/Research_Projects/MiddleAge_LARA/STATS/1_PreliminaryAnalysis/data_wrangling/AGE_COG_data_155Subj.csv")$gender_code

latent_scores <- rio::import("latent_scores.csv") %>%
  mutate(
    Lx1_25 = Lx1_25 * (-1), # To put the 56-60 age group in the positive values
    Ly1_25 = Ly1_25 * (-1),
    net_Lx1_25 = net_Lx1_25 * (-1), # To put the 56-60 age group in the positive values
    net_Ly1_25 = net_Ly1_25 * (-1)
  ) %>%
  cbind(., age) %>%
  mutate(cog_age_groups = ifelse(age < 51, "young",
    ifelse(age < 56 & age >= 51, "middle",
      "old"
    )
  ))

# DATA ANALYSIS -----------------------------------------------------------

# LC1 for voxel
latent_scores %>%
  ggplot(aes(x = Lx1_25, y = Ly1_25, color = cog_age_groups)) +
  geom_smooth(method = "lm", color = "black", size = 2, alpha = 0.3) +
  geom_point(aes(color = cog_age_groups), size = 5) +
  scale_color_manual(
    values = c("#E66101", "darkviolet", "#FDB863"),
    labels = c("45-50", "51-55", "56-60")
  ) +
  scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 2)) +
  # scale_color_viridis(name = 'Cognitive age', discrete = F, option = "magma")+
  labs(
    x = "Latent voxel",
    y = "Latent cognitive",
    color = "Age groups",
    title = "LC1 - Voxel-level PLS"
  ) +
  theme_pubr(base_size = 16)

# LC1 for network
latent_scores %>%
  ggplot(aes(x = net_Lx1_25, y = net_Ly1_25, color = cog_age_groups)) +
  geom_smooth(method = "lm", color = "black", size = 2, alpha = 0.3) +
  geom_point(aes(color = cog_age_groups), size = 5) +
  scale_color_manual(
    values = c("#E66101", "darkviolet", "#FDB863"),
    labels = c("45-50", "51-55", "56-60")
  ) +
  scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 2)) +
  # scale_color_viridis(name = 'Cognitive age', discrete = F, option = "magma")+
  labs(
    x = "Latent network",
    y = "Latent cognitive",
    title = "LC1 - Network-level PLS"
  ) +
  theme_pubr(base_size = 16)

# LC2 voxel
latent_scores %>%
  ggplot(aes(x = Lx2_25, y = Ly2_25, color = cog_age_groups)) +
  geom_smooth(method = "lm", color = "black", size = 2, alpha = 0.3) +
  geom_point(aes(color = cog_age_groups), size = 5) +
  scale_color_manual(
    values = c("#E66101", "darkviolet", "#FDB863"),
    labels = c("45-50", "51-55", "56-60")
  ) +
  scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 2)) +
  # scale_color_viridis(name = 'Cognitive age', discrete = F, option = "magma")+
  labs(
    x = "Latent network",
    y = "Latent cognitive",
    title = "LC2 - Voxel-level PLS"
  ) +
  theme_pubr(base_size = 16)




# Combining both latent components ----

library(FactoMineR)
library(ggpubr)

latent_scores$pca_L1 <- (FactoMineR::PCA(latent_scores[, c(1, 2)]))$ind$coord[, 1] * (-1)
latent_scores$pca_netL1 <- (FactoMineR::PCA(latent_scores[, c(5, 6)]))$ind$coord[, 1] * (-1)
latent_scores$pca_L2 <- (FactoMineR::PCA(latent_scores[, c(3, 4)]))$ind$coord[, 1]

mod <- mgcv::gam(pca_L2 ~ s(age), data = latent_scores)
summary(mod)
plot(mod)

flexplot::flexplot(pca_netL1 ~ age, data = latent_scores, method = "loess")
mod <- mgcv::gam(pca_netL1 ~ s(age), data = latent_scores)
summary(mod)
plot(mod)

flexplot::flexplot(pca_L1 ~ age, data = latent_scores, method = "loess")
mod <- mgcv::gam(pca_L1 ~ s(age), data = latent_scores)
summary(mod)
plot(mod)

# Figure 3A ----
latent_scores %>%
  ggplot(aes(x = pca_L2, y = pca_L1, color = cog_age_groups)) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  # geom_smooth(method = "loess", color = 'black', size = 2, alpha = 0.3) +
  geom_point(aes(color = cog_age_groups), size = 5) +
  scale_color_manual(values = c("#E66101", "darkviolet", "#FDB863")) +
  # scale_color_viridis(name = 'Cognitive age', discrete = T, option = "magma")+
  theme_pubclean(base_size = 16) +
  labs(
    y = "LC1 (domain-general)",
    x = "LC2 (language-specific)"
  )



# Figure 5 ----
library(ggpubr)

latent_scores$interplay <- (FactoMineR::PCA(latent_scores[, c(9, 11)]))$ind$coord[, 1]

latent_scores_plot <- latent_scores %>%
  dplyr::select(age, pca_L2, pca_L1, interplay) %>%
  pivot_longer(c("pca_L2", "pca_L1", "interplay"), names_to = "latent", values_to = "Values")

latent_scores_plot$latent <- factor(latent_scores_plot$latent, levels = c("pca_L2", "pca_L1", "interplay"))

latent_scores_plot %>%
  ggplot(aes(age, Values, color = latent)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(height = 0.05, alpha = 0) +
  geom_smooth(linewidth = 3, method = "gam", formula = y ~ s(x, k = 3), alpha = 0) +
  scale_x_continuous(breaks = seq(45, 60, 1)) +
  scale_y_continuous(breaks = seq(-1.5, 1, 0.5)) +
  coord_cartesian(ylim = c(-1.5, 1)) +
  scale_color_manual(values = c("#FDCC8A", "#08519C", "red", "black")) +
  theme_pubr(
    base_size = 12,
    legend = "none",
  ) +
  theme(plot.title.position = "plot") +
  labs(y = "", x = "") +
  # labs(y = "Latent component score", x = "Age") +
  ggtitle("")


# Derivatives ----

library(gratia)
ub <- 0.01
lb <- -0.05
this_font_size <- 12
modobj <- mgcv::gam(interplay ~ s(age, k = 3),
  data = latent_scores,
  method = "REML"
)

# get model derivatives
derv <- derivatives(modobj, interval = "confidence", unconditional = F, order = 2L, eps = 1e-6) # from gratia. "confidence" for point-wise intervals

# add significance variable (true or false)
derv <- derv %>%
  mutate(sig = !(0 > lower & 0 < upper)) # derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
# new variable with only significant derivatives (non-sig. ones are set to 0)
derv$sig_deriv <- derv$derivative * derv$sig

# print changes range if significant
if (all(derv$sig_deriv == 0)) {
  cat(sprintf("\n No significant change in %s \n", nmf_network))
} else {
  cat(sprintf("\nSig change: %1.2f - %1.2f\n", min(derv$data[derv$sig == T]), max(derv$data[derv$sig == T])))
}

# plot change
derv[derv == 0] <- NA
if (is.null(lb) & is.null(ub)) {
  d <- ggplot(data = derv) +
    geom_tile(aes(x = data, y = .5, fill = sig_deriv)) +
    scale_fill_gradient(
      low = "darkblue", high = "darkorange", na.value = "white",
      limits = c(min(derv$sig_deriv), max(derv$sig_deriv))
    )
} else {
  d <- ggplot(data = derv) +
    geom_tile(aes(x = data, y = .5, fill = sig_deriv)) +
    scale_fill_gradient(
      low = "darkblue", high = "darkorange", na.value = "white",
      limits = c(min(derv$sig_deriv), max(derv$sig_deriv))
    )
}


d +
  coord_cartesian(xlim = c(45, 60)) +
  labs(x = "", fill = "") +
  scale_x_continuous(breaks = seq(45, 60, 1)) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = this_font_size),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    text = element_text(size = this_font_size),
    legend.text = element_text(size = this_font_size),
    legend.title = element_text(size = this_font_size),
    axis.title = element_text(size = this_font_size),
    legend.key.width = unit(1, "cm"),
    legend.position = "right",
    plot.margin = unit(c(0, 0, 0.5, 0), "cm"), # Top, left,Bottom, right
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(colour = "black", size = 1.5),
    axis.line.y = element_line(colour = "black", size = 1.5),
    axis.ticks.length = unit(.25, "cm"),
    axis.text = element_text(size = 50)
  ) +
  guides(fill = guide_colorbar(
    ticks = T,
    ticks.linewidth = 1,
    ticks.colour = "black",
    draw.ulim = T,
    frame.colour = "black",
    frame.linetype = 1,
    frame.linewidth = 1,
    reverse = T,
    direction = "horizontal",
    title.position = "top"
  )) +
  geom_rect(aes(ymin = 0, ymax = 1, xmin = min(data), xmax = max(data)), color = "black", fill = "white", alpha = 0)




# Graphical abstract ----
library(ggpubr)

latent_scores$interplay <- (FactoMineR::PCA(latent_scores[, c(9, 11)]))$ind$coord[, 1]

latent_scores_plot <- latent_scores %>%
  dplyr::select(age, pca_L2, pca_L1, interplay) %>%
  pivot_longer(c("pca_L2", "pca_L1", "interplay"), names_to = "latent", values_to = "Values")

latent_scores_plot$latent <- factor(latent_scores_plot$latent, levels = c("pca_L2", "pca_L1", "interplay"))

latent_scores_plot %>%
  filter(latent == "pca_L2") %>%
  ggplot(aes(age, Values, color = latent)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(height = 0.05, alpha = 0) +
  geom_smooth(linewidth = 3, method = "gam", formula = y ~ s(x, k = 3), alpha = 0) +
  scale_x_continuous(breaks = seq(45, 60, 1)) +
  scale_y_continuous(breaks = seq(-1.5, 1, 0.5)) +
  coord_cartesian(ylim = c(-1.5, 1)) +
  scale_color_manual(values = c("#FDCC8A")) +
  theme_pubr(
    base_size = 14,
    legend = "none",
  ) +
  theme(plot.title.position = "plot") +
  labs(y = "", x = "") +
  # labs(y = "Latent component score", x = "Age") +
  ggtitle("")

latent_scores_plot %>%
  filter(latent == "pca_L1") %>%
  ggplot(aes(age, Values, color = latent)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(height = 0.05, alpha = 0) +
  geom_smooth(linewidth = 3, method = "gam", formula = y ~ s(x, k = 3), alpha = 0) +
  scale_x_continuous(breaks = seq(45, 60, 1)) +
  scale_y_continuous(breaks = seq(-1.5, 1, 0.5)) +
  coord_cartesian(ylim = c(-1.5, 1)) +
  scale_color_manual(values = c("#08519C")) +
  theme_pubr(
    base_size = 14,
    legend = "none",
  ) +
  theme(plot.title.position = "plot") +
  labs(y = "", x = "") +
  # labs(y = "Latent component score", x = "Age") +
  ggtitle("")

latent_scores_plot %>%
  filter(latent == "interplay") %>%
  ggplot(aes(age, Values, color = latent)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(height = 0.05, alpha = 0) +
  geom_smooth(linewidth = 3, method = "gam", formula = y ~ s(x, k = 3), alpha = 0) +
  scale_x_continuous(breaks = seq(45, 60, 1)) +
  scale_y_continuous(breaks = seq(-1.5, 1, 0.5)) +
  coord_cartesian(ylim = c(-1.5, 1)) +
  scale_color_manual(values = c("black")) +
  theme_pubr(
    base_size = 14,
    legend = "none",
  ) +
  theme(plot.title.position = "plot") +
  labs(y = "", x = "") +
  # labs(y = "Latent component score", x = "Age") +
  ggtitle("")



# SUPPLEMENTARY RESULTS -----

# GENDER ----
latent_scores %>%
  mutate(gender = as.factor(gender)) %>%
  ggplot(aes(x = pca_L2, y = pca_L1, color = gender)) +
  geom_hline(yintercept = 0, size = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, size = 1, linetype = "dashed") +
  geom_point(aes(color = gender), size = 5) +
  scale_color_manual(
    values = c("-0.5" = "#00AFBB", "0.5" = "#FC4E07"),
    labels = c("Male", "Female")
  ) +
  theme_pubclean(base_size = 16) +
  labs(
    y = "LC1 (domain-general)",
    x = "LC2 (language-specific)",
    color = "Gender"
  )



# MALES
latent_score_gender <- latent_scores[, 1:8] %>%
  cbind(., gender)

latent_score_gender$pca_L1 <- (FactoMineR::PCA(latent_score_gender[, c(1, 2)]))$ind$coord[, 1] * (-1)
latent_score_gender$pca_L2 <- (FactoMineR::PCA(latent_score_gender[, c(3, 4)]))$ind$coord[, 1]
latent_score_gender$LP <- (FactoMineR::PCA(latent_score_gender[, c(10, 11)]))$ind$coord[, 1]

latent_score_gender_plot <- latent_score_gender %>%
  dplyr::select(age, gender, pca_L2, pca_L1, LP) %>%
  pivot_longer(c("pca_L2", "pca_L1", "LP"), names_to = "latent", values_to = "Values")

latent_score_gender_plot$latent <- factor(latent_score_gender_plot$latent, levels = c("pca_L2", "pca_L1", "LP"))


male_plot <- latent_score_gender_plot %>%
  filter(gender == -0.5) %>%
  ggplot(aes(age, Values, color = latent)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(height = 0.05, alpha = 0) +
  geom_smooth(linewidth = 2, method = "gam", formula = y ~ s(x, k = 3), alpha = 0.1) +
  scale_x_continuous(breaks = seq(45, 60, 1)) +
  scale_y_continuous(breaks = seq(-1.5, 1.7, 0.8)) +
  coord_cartesian(ylim = c(-1.5, 1.7)) +
  scale_color_manual(values = c("#FDCC8A", "#08519C", "red")) +
  theme_pubr(
    base_size = 16,
    legend = "none",
  ) +
  theme(plot.title.position = "plot") +
  labs(y = "", x = "") +
  # labs(y = "Latent component score", x = "Age") +
  ggtitle("")

male_plot

# FEMALES
female_plot <- latent_score_gender_plot %>%
  filter(gender == 0.5) %>%
  ggplot(aes(age, Values, color = latent)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_jitter(height = 0.05, alpha = 0) +
  geom_smooth(linewidth = 2, method = "gam", formula = y ~ s(x, k = 3), alpha = 0.1) +
  scale_x_continuous(breaks = seq(45, 60, 1)) +
  scale_y_continuous(breaks = seq(-1.5, 1.7, 0.8)) +
  coord_cartesian(ylim = c(-1.5, 1.7)) +
  scale_color_manual(values = c("#FDCC8A", "#08519C", "red")) +
  theme_pubr(
    base_size = 16,
    legend = "none",
  ) +
  theme(plot.title.position = "plot") +
  labs(y = "", x = "") +
  # labs(y = "Latent component score", x = "Age") +
  ggtitle("")

female_plot
