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

# Plot BSR ----

# For LC1
loading_names <- c(
  "Cattell",
  "Proverb",
  "Naming",
  "Multitasking",
  "Sentence\n Comprehension"
)

# For voxel
U <- c(
  -0.171761393844752,
  -0.106242562735677,
  -0.162321503586366,
  -0.200400870884538,
  -0.126095272814014
)

U_std <- c(
  0.0273412599196581,
  0.0319785714537446,
  0.0260943183570321,
  0.0304235632493853,
  0.0290905078309760
)

# For LC2
loading_names <- c(
  "Cattell",
  "Proverb",
  "Naming",
  "Sentence\n Comprehension",
  "Story Recall"
)

# For voxel
U <- c(
  -0.119555711579538,
  0.420042799340638,
  0.356627711158160,
  0.338764627731605,
  0.149715296509220
)

U_std <- c(
  0.0284520083340552,
  0.0355762835997668,
  0.0285196718663805,
  0.0300306486891984
)



plot_loadings <- cbind(U, U_std, loading_names) %>%
  as.data.frame() %>%
  mutate(score = as.numeric(U)) %>%
  mutate(
    ymin = score - as.numeric(U_std),
    ymax = score + as.numeric(U_std)
  )

plot_loadings$loading_names <- factor(plot_loadings$loading_names) %>%
  fct_reorder(plot_loadings$score, .desc = F)



library(ggchicklet)
ggplot(
  plot_loadings,
  aes(x = factor(loading_names), y = score)
) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .2) +
  geom_chicklet(radius = grid::unit(5, "mm"), aes(fill = score, alpha = .8, color = score)) +
  scale_fill_fermenter(palette = "Oranges", direction = 1) +
  # scale_y_continuous(limits = c(-0.25, -0.1), breaks = seq(-0.25, -0.1, 0.05)) +
  coord_flip() +
  # geom_text(aes(y = -0.1, label = loading_names, family = 'arial', fontface = "bold", hjust = -.5), size = 7) +
  theme_pubr(
    base_size = 16,
    legend = "none",
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank()
  )
