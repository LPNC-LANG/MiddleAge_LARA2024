################################################################################
# Written by Clément Guichet, PhD Student
# LPNC - CNRS UMR 5105
# 2024

################################################################################

library(tidyverse)
library(rio)
library(stringr)
require(preprocessCore)
require(ggpubr)
library(readxl)

rm(list = ls())
source("E:/Research_Projects/MiddleAge_LARA/STATS/1_PreliminaryAnalysis/analysis/lifespan.R")
rm(list = ls()[!ls() %in% c("data_lifespan_full")])
# 613 subjects after excluding participants with >3 missing values
# see SENECA project for details on exclusion procedure

# Mass-univariate GAM -----------------------------------------------------
library(mgcv)
library(broom)

cognitive_var <- data_lifespan_full %>%
  dplyr::select(Cattell:Verbal_Fluency) %>%
  scale() %>%
  as.data.frame()

data_fitted <- list()
data_pvalue <- list()
data_gam_list <- list()
# 2 - Compute CURVATURE
derv_list <- list()
this_font_size <- 12

for (i in 1:ncol(cognitive_var)) {
  # Load the tract TW-FA values
  tmp <- cognitive_var[, i] %>% as.data.frame()
  colnames(tmp) <- colnames(cognitive_var)[i]
  print(colnames(tmp)) # Make sure the loop is working as expected
  # Fit the GAM model
  data_gam <- cbind(
    Age_Cog = data_lifespan_full$Age_Cog,
    tmp,
    handedness = data_lifespan_full$hand,
    MMSE = data_lifespan_full$MMSE,
    gender = data_lifespan_full$gender,
    TIV = data_lifespan_full$tiv_cubicmm
  ) %>% as.data.frame()
  data_gam_list[[i]] <- data_gam
  gam_tmp <- mgcv::gam(
    data_gam[, 2] ~ s(Age_Cog, bs = "ps") +
      handedness +
      MMSE +
      gender +
      TIV,
    data = data_gam,
    method = "REML"
  )
  # Store the fitted values
  data_fitted[[i]] <- broom::augment(gam_tmp) %>%
    as.data.frame() %>%
    dplyr::select(.fitted)
  data_pvalue[[i]] <- broom::tidy(gam_tmp) %>%
    as.data.frame() %>%
    dplyr::select(p.value)

  derv <- gratia::derivatives(gam_tmp, order = 2, interval = "confidence", unconditional = F)
  # add significance variable (true or false)
  derv <- derv %>%
    mutate(sig = !(0 > lower & 0 < upper)) # derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  # new variable with only significant derivatives (non-sig. ones are set to 0)
  derv$sig_deriv <- derv$derivative * derv$sig

  derv_list[[i]] <- (derv %>% subset(sig_deriv != 0))[, c("data", "derivative", "lower", "upper", "sig_deriv")]

  library(viridis)
  derv %>%
    ggplot(aes(x = data, y = sig_deriv)) +
    geom_point(size = 2) +
    geom_smooth() +
    scale_color_viridis(name = "Cognitive age", discrete = F, option = "magma") +
    theme_minimal(base_size = 18) +
    labs(
      x = "Age",
      y = "Cognitive decline"
    )

  # print changes range if significant
  if (all(derv$sig_deriv == 0)) {
    cat(sprintf("\n No significant change in %s \n", cog_var))
  } else {
    cat(sprintf("\nSig change: %1.2f - %1.2f\n", min(derv$data[derv$sig == T]), max(derv$data[derv$sig == T])))
  }

  # plot change
  derv[derv == 0] <- NA
  d <- ggplot(data = derv) +
    geom_tile(aes(x = data, y = .5, fill = sig_deriv)) +
    scale_x_continuous(breaks = seq(20, 90, 5)) +
    scale_fill_gradient(
      low = "darkblue", high = "darkorange", na.value = "white",
      limits = c(min(derv$sig_deriv), max(derv$sig_deriv))
    )

  d +
    labs(x = "Age (years)", fill = "Change \nper year") +
    theme(
      axis.title.y = element_blank(),
      # axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = this_font_size),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      text = element_text(size = this_font_size),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 30),
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
      ticks.linewidth = 2,
      ticks.colour = "black",
      draw.ulim = F,
      frame.colour = "black",
      frame.linetype = 1,
      frame.linewidth = 3,
      reverse = F,
      direction = "vertical",
      title.position = "top"
    )) +
    geom_rect(aes(ymin = 0, ymax = 1, xmin = min(data), xmax = max(data)), color = "black", fill = "white", alpha = 0)
}

# 1 - FDR correction
data_fitted_unlisted <- as.data.frame(do.call(cbind, data_fitted))
data_pvalue_unlisted <- as.data.frame(do.call(rbind, data_pvalue))
p.adjust(data_pvalue_unlisted[, 1], method = "fdr")

# 2 - Plot CURVATURE
derv_unlisted <- rbindlist(derv_list, fill = T, idcol = T) %>%
  dplyr::select(-.id) %>%
  group_by(data) %>%
  mutate_at(vars(derivative, upper, lower, sig_deriv), funs(mean(.))) %>%
  distinct(data, .keep_all = T)

ggplot(derv_unlisted, aes(x = data, y = derivative)) +
  geom_point(size = 2, colour = "black") +
  geom_ribbon(aes(x = data, y = derivative, ymin = lower, ymax = upper), inherit.aes = F, alpha = .7, color = "black", fill = "orange") +
  # geom_line(aes(x = data, y = derivative), inherit.aes = FALSE) +
  geom_line() +
  geom_smooth(color = "darkviolet", alpha = .7, size = 1.25, method = "gam", formula = y ~ s(x, bs = "ps")) +
  xlab("Age") +
  ylab("Curvature") +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = seq(45, 60, .02), color = "darkviolet", alpha = .04) +
  geom_jitter(height = 0.05, alpha = 0.08) +
  scale_x_continuous(breaks = seq(20, 90, 5)) +
  coord_cartesian(xlim = c(20, 90), ylim = c(-0.05, 0.01)) +
  scale_color_brewer(palette = "Paired") +
  theme_pubclean(base_size = 16)
