################################################################################
# Written by Clément Guichet, PhD Student
# LPNC - CNRS UMR 5105
# 2024

################################################################################

library(tidyverse)
library(readxl)
library(rio)
library(rlang)
library(mgcv)
library(stringr)
library(fitdistrplus)
library(ggpubr)
library(mgcViz)
library(cowplot)
library(dplyr)
library(olsrr)

rm(list = ls())
setwd("E:/Research_Projects/MiddleAge_LARA/STATS/2_Neuroanatomical_NMF/output/NMF_results")

# DATA WRANGLING ----
expansion_coeffs <- rio::import("expansion_coeffs.csv") %>%
  t() %>%
  as.data.frame()
demographics <- rio::import("E:/Research_Projects/MiddleAge_LARA/STATS/1_PreliminaryAnalysis/data_wrangling/AGE_COG_data_155Subj.csv")

# Confounds
TIV <- read_excel("E:/Research_Projects/MiddleAge_LARA/DWI/participant_data_DTI.xlsx", sheet = "TIV", col_names = F)
colnames(TIV) <- "TIV"

confounds_df <- cbind(demographics[, c(6:8)], TIV)
colnames(confounds_df) <- c("handedness", "gender", "MMSE", "TIV")


# GAM ----
library(mgcv)
library(broom)

data_fitted <- list() # encodes the fitted values for each network-specific GAM model
data_pvalue <- list() # same for the p value
data_gam_list <- list() # encodes the dataset used for the model

for (i in 1:ncol(expansion_coeffs)) {
  # Load the tract TW-FA values
  tmp <- expansion_coeffs[, i] %>% as.data.frame()
  colnames(tmp) <- colnames(expansion_coeffs)[i]
  print(colnames(tmp)) # Make sure the loop is working as expected
  # Fit the GAM model
  data_gam <- cbind(age = demographics$age, tmp) %>%
    scale() %>%
    as.data.frame()
  data_gam_list[[i]] <- data_gam
  gam_tmp <- mgcv::gam(
    data_gam[, 2] ~ s(age, k = 3)
      + gender
      + TIV
      + handedness
      + MMSE,
    data = cbind(data_gam, confounds_df), # combine with the covariates
    method = "REML"
  )
  # Store the fitted values
  data_fitted[[i]] <- broom::augment(gam_tmp) %>%
    as.data.frame() %>%
    dplyr::select(.fitted)
  data_pvalue[[i]] <- broom::tidy(gam_tmp) %>%
    as.data.frame() %>%
    dplyr::select(p.value)
}

data_fitted_unlisted <- as.data.frame(do.call(cbind, data_fitted))
data_pvalue_unlisted <- as.data.frame(do.call(rbind, data_pvalue))
# FDR correction
p.adjust(data_pvalue_unlisted[, 1], method = "fdr")

# 4 5 7 11 14 16


# Get Partial R2 ----
# Based on Chenying script (https://github.com/PennLINC/ModelArray_paper/blob/enh/figures/notebooks/utils.R#L12)

expansion_coeff_and_age <- cbind(age = demographics$age, confounds_df, expansion_coeffs %>% as.data.frame() %>% janitor::clean_names()) %>%
  dplyr::select(age:TIV, x1:x16)

# Get components numbering
comp_names <- expansion_coeff_and_age %>% dplyr::select(x1:x16)
Components <- names(comp_names)


# Full model
gamModels_age <- lapply(Components, function(x) {
  gam(substitute(i ~ s(age, k = 3) + gender + MMSE + TIV + handedness, list(i = as.name(x))), method = "REML", data = expansion_coeff_and_age)
})
# Reduced models
redmodel <- lapply(Components, function(x) {
  gam(substitute(i ~ gender + MMSE + TIV + handedness, list(i = as.name(x))), method = "REML", data = expansion_coeff_and_age)
})

partialRsq <- function(fullmodel, redmodel) {
  # calculating SSE: used observed y (i.e. excluding observations with NA), and fitted values, directly from model object

  sse.full <- sum((fullmodel$y - fullmodel$fitted.values)^2)
  sse.red <- sum((redmodel$y - redmodel$fitted.values)^2)

  partialRsq <- (sse.red - sse.full) / sse.red

  toReturn <- list(
    partialRsq = partialRsq,
    sse.full = sse.full,
    sse.red = sse.red
  )
  return(toReturn)
}

partialR2_V1 <- partialRsq(gamModels_age[[1]], redmodel[[1]])
partialR2_V2 <- partialRsq(gamModels_age[[2]], redmodel[[2]])
partialR2_V3 <- partialRsq(gamModels_age[[3]], redmodel[[3]])
partialR2_V4 <- partialRsq(gamModels_age[[4]], redmodel[[4]])
partialR2_V5 <- partialRsq(gamModels_age[[5]], redmodel[[5]])
partialR2_V6 <- partialRsq(gamModels_age[[6]], redmodel[[6]])
partialR2_V7 <- partialRsq(gamModels_age[[7]], redmodel[[7]])
partialR2_V8 <- partialRsq(gamModels_age[[8]], redmodel[[8]])
partialR2_V9 <- partialRsq(gamModels_age[[9]], redmodel[[9]])
partialR2_V10 <- partialRsq(gamModels_age[[10]], redmodel[[10]])
partialR2_V11 <- partialRsq(gamModels_age[[11]], redmodel[[11]])
partialR2_V12 <- partialRsq(gamModels_age[[12]], redmodel[[12]])
partialR2_V13 <- partialRsq(gamModels_age[[13]], redmodel[[13]])
partialR2_V14 <- partialRsq(gamModels_age[[14]], redmodel[[14]])
partialR2_V15 <- partialRsq(gamModels_age[[15]], redmodel[[15]])
partialR2_V16 <- partialRsq(gamModels_age[[16]], redmodel[[16]])

# Merge Partial R2 values with F-stats and p-values
partialR2 <- as.data.frame(cbind(
  partialR2_V1[[1]], partialR2_V2[[1]], partialR2_V3[[1]], partialR2_V4[[1]], partialR2_V5[[1]],
  partialR2_V6[[1]], partialR2_V7[[1]], partialR2_V8[[1]], partialR2_V9[[1]], partialR2_V10[[1]],
  partialR2_V11[[1]], partialR2_V12[[1]], partialR2_V12[[1]], partialR2_V14[[1]], partialR2_V15[[1]],
  partialR2_V16[[1]]
))

partialR2 <- as.data.frame(t(partialR2))


# Partial residual plots

network_order <- c(4, 5, 7, 11, 14, 16)
this_font_size <- 18

resid_plot <- function(i, add.intercept = FALSE) {
  data <- cbind(expansion_coeff_and_age[1:5], expansion_coeff_and_age[5 + network_order[i]])

  # covariate gam
  modobj <- mgcv::gam(
    data[, 6] ~ s(age, k = 3)
      + gender
      + MMSE
      + TIV
      + handedness,
    data = data,
    method = "REML"
  )
  df <- modobj$model
  mod.intercept <- modobj$coefficients["(Intercept)"]
  pterms <- predict(modobj, type = "terms", se.fit = TRUE)

  if (add.intercept == TRUE) {
    pterms.fit <- pterms$fit + mod.intercept
  } else {
    pterms.fit <- pterms$fit
  }
  pterms.sefit <- pterms$se.fit

  colnames(pterms.fit) <- gsub(x = colnames(pterms.fit), pattern = "s\\(", replacement = "") %>%
    gsub(pattern = "\\)", replacement = "")
  colnames(pterms.sefit) <- gsub(x = colnames(pterms.sefit), pattern = "s\\(", replacement = "") %>%
    gsub(pattern = "\\)", replacement = "")

  pterms.df <- data.frame(pterms.fit) %>%
    dplyr::select(age) %>%
    plyr::rename(c("age" = "fit")) %>%
    cbind(data.frame(pterms.sefit) %>%
      dplyr::select(age) %>%
      plyr::rename(c("age" = "se.fit"))) %>%
    mutate(
      upr = fit + 1.96 * se.fit,
      lwr = fit - 1.96 * se.fit
    )

  partial.residuals <- data.frame(pterms.fit) %>%
    mutate(across(.cols = everything(), .fns = function(x) {
      x + resid(modobj)
    })) %>%
    cbind(rawdata = df[, "age"])

  plot.df <- cbind(partial.residuals, pterms.df)

  # default plot
  ggplot(plot.df, aes(x = rawdata, y = fit)) +
    geom_point(size = 2, colour = "gray56") +
    geom_ribbon(aes(x = rawdata, y = fit, ymin = lwr, ymax = upr), inherit.aes = FALSE, alpha = .5) +
    geom_line(aes(x = rawdata, y = fit), inherit.aes = FALSE) +
    xlab("") +
    ylab("Network loading \n(Average TW-FA)") +
    scale_x_continuous(breaks = seq(45, 60, 1)) +
    theme_pubr(base_size = this_font_size)
}

get_derivs_and_plot <- function(i, low_color = NULL, hi_color = NULL, deriv_order = NULL, eps, lb, ub) {
  library(gratia)
  data <- cbind(expansion_coeff_and_age[1:5], expansion_coeff_and_age[5 + network_order[i]])

  # covariate gam
  modobj <- mgcv::gam(
    data[, 6] ~ s(age, k = 3)
      + gender
      + MMSE
      + TIV
      + handedness,
    data = data,
    method = "REML"
  )


  nmf_network <- Components[i] # get the network name

  # get model derivatives
  derv <- derivatives(modobj, interval = "confidence", unconditional = F, order = deriv_order, eps = eps) # from gratia. "confidence" for point-wise intervals

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
        limits = c(lb, ub)
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
}


get_derivs_and_plot(2, deriv_order = 1L, eps = .01, lb = -.5, ub = 0.25) # network 5
get_derivs_and_plot(6, deriv_order = 1L, eps = .01, lb = -.5, ub = 0.25) # network 16
get_derivs_and_plot(4, deriv_order = 1L, eps = .01, lb = -.5, ub = 0.25) # network 11
get_derivs_and_plot(1, deriv_order = 1L, eps = .01, lb = -.5, ub = 0.25) # network 4
get_derivs_and_plot(5, deriv_order = 1L, eps = .01, lb = -.5, ub = 0.25) # network 14
get_derivs_and_plot(3, deriv_order = 1L, eps = .01, lb = -.5, ub = 0.25) # network 7


# Figure 2 ----
dev_plots_list <- lapply(X = seq(1:6), FUN = resid_plot, add.intercept = TRUE)
dev_plots_list
bar_plots_list <- lapply(X = seq(1:6), FUN = get_derivs_and_plot, deriv_order = 1L, eps = .01, lb = NULL, ub = NULL)
bar_plots_list
