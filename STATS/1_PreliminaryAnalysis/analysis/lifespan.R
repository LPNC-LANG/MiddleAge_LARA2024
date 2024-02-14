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
setwd("E:/Research_Projects/MiddleAge_LARA/STATS/1_PreliminaryAnalysis/analysis/")

# GAMs --------------------------------------------------------------------
source("E:/Research_Projects/MiddleAge_LARA/STATS/1_PreliminaryAnalysis/data_wrangling/AGE_COG_data_for_DTI.R")
setwd("E:/Research_Projects/MiddleAge_LARA/STATS/1_PreliminaryAnalysis/analysis/")


data_lifespan <- CAMCAN_cognitive_data_imputed %>%
  mutate(
    ToT_Ratio_modified = 1 - ToT_Ratio,
    Hotel_Task_modified = 1 - log(Hotel_Task)
  )

TIV_lifespan <- rio::import("../data_wrangling/meta_data/participant_data_T1.xlsx") %>%
  filter(Observations %in% data_lifespan$Observations) %>%
  dplyr::select(Observations, tiv_cubicmm)

data_lifespan_full <- merge(data_lifespan, TIV_lifespan, by = "Observations") %>%
  mutate_at(vars(hand, tiv_cubicmm), funs(as.numeric(.)))


data_lifespan_full %>%
  pivot_longer(
    c(Cattell:Naming, Sentence_Comprehension_c:Hotel_Task_modified),
    names_to = "Cognitive_assessment",
    values_to = "performance"
  ) %>%
  mutate(Cognitive_assessment = ifelse(Cognitive_assessment == "Verbal_Fluency", "Verbal Fluency",
    ifelse(Cognitive_assessment == "Sentence_Comprehension_c", "Sentence Comp",
      ifelse(Cognitive_assessment == "Story_Recall", "Story Recall",
        ifelse(Cognitive_assessment == "ToT_Ratio_modified", "Tip-of-the-tongue",
          ifelse(Cognitive_assessment == "Hotel_Task_modified", "Hotel Task",
            Cognitive_assessment
          )
        )
      )
    )
  )) %>%
  group_by(Cognitive_assessment) %>%
  # perform quantile normalization   require(preprocessCore)
  mutate(performance = as.numeric(scale(normalize.quantiles(as.matrix(performance))))) %>%
  ggplot(aes(age, performance, color = Cognitive_assessment)) +
  geom_hline(yintercept = 0, color = "red") +
  geom_jitter(height = 0.05, alpha = 0.08) +
  geom_smooth(linewidth = 2, method = "gam", formula = y ~ s(x, bs = "ps"), alpha = .5) +
  scale_x_continuous(breaks = seq(20, 90, 10)) +
  coord_cartesian(ylim = c(-1, 1), xlim = c(20, 90)) + # Specify the age range for the study
  scale_color_brewer(palette = "Paired") +
  theme_pubr(
    base_size = 12,
    legend = "none"
  ) +
  theme(
    plot.title.position = "plot",
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank()
  ) +
  labs(y = "Normalized score") +
  facet_wrap(~Cognitive_assessment, scale = "free", ncol = 4)



# # Naming
# data_lifespan_full %>%
#   pivot_longer(
#     c(Cattell:Naming, Sentence_Comprehension_c:Hotel_Task_modified),
#     names_to = "Cognitive_assessment",
#     values_to = "performance"
#   ) %>%
#   mutate(Cognitive_assessment = ifelse(Cognitive_assessment == "Verbal_Fluency", "Verbal Fluency",
#                                        ifelse(Cognitive_assessment == "Sentence_Comprehension_c", "Sentence Comp",
#                                               ifelse(Cognitive_assessment == "Story_Recall", "Story Recall",
#                                                      ifelse(Cognitive_assessment == "ToT_Ratio_modified", "Tip-of-the-tongue",
#                                                             ifelse(Cognitive_assessment == "Hotel_Task_modified", "Hotel Task",
#                                                                    Cognitive_assessment)))))) %>%
#   group_by(Cognitive_assessment) %>%
#   #perform quantile normalization   require(preprocessCore)
#   mutate(performance = as.numeric(scale(normalize.quantiles(as.matrix(performance))))) %>%
#   filter(Cognitive_assessment == "Naming") %>%
#   ggplot(aes(age, performance)) +
#   geom_jitter(height = 0.05, alpha = 0) +
#   geom_smooth(linewidth = 2, method = "gam", formula = y ~ s(x, bs = 'ps'), alpha = 0, color = 'red') +
#   scale_x_continuous(breaks = seq(20,90,5)) +
#   coord_cartesian(ylim = c(-1, 1), xlim = c(20,90)) +  # Specify the age range for the study
#   scale_color_brewer(palette = "Paired") +
#   theme_pubr(base_size = 16,
#              legend = "none") +
#   theme(plot.title.position = "plot",
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_blank()) +
#   labs(y = "Normalized score")
#
# # Sen. Compr
# data_lifespan_full %>%
#   pivot_longer(
#     c(Cattell:Naming, Sentence_Comprehension_c:Hotel_Task_modified),
#     names_to = "Cognitive_assessment",
#     values_to = "performance"
#   ) %>%
#   mutate(Cognitive_assessment = ifelse(Cognitive_assessment == "Verbal_Fluency", "Verbal Fluency",
#                                        ifelse(Cognitive_assessment == "Sentence_Comprehension_c", "Sentence Comp",
#                                               ifelse(Cognitive_assessment == "Story_Recall", "Story Recall",
#                                                      ifelse(Cognitive_assessment == "ToT_Ratio_modified", "Tip-of-the-tongue",
#                                                             ifelse(Cognitive_assessment == "Hotel_Task_modified", "Hotel Task",
#                                                                    Cognitive_assessment)))))) %>%
#   group_by(Cognitive_assessment) %>%
#   #perform quantile normalization   require(preprocessCore)
#   mutate(performance = as.numeric(scale(normalize.quantiles(as.matrix(performance))))) %>%
#   filter(Cognitive_assessment == "Sentence Comp") %>%
#   ggplot(aes(age, performance)) +
#   geom_jitter(height = 0.05, alpha = 0) +
#   geom_smooth(linewidth = 2, method = "gam", formula = y ~ s(x, bs = 'ps'), alpha = 0, color = 'darkblue') +
#   scale_x_continuous(breaks = seq(20,90,5)) +
#   coord_cartesian(ylim = c(-1, 1), xlim = c(20,90)) +  # Specify the age range for the study
#   scale_color_brewer(palette = "Paired") +
#   theme_pubr(base_size = 16,
#              legend = "none") +
#   theme(plot.title.position = "plot",
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_blank()) +
#   labs(y = "Normalized score")
#
# # Multitasking
# data_lifespan_full %>%
#   pivot_longer(
#     c(Cattell:Naming, Sentence_Comprehension_c:Hotel_Task_modified),
#     names_to = "Cognitive_assessment",
#     values_to = "performance"
#   ) %>%
#   mutate(Cognitive_assessment = ifelse(Cognitive_assessment == "Verbal_Fluency", "Verbal Fluency",
#                                        ifelse(Cognitive_assessment == "Sentence_Comprehension_c", "Sentence Comp",
#                                               ifelse(Cognitive_assessment == "Story_Recall", "Story Recall",
#                                                      ifelse(Cognitive_assessment == "ToT_Ratio_modified", "Tip-of-the-tongue",
#                                                             ifelse(Cognitive_assessment == "Hotel_Task_modified", "Hotel Task",
#                                                                    Cognitive_assessment)))))) %>%
#   group_by(Cognitive_assessment) %>%
#   #perform quantile normalization   require(preprocessCore)
#   mutate(performance = as.numeric(scale(normalize.quantiles(as.matrix(performance))))) %>%
#   filter(Cognitive_assessment == "Hotel Task") %>%
#   ggplot(aes(age, performance)) +
#   geom_jitter(height = 0.05, alpha = 0) +
#   geom_smooth(linewidth = 2, method = "gam", formula = y ~ s(x, bs = 'ps'), alpha = 0, color = 'orange') +
#   scale_x_continuous(breaks = seq(20,90,5)) +
#   coord_cartesian(ylim = c(-1, 1), xlim = c(20,90)) +  # Specify the age range for the study
#   scale_color_brewer(palette = "Paired") +
#   theme_pubr(base_size = 16,
#              legend = "none") +
#   theme(plot.title.position = "plot",
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_blank()) +
#   labs(y = "Normalized score")
#
# # Fluid intelligence
# data_lifespan_full %>%
#   pivot_longer(
#     c(Cattell:Naming, Sentence_Comprehension_c:Hotel_Task_modified),
#     names_to = "Cognitive_assessment",
#     values_to = "performance"
#   ) %>%
#   mutate(Cognitive_assessment = ifelse(Cognitive_assessment == "Verbal_Fluency", "Verbal Fluency",
#                                        ifelse(Cognitive_assessment == "Sentence_Comprehension_c", "Sentence Comp",
#                                               ifelse(Cognitive_assessment == "Story_Recall", "Story Recall",
#                                                      ifelse(Cognitive_assessment == "ToT_Ratio_modified", "Tip-of-the-tongue",
#                                                             ifelse(Cognitive_assessment == "Hotel_Task_modified", "Hotel Task",
#                                                                    Cognitive_assessment)))))) %>%
#   group_by(Cognitive_assessment) %>%
#   #perform quantile normalization   require(preprocessCore)
#   mutate(performance = as.numeric(scale(normalize.quantiles(as.matrix(performance))))) %>%
#   filter(Cognitive_assessment == "Cattell") %>%
#   ggplot(aes(age, performance)) +
#   geom_jitter(height = 0.05, alpha = 0) +
#   geom_smooth(linewidth = 2, method = "gam", formula = y ~ s(x, bs = 'ps'), alpha = 0, color = 'orange') +
#   scale_x_continuous(breaks = seq(20,90,5)) +
#   coord_cartesian(ylim = c(-1, 1), xlim = c(20,90)) +  # Specify the age range for the study
#   scale_color_brewer(palette = "Paired") +
#   theme_pubr(base_size = 16,
#              legend = "none") +
#   theme(plot.title.position = "plot",
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         axis.title.x = element_blank()) +
#   labs(y = "Normalized score")
