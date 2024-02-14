################################################################################
# Written by Clément Guichet, PhD Student
# LPNC - CNRS UMR 5105
# 2024

################################################################################

library(tidyverse)
library(data.table)
library(plyr)
library(janitor)
library(jsonlite)
library(ggpubr)
library(wordcloud2)

rm(list = ls())
setwd("E:/Research_Projects/MiddleAge_LARA/STATS/3_Neurocognitive_PLS/voxel_level_PLS/output/postproc/LC_composition/NIFTI")

# Fiber connectivity ------------------------------------------------------
track_list <- rio::import("list_bundles.txt")

track_list_clean <- rbind(colnames(track_list), track_list)
colnames(track_list_clean) <- "track_name"

# Change accordingly ----
setwd("./LC_x_tsv_files")

listfile <- list.files(getwd(), pattern = "LC1_percentage_tracks_cluster1.tsv")
# listfile <- list.files(getwd(), pattern = "LC1_percentage_tracks_cluster2.tsv")
# listfile <- list.files(getwd(), pattern = "LC1_percentage_tracks_cluster3.tsv")
#
# listfile <- list.files(getwd(), pattern = "LC2_percentage_tracks_cluster1.tsv")
# listfile <- list.files(getwd(), pattern = "LC2_percentage_tracks_cluster2.tsv")

files <- ldply(listfile, read.table, header = T, sep = "\t") %>% janitor::clean_names()


tmp <- files %>%
  as.data.frame() %>%
  na.omit()
colnames(tmp) <- colnames(track_list)

tmp_clean <-
  # Remove the row of total count
  rbind(colnames(tmp), tmp) %>%
  mutate(helper_vector = rep(seq(3), times = 72)) %>%
  subset(helper_vector != 3) %>%
  # Align the track name with the count
  mutate(track_file = lag(.[, 1])) %>%
  # Remove rows every two steps (those who are unaligned)
  .[-seq(1, nrow(.), 2), c(1:3)] %>%
  # Put back clean name
  cbind(., track_list_clean) %>%
  # Extract numbers
  mutate(count = parse_number(.[, 1])) %>%
  dplyr::select(track_name, count)

track_overlap <- tmp_clean %>%
  mutate(percentage_conn = scale((as.numeric(count) / 2000) * 100)) %>%
  dplyr::select(track_name, percentage_conn)


# Spatial overlap bundle-focused ---------------------------------------------------------
library(readr)
tmp <- read_table("LC1_spatial_percentage_tracks_cluster1.tsv", col_names = FALSE) %>% dplyr::select(X9)
# tmp <- read_table("LC1_spatial_percentage_tracks_cluster2.tsv", col_names = FALSE) %>% dplyr::select(X9)
# tmp <- read_table("LC1_spatial_percentage_tracks_cluster3.tsv", col_names = FALSE) %>% dplyr::select(X9)
#
# tmp <- read_table("LC2_spatial_percentage_tracks_cluster1.tsv", col_names = FALSE) %>% dplyr::select(X9)
# tmp <- read_table("LC2_spatial_percentage_tracks_cluster2.tsv", col_names = FALSE) %>% dplyr::select(X9)

spatial_overlap_bundle <- tmp %>%
  mutate(percentage_spatial_bundle = X9 / lag(X9) * 100) %>%
  .[-seq(1, nrow(.), 2), 2] %>%
  scale() %>%
  as.data.frame() %>%
  mutate(track_name = track_list_clean$track_name)

# Spatial overlap network-focused ---------------------------------------------------------
library(readr)
tmp <- read_table("LC1_spatial_percentage_network_cluster1.tsv", col_names = FALSE) %>% dplyr::select(X9)
# tmp <- read_table("LC1_spatial_percentage_network_cluster2.tsv", col_names = FALSE) %>% dplyr::select(X9)
# tmp <- read_table("LC1_spatial_percentage_network_cluster3.tsv", col_names = FALSE) %>% dplyr::select(X9)
#
# tmp <- read_table("LC2_spatial_percentage_network_cluster1.tsv", col_names = FALSE) %>% dplyr::select(X9)
# tmp <- read_table("LC2_spatial_percentage_network_cluster2.tsv", col_names = FALSE) %>% dplyr::select(X9)

spatial_overlap_network <- tmp %>%
  mutate(percentage_spatial_network = X9 / lag(X9) * 100) %>%
  .[-seq(1, nrow(.), 2), 2] %>%
  scale() %>%
  as.data.frame() %>%
  mutate(track_name = track_list_clean$track_name)

merge_overlap <- merge(spatial_overlap_bundle, spatial_overlap_network, by = "track_name") %>%
  merge(., track_overlap, by = "track_name") %>%
  mutate(composite = (percentage_spatial_bundle + percentage_spatial_network + percentage_conn) / 3) %>%
  relocate(composite, .after = track_name)
