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

rm(list=ls())


# Fiber connectivity ------------------------------------------------------

setwd("E:/Research_Projects/MiddleAge_LARA/STATS/2_Neuroanatomical_NMF/output/NIFTI/networks_nifti_cluster")

track_list <- rio::import('list_bundles.txt')
track_list_clean <- rbind(colnames(track_list), track_list)
colnames(track_list_clean) <- 'track_name'


listfile <- list.files(getwd(), pattern = "*connectivity_overlap.tsv")
# reorder the files list to be in the increasing order
listfile
network_order = c(4, 5, 7, 11, 14, 16)

files <- ldply(listfile, read.table, header = T, sep = '\t') %>% janitor::clean_names()
files <- cbind(files[,c(4:6)], files[,c(1:3)])
  
df_networks_list = list()
for (i in 1:6){
  tmp <- files[,i] %>% as.data.frame() %>% na.omit()
  colnames(tmp) <- colnames(track_list)
  
  tmp_clean <- 
    # Remove the row of total count
    rbind(colnames(tmp), tmp) %>% 
    mutate(helper_vector = rep(seq(3), times = 72)) %>% 
    subset(helper_vector != 3) %>% 
    # Align the track name with the count
    mutate(track_file = lag(.[,1]))  %>% 
    # Remove rows every two steps (those who are unaligned)
    .[- seq(1, nrow(.), 2),c(1:3)] %>% 
    # Put back clean name
    cbind(., track_list_clean) %>% 
    # Extract numbers
    mutate(count = parse_number(.$AF_left)) %>% 
    dplyr::select(track_name, count)
  
  df_networks_list[[i]] <- tmp_clean %>% 
    mutate(percentage = scale((as.numeric(count)/2000)*100)) %>%
    # coefficient for composite score
    mutate(percentage = percentage*(5/4)) %>%
    plyr::rename(c("count" = sprintf("network_%d", network_order[i]))) %>%
    plyr::rename(c("percentage" = sprintf("percentage_conn_for_network_%d", network_order[i])))
    
}

df_networks_unlisted <-  as.data.frame(do.call(cbind, df_networks_list)) # unscaled percentages


# Spatial overlap bundle-focused ---------------------------------------------------------
library(readr)

listfile <- list.files(getwd(), pattern = "*bundlefoc.tsv")
# reorder the files list to be in the increasing order
listfile
network_order = c(4, 5, 7, 11, 14, 16)

df_networks_list2 = list()
for (i in 1:6){
  tmp <- read_table(sprintf("%d_spatial_overlap_bundlefoc.tsv", network_order[i]), col_names = FALSE) %>% dplyr::select(X9)
  
  df_networks_list2[[i]] <- tmp %>% 
    mutate(percentage_spatial_bundle = X9/lag(X9)*100) %>% 
    .[-seq(1,nrow(.),2),2] %>% 
    scale() %>% 
    as.data.frame() %>% 
    mutate(track_name = track_list_clean$track_name) %>% 
    plyr::rename(c("percentage_spatial_bundle" = sprintf("percentage_spatial_bundlefocused_for_network_%d", network_order[i])))
}

df_networks_unlisted2 <-  as.data.frame(do.call(cbind, df_networks_list2))


# Spatial overlap network-focused ---------------------------------------------------------
library(readr)

listfile <- list.files(getwd(), pattern = "*netfoc.tsv")
# reorder the files list to be in the increasing order
listfile
network_order = c(4, 5, 7, 11, 14, 16)

df_networks_list3 = list()
for (i in 1:6){
  tmp <- read_table(sprintf("%d_spatial_overlap_netfoc.tsv", network_order[i]), col_names = FALSE) %>% dplyr::select(X9)

  df_networks_list3[[i]] <- tmp %>%
    mutate(percentage_spatial_net = X9/lag(X9)*100) %>%
    .[-seq(1,nrow(.),2),2] %>%
    scale() %>%
    as.data.frame() %>%
    # coefficient for composite score
    mutate(percentage_spatial_net = percentage_spatial_net*(1/2)) %>% # minimizing the influence or very large pathways like the corpus callosum
    mutate(track_name = track_list_clean$track_name) %>%
    plyr::rename(c("percentage_spatial_net" = sprintf("percentage_spatial_netfocused_for_network_%d", network_order[i])))
}

df_networks_unlisted3 <-  as.data.frame(do.call(cbind, df_networks_list3))

# Merge Spatial overlap and Connectivity ----------------------------------

network_composition_list <- list()
for (i in seq(1:6)){
  percentage_conn <- df_networks_unlisted %>% 
    janitor::clean_names() %>% 
    arrange(track_name) %>% 
    .[,c(1, 3*i)]
  
  percentage_spatial_bundlefocused <- df_networks_unlisted2 %>% 
    janitor::clean_names() %>%
    arrange(track_name) %>% 
    .[,c(seq(1,12,2))] %>% 
    cbind(., track_name = df_networks_unlisted2[,2]) %>% 
    .[,c(i, 7)]
  
  percentage_spatial_netfocused <- df_networks_unlisted3 %>%
    janitor::clean_names() %>%
    arrange(track_name) %>%
    .[,c(seq(1,12,2))] %>%
    cbind(., track_name = df_networks_unlisted3[,2]) %>%
    .[,c(i, 7)]

  network_composition_list[[i]] <- merge(percentage_conn, percentage_spatial_bundlefocused, by = "track_name") %>% 
    merge(., percentage_spatial_netfocused, by = "track_name") %>%
    dplyr::mutate(composite = rowSums(across(where(is.numeric)))/3) %>% 
    relocate(composite, .after = track_name) %>% 
    plyr::rename(c("composite" = sprintf("composite_%d", network_order[i])))
}

network_composition_unlisted = as.data.frame(do.call(cbind,network_composition_list))


# library(writexl)
# write_xlsx(network_composition_unlisted, "NMF_networks_composition.xlsx") 
