################################################################################
# Written by Clément Guichet, PhD Student
# LPNC - CNRS UMR 5105 
# 2024

################################################################################

library(rio)
library(tidyverse)

rm(list=ls())
set.seed(1551)

df <- rio::import("../data/template_data_age.xlsx", sheet = "participant_data_DTI")[,1:2]

# save mean and sd the subsamples should have
aimed_mean <- mean(df$age)
aimed_sd <- sd(df$age)

# set number of replications / iterations
n_replication <- 10000

# set size of sample
size_sample <- 35

# set desired number of samples
n_sample <- 1

# set deviation from mean and sd you can accept
deviation_mean <- .01*aimed_mean
deviation_sd <- .01*aimed_sd

# create empty container for resulting samples
samples <- list(n_replication)

# Repeatedly sample from length
i <- 0
sample_count <- 0

repeat {
  
  i <- i+1
  
  # take a sample from length
  sample_obs <- as.data.frame(sample(df$Observations, size_sample, replace = F)) %>% dplyr::rename(Observations = `sample(df$Observations, size_sample, replace = F)`)
  sample_length <- left_join(sample_obs, df, by = "Observations")
  
  
  # keep the sample when is is close enough
  if(abs(aimed_mean - mean(sample_length$age)) < deviation_mean &
     abs(aimed_sd - sd(sample_length$age)) < deviation_sd){
    
    samples[[i]] <- sample_length$Observations
    sample_count <- sample_count + 1
    
  }
  
  if(i == n_replication | sample_count == n_sample){
    break
  }
  
}

# your samples
sampled_df <- unlist(samples) %>% as.data.frame()
colnames(sampled_df) <- "subjects"

# write_tsv(sampled_df %>% arrange(subjects), "template_sample.txt")

