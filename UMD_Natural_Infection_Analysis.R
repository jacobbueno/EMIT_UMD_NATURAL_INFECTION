# EMIT UMD Natural Infection Analysis
# Program Objective: Run analysis on EMIT Natural Infection datasets
# Author: Jacob Bueno de Mesquita using material from Jing Yan and Don Milton
# Date: December 14, 2018 - February 2019
# Summary: While most of the analysis for the EMIT_UMD study was done in SAS (especially that as part of the PNAS manuscript), here are some others analyses that we will include as well. 

#### Load required packages and set working directory ####

library(tidyverse)
library(RcppRoll)
library(readxl)
library(knitr)
library(data.table)
library(lubridate)
library(arsenal)
library(lme4)

setwd("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection")

sessionInfo() # for reproducibility

#### Estimate the mean and variance for the ratio of ffu/copy number in fine aerosols from EMIT ####
# Use a mixed model with random effect of person and fixed effect of day since onset of symptoms. 
# First have to cut the right df

PNAS_data_full <- read.csv("Curated Data/Analytical Datasets/PNAS_data_full.csv") %>%
  select(-X)
PNAS_data_full$date.visit <- as.Date(PNAS_data_full$date.visit)

# Find each sample where we have PCR and focus assay

focus_subjects <- PNAS_data_full %>%
  filter(!is.na(focus.ct)) %>%
  distinct(subject.id, date.visit, .keep_all = TRUE) %>%
  distinct(subject.id)

focus <- PNAS_data_full %>%
  filter(!is.na(focus.ct))

focus_pos_subjects <- focus %>%
  filter(focus.ct > 0) %>%
  distinct(subject.id, date.visit, .keep_all = TRUE) %>%
  distinct(subject.id)

focus_pos <- focus %>%
  filter(focus.ct > 0)

focus_pos_fine_subjects <- focus_pos %>%
  filter(sample.type == "GII condensate NO mask") %>%
  distinct(subject.id, date.visit) %>%
  distinct(subject.id)
print(nrow(focus_pos_fine_subjects))

focus_pos_person_visits <- focus_pos %>%
  distinct(subject.id, date.visit)
print(nrow(focus_pos_person_visits))

focus_pos_fine <- focus_pos %>%
  filter(sample.type == "GII condensate NO mask")

# (Could also find the number of fine positive where focus was negative (make sure this uses the set where focus was taken))

# Find the number focus positive where fine was negative (both replicates)
focus_pos_fine_ever_pos_person_visits <- focus_pos_fine %>%
  filter(!is.na(final.copies)) # %>%
# distinct(subject.id, date.visit)
# print(nrow(focus_pos_fine_ever_pos_person_visits))

focus_pos_fine_ever_pos_subjects <- focus_pos_fine %>%
  filter(!is.na(final.copies)) %>%
  distinct(subject.id)
print(nrow(focus_pos_fine_ever_pos_subjects))

focus_pos_fine_neg <- focus_pos_fine %>%
  anti_join(focus_pos_fine_ever_pos_subjects) %>%
  distinct(subject.id, date.visit)
print(nrow(focus_pos_fine_neg))

focus_pos_fine_pos <- focus_pos_fine %>%
  anti_join(focus_pos_fine_neg)

focus_pos_fine_pos_person_visits <- focus_pos_fine %>%
  anti_join(focus_pos_fine_neg) %>%
  distinct(subject.id, date.visit)
print(nrow(focus_pos_fine_pos_person_visits))

# Find the set of sample instances where focus was positive and fine was 1 positive and 1 negative
# Use the LOD/sqrt(2) to impute a value for these nondetects. 
# Looks like there are 6 instances where there are NA for final.copies (5 fluA and 1 fluB)
# focus_pos_fine_pos$final.copies[is.na(focus_pos_fine_pos$final.copies)] <- (500/(2)^(.5))

# Now we are changing the plan and will exclude pcr results if one of the replicates was a nondetect
one_or_both_replicates_negative <- focus_pos_fine %>%
  filter(is.na(final.copies)) %>%
  distinct(subject.id, date.visit)

focus_pos_fine_pos <- focus_pos_fine %>%
  anti_join(one_or_both_replicates_negative)

focus_pos_fine_pos_averaged_replicates <- focus_pos_fine_pos %>%
  group_by(subject.id, date.visit) %>%
  mutate(avg.final.copies = mean(final.copies)) %>%
  mutate(avg.focus.ct = mean(focus.ct)) %>%
  mutate(copy_number_ffu_ratio = avg.final.copies / avg.focus.ct) %>%
  mutate(log10_avg.final.copies = log10(avg.final.copies)) %>%
  mutate(log10_avg.focus.ct = log10(avg.focus.ct)) %>%
  mutate(log10_copy_number_log10_ffu_ratio = log10_avg.final.copies / log10_avg.focus.ct) %>%
  mutate(log10_ratio = log10(copy_number_ffu_ratio)) %>%
  distinct(subject.id, date.visit, .keep_all = TRUE)

## Estimate the mean and variance for the ratio of ffu/copy number in fine aerosols.
## Use a mixed model with random effect of person and fixed effect of day since onset of symptoms. 

# Random intercept (to take by-subject variability into account) with fixed mean model 
# (assumes the effect of dpo to be the same across all subjects)
model_1 <- lmer(copy_number_ffu_ratio ~ dpo + (1|subject.id), data = focus_pos_fine_pos_averaged_replicates)
summary(model_1)
coef(model_1)

summary(focus_pos_fine_pos_averaged_replicates$copy_number_ffu_ratio)
# This is quite skewed

model_1_log_scale <- lmer(log10_copy_number_log10_ffu_ratio ~ dpo + (1|subject.id), data = focus_pos_fine_pos_averaged_replicates)
summary(model_1_log_scale)

summary(focus_pos_fine_pos_averaged_replicates$log10_copy_number_log10_ffu_ratio)
# This is less skewed but the team decided to go with a different plan - see below.

# But if we want the fixed effect of dpo to not be continuous, then we will convert this dpo variable to character
# focus_pos_fine_pos_averaged_replicates$dpo <- ordered(focus_pos_fine_pos_averaged_replicates$dpo)
# Could use dpo as ordinal - as done above, but decided to go with using dpo as a character variable for now. 

focus_pos_fine_pos_averaged_replicates$dpo <- as.character(focus_pos_fine_pos_averaged_replicates$dpo)

# Also, in fitting this second model, we want to use the log10(ratio) because the model was somewhat skewed before
model_2 <- lmer(log10_ratio ~ dpo + (1|subject.id), data = focus_pos_fine_pos_averaged_replicates)
summary(model_2)

summary(focus_pos_fine_pos_averaged_replicates$log10_ratio)

table(focus_pos_fine_pos_averaged_replicates$dpo)

# writing this out so that Charles Ma can work on it related to the power calculation for new RO1 grants due Februrary 5, 2019
# write.csv(focus_pos_fine_pos_averaged_replicates, "/Users/jbueno/Desktop/focus_pos_fine_pos.csv")
