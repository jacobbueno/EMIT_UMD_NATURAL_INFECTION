# EMIT_UMD_Natural_Infection_Analysis.R
# Program Objective: Run analysis on EMIT Natural Infection datasets
# Author: Jacob Bueno de Mesquita using material from Jing Yan and Don Milton
# Date: December 14, 2018 -- February 2019, August 2019
# Summary: While most of the analysis for the EMIT_UMD study was done in SAS (especially that as part of the PNAS manuscript), here are some others analyses that we will include as well. Here, we estimate the mean and variance for the ratio of ffu/copy number in fine aerosols from EMIT using a mixed model with random effect of person and fixed effect of day since onset of symptoms. 

#### Load required packages and set working directory ####

library(tidyverse)
library(RcppRoll)
library(readxl)
library(knitr)
library(data.table)
library(lubridate)
library(arsenal)
library(lme4)

# setwd("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection")
# Commented out the setwd() command in order to properly run the markdown report compilation.

sessionInfo() # for reproducibility

#### Estimate the mean and variance for the ratio of ffu/copy number in fine aerosols from EMIT ####
# Use a mixed model with random effect of person and fixed effect of day since onset of symptoms. 
# First have to cut the right df

PNAS_data_full <- read.csv("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Curated Data/Analytical Datasets/PNAS_data_full.csv")
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



## For Gene Tan JCVI Sequencing ####
# What were the ct values of the influenza B virus samples.

# Flu B NPS

flub_nps <- PNAS_data_full %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  mutate(mean_ct = mean(Ct, na.rm = TRUE))


# Flu B fine aerosol

flub_fine <- PNAS_data_full %>%
  filter(type == "B") %>%
  filter(sample.type == "GII condensate NO mask") %>%
  mutate(mean_ct = mean(Ct, na.rm = TRUE))


# Flu B coarse aerosol

flub_coarse <- PNAS_data_full %>%
  filter(type == "B") %>%
  filter(sample.type == "Impactor 5 um NO mask") %>%
  mutate(mean_ct = mean(Ct, na.rm = TRUE))


# Flu A NPS

flua_nps <- PNAS_data_full %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  mutate(mean_ct = mean(Ct, na.rm = TRUE))


# Flu A fine aerosol

flua_fine <- PNAS_data_full %>%
  filter(type == "A") %>%
  filter(sample.type == "GII condensate NO mask") %>%
  mutate(mean_ct = mean(Ct, na.rm = TRUE))


# Flu A coarse aerosol

flua_coarse <- PNAS_data_full %>%
  filter(type == "A") %>%
  filter(sample.type == "Impactor 5 um NO mask") %>%
  mutate(mean_ct = mean(Ct, na.rm = TRUE))


## Checking JCVI sequencing metadata with Todd Treagen - August 2019

## Read in the batch data in the files from JCVI

UMD_batch1 <- read_xls("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Sequence_Data/JCVI/20190306_UMDA_batch1_final_stats.xls", sheet = "UMDA_batch1")

UMD_batch1_draft_status <- read_xls("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Sequence_Data/JCVI/20190306_UMDA_batch1_final_stats.xls", sheet = "Draft Status", col_names = FALSE)

UMD_batch2 <- read_xls("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Sequence_Data/JCVI/20190411_UMDA_batch2_final_stats.xls", sheet = "UMDA_batch2")

UMD_batch2_draft_status <- read_xls("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Sequence_Data/JCVI/20190411_UMDA_batch2_final_stats.xls", sheet = "Draft Status", col_names = FALSE)


## Deal with the formatting in the "Draft Status" sheets -- fill in the JCVI Bac ID for each place where it is missing (in the UMD_batch1_draft_status and UMD_batch2_draft_status dfs)

## For batch1 and batch2, set the column names, remove the blank rows, then fill in the JCVI Bac IDs where needed

# batch1...

UMD_batch1_draft_status_cols <- UMD_batch1_draft_status %>%
  rename(JCVI_Bac_ID = X__1,
         segment = X__2,
         note = X__3) %>%
  filter(!is.na(segment))

for (row in 2:length(UMD_batch1_draft_status_cols$JCVI_Bac_ID)) { # 2 to not affect column names
  if(is.na(UMD_batch1_draft_status_cols$JCVI_Bac_ID[row])) {  
    UMD_batch1_draft_status_cols$JCVI_Bac_ID[row] = UMD_batch1_draft_status_cols$JCVI_Bac_ID[row - 1] 
  }
}


# batch2...


UMD_batch2_draft_status_cols <- UMD_batch2_draft_status %>%
  rename(JCVI_Bac_ID = X__1,
         segment = X__2,
         note = X__3) %>%
  filter(!is.na(segment))

for (row in 2:length(UMD_batch2_draft_status_cols$JCVI_Bac_ID)) { # 2 to not affect column names
  if(is.na(UMD_batch2_draft_status_cols$JCVI_Bac_ID[row])) {  
    UMD_batch2_draft_status_cols$JCVI_Bac_ID[row] = UMD_batch2_draft_status_cols$JCVI_Bac_ID[row - 1] 
  }
}


## Merge batch data with the draft status for batches 1 and 2

UMD_batch1_with_status <- UMD_batch1 %>%
  left_join(UMD_batch1_draft_status_cols, by = c("JCVI Bac ID" = "JCVI_Bac_ID"))

UMD_batch2_with_status <- UMD_batch2 %>%
  left_join(UMD_batch2_draft_status_cols, by = c("JCVI Bac ID" = "JCVI_Bac_ID"))


## Read in the sample metadata produced at UMD

metadata_sample_codes <- read_csv("/Users/jbueno/Box Sync/EMIT/EMIT_Data_Analysis_Jake/EMIT_UMD_Natural_Infection/Sequence_Data/JCVI/EMIT_seq_metadata_simple_data_structure.csv")


## Row bind the batch 1 and batch 2 together and then add the merge in the metadata. 

UMD_seq_batches_1_and_2 <- UMD_batch1_with_status %>%
  bind_rows(UMD_batch2_with_status) %>%
  left_join(metadata_sample_codes, by = c("Blinded Number" = "JCVI Sticker Code"))
  
batch_without_match <- UMD_seq_batches_1_and_2 %>%
  filter(is.na(`Samples Picked for JCVI Sequencing`)) %>%
  distinct(`JCVI Bac ID`, `Blinded Number`, `Organism Name provided by collaborator`, `Updated Organism Name (names that changed are in blue)`, `Special Notes`, `Complete/Draft`)


UMD_metadata_seq_batches_1_and_2 <- UMD_batch1_with_status %>%
  bind_rows(UMD_batch2_with_status) %>%
  right_join(metadata_sample_codes, by = c("Blinded Number" = "JCVI Sticker Code"))
  
metadata_without_match <- UMD_metadata_seq_batches_1_and_2 %>%
  filter(is.na(`JCVI Bac ID`)) %>%
  distinct(`Blinded Number`, date.visit, pcrTest, `Samples Picked for JCVI Sequencing`, `Sample Type`, `Location in UMD Box`)


full_seq_plus_metadata <- UMD_batch1_with_status %>%
  bind_rows(UMD_batch2_with_status) %>%
  full_join(metadata_sample_codes, by = c("Blinded Number" = "JCVI Sticker Code"))
  



## Checking Nancy Leung's 3-Climate Paper data ####

# It appears that Nancy may have used just the first day of G2 sampling from each participant. To test this, I will compute the number and percent of G2 positive samples, the GM and GSD, and the the range for the first day of G2 sampling for each participant, as well as for all days of G2 sampling combined (over all 3 potential study days). 

# We can filter the data to a subset where g2.run == 1 and then compute for NPS, coarse, and fine

PNAS_data_full_g2_run_1 <- PNAS_data_full %>%
  filter(g2.run == 1)

PNAS_data_full_g2_run_1_A_NPS_1_obs <- PNAS_data_full_g2_run_1 %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count == 1) %>%
  select(-count) %>%
  left_join(PNAS_data_full_g2_run_1) %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies)
  

PNAS_data_full_g2_run_1_A_NPS_2_plus_obs <- PNAS_data_full_g2_run_1 %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count >= 2) %>%
  select(-count) %>%
  left_join(PNAS_data_full_g2_run_1) %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies) %>%
  mutate(final.copies = if_else(is.na(final.copies), 2000*1/sqrt(2), final.copies)) %>% 
  group_by(subject.id, dpo) %>%
  mutate(final.copies = mean(final.copies)) %>%
  distinct(subject.id, dpo, final.copies)


NPS_A <- PNAS_data_full_g2_run_1_A_NPS_1_obs %>%
  bind_rows(PNAS_data_full_g2_run_1_A_NPS_2_plus_obs) %>%
  filter(!is.na(final.copies)) %>%
  mutate(ln.final.copies = log(final.copies)) %>%
  ungroup() %>%
  summarize(A_Fine_Positive_Samples_Geom_Mean = exp(mean(ln.final.copies, na.rm = TRUE)),
            A_Fine_Positive_Samples_GSD = exp(sd(ln.final.copies, na.rm = TRUE)),
            A_Fine_Positive_Samples_Max = exp(max(ln.final.copies, na.rm = TRUE)),
            A_Fine_Positive_Samples_n = n())
NPS_A$A_Fine_Positive_Samples_Geom_Mean <- 
  format(NPS_A$A_Fine_Positive_Samples_Geom_Mean, scientific = TRUE)
NPS_A$A_Fine_Positive_Samples_Max <- 
  format(NPS_A$A_Fine_Positive_Samples_Max, scientific = TRUE)

### Repeating the above for flu B using only the first day of g2 runs
PNAS_data_full_g2_run_1 <- PNAS_data_full %>%
  filter(g2.run == 1)

PNAS_data_full_g2_run_1_B_NPS_1_obs <- PNAS_data_full_g2_run_1 %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count == 1) %>%
  select(-count) %>%
  left_join(PNAS_data_full_g2_run_1) %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies)


PNAS_data_full_g2_run_1_B_NPS_2_plus_obs <- PNAS_data_full_g2_run_1 %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count >= 2) %>%
  select(-count) %>%
  left_join(PNAS_data_full_g2_run_1) %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies) %>%
  mutate(final.copies = if_else(is.na(final.copies), 2000*1/sqrt(2), final.copies)) %>% 
  group_by(subject.id, dpo) %>%
  mutate(final.copies = mean(final.copies)) %>%
  distinct(subject.id, dpo, final.copies)


NPS_B <- PNAS_data_full_g2_run_1_B_NPS_1_obs %>%
  bind_rows(PNAS_data_full_g2_run_1_B_NPS_2_plus_obs) %>%
  filter(!is.na(final.copies)) %>%
  mutate(ln.final.copies = log(final.copies)) %>%
  ungroup() %>%
  summarize(B_Fine_Positive_Samples_Geom_Mean = exp(mean(ln.final.copies, na.rm = TRUE)),
            B_Fine_Positive_Samples_GSD = exp(sd(ln.final.copies, na.rm = TRUE)),
            B_Fine_Positive_Samples_Max = exp(max(ln.final.copies, na.rm = TRUE)),
            B_Fine_Positive_Samples_n = n())
NPS_B$B_Fine_Positive_Samples_Geom_Mean <- 
  format(NPS_B$B_Fine_Positive_Samples_Geom_Mean, scientific = TRUE)
NPS_B$B_Fine_Positive_Samples_Max <- 
  format(NPS_B$B_Fine_Positive_Samples_Max, scientific = TRUE)





### Repeating the above but using all the data, not just the first day of g2 sampling
# For flu A
PNAS_data_full_g2_run_1_A_NPS_1_obs <- PNAS_data_full %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count == 1) %>%
  select(-count) %>%
  left_join(PNAS_data_full) %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies)


PNAS_data_full_g2_run_1_A_NPS_2_plus_obs <- PNAS_data_full %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count >= 2) %>%
  select(-count) %>%
  left_join(PNAS_data_full) %>%
  filter(type == "A") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies) %>%
  mutate(final.copies = if_else(is.na(final.copies), 2000*1/sqrt(2), final.copies)) %>% 
  group_by(subject.id, dpo) %>%
  mutate(final.copies = mean(final.copies)) %>%
  distinct(subject.id, dpo, final.copies)


NPS_A <- PNAS_data_full_g2_run_1_A_NPS_1_obs %>%
  bind_rows(PNAS_data_full_g2_run_1_A_NPS_2_plus_obs) %>%
  filter(!is.na(final.copies)) %>%
  mutate(ln.final.copies = log(final.copies)) %>%
  ungroup() %>%
  summarize(A_Fine_Positive_Samples_Geom_Mean = exp(mean(ln.final.copies, na.rm = TRUE)),
            A_Fine_Positive_Samples_GSD = exp(sd(ln.final.copies, na.rm = TRUE)),
            A_Fine_Positive_Samples_Max = exp(max(ln.final.copies, na.rm = TRUE)),
            A_Fine_Positive_Samples_n = n())
NPS_A$A_Fine_Positive_Samples_Geom_Mean <- 
  format(NPS_A$A_Fine_Positive_Samples_Geom_Mean, scientific = TRUE)
NPS_A$A_Fine_Positive_Samples_Max <- 
  format(NPS_A$A_Fine_Positive_Samples_Max, scientific = TRUE)

# For flu B

PNAS_data_full_g2_run_1_B_NPS_1_obs <- PNAS_data_full %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count == 1) %>%
  select(-count) %>%
  left_join(PNAS_data_full) %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies)


PNAS_data_full_g2_run_1_B_NPS_2_plus_obs <- PNAS_data_full %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  group_by(subject.id, dpo) %>%
  mutate(sample_copy_mean = mean(final.copies, na.rm = TRUE)) %>% # doing this because there were some instances where there were repeated NPS pcr assays - taking the mean of these repeats here. %>%
  distinct(subject.id, dpo, final.copies) %>%
  summarize(count = n()) %>%
  filter(count >= 2) %>%
  select(-count) %>%
  left_join(PNAS_data_full) %>%
  filter(type == "B") %>%
  filter(sample.type == "Nasopharyngeal swab") %>%
  select(subject.id, dpo, final.copies) %>%
  mutate(final.copies = if_else(is.na(final.copies), 2000*1/sqrt(2), final.copies)) %>% 
  group_by(subject.id, dpo) %>%
  mutate(final.copies = mean(final.copies)) %>%
  distinct(subject.id, dpo, final.copies)


NPS_B <- PNAS_data_full_g2_run_1_B_NPS_1_obs %>%
  bind_rows(PNAS_data_full_g2_run_1_B_NPS_2_plus_obs) %>%
  filter(!is.na(final.copies)) %>%
  mutate(ln.final.copies = log(final.copies)) %>%
  ungroup() %>%
  summarize(B_Fine_Positive_Samples_Geom_Mean = exp(mean(ln.final.copies, na.rm = TRUE)),
            B_Fine_Positive_Samples_GSD = exp(sd(ln.final.copies, na.rm = TRUE)),
            B_Fine_Positive_Samples_Max = exp(max(ln.final.copies, na.rm = TRUE)),
            B_Fine_Positive_Samples_n = n())
NPS_B$B_Fine_Positive_Samples_Geom_Mean <- 
  format(NPS_B$B_Fine_Positive_Samples_Geom_Mean, scientific = TRUE)
NPS_B$B_Fine_Positive_Samples_Max <- 
  format(NPS_B$B_Fine_Positive_Samples_Max, scientific = TRUE)









